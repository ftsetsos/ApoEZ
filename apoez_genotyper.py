#!/usr/bin/env python3
"""Call APOE epsilon genotypes from rs429358 and rs7412 in a VCF.

GRCh38 loci:
  rs429358  chr19:44908684 T>C
  rs7412    chr19:44908822 C>T

The script uses phased haplotypes when both records contain phased diploid GTs.
For unphased calls, the rs429358 C/T + rs7412 C/T combination is reported as
epsilon1/epsilon3 OR epsilon2/epsilon4 instead of being guessed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path


LOCI = {
    "rs429358": {"chrom": "19", "pos": 44908684, "ref": "T", "alts": {"C"}},
    "rs7412": {"chrom": "19", "pos": 44908822, "ref": "C", "alts": {"T"}},
}

# Only loci with retained genotype-specific magnitude scores. The two APOE
# defining SNPs remain inputs to the epsilon diplotype and are not scored again
# individually, which avoids double-counting APOE in the overall magnitude.
GENOTYPE_SCORES = {
    "rs449647": {
        "pos": 44905307,
        "scores": {"AA": 2.0, "AT": 1.2, "TT": 1.5},
    },
    "rs267606664": {
        "pos": 44908730,
        "scores": {"AG": 4.0, "GG": 0.0},
    },
    "rs121918393": {
        "pos": 44908756,
        "scores": {"AA": 6.0, "AC": 5.0, "CC": 0.0},
        "name": "APOE Christchurch",
    },
    "rs387906567": {
        "pos": 44908774,
        "scores": {"CC": 0.0, "CT": 3.0},
    },
    "rs121918394": {
        "pos": 44908786,
        "scores": {"AA": 0.0, "AG": 5.0},
    },
    "rs199768005": {
        "pos": 44909057,
        "scores": {"AT": 2.1},
    },
    "rs4420638": {
        "pos": 44919689,
        "scores": {"AA": 0.0, "AG": 2.0, "GG": 3.0},
        "name": "APOC1 marker",
    },
}

MARKERS = {
    **{name: {"pos": spec["pos"]} for name, spec in LOCI.items()},
    **{name: {"pos": spec["pos"]} for name, spec in GENOTYPE_SCORES.items()},
}

HAPLOTYPE_TO_EPSILON = {
    ("C", "T"): "ε1",
    ("T", "T"): "ε2",
    ("T", "C"): "ε3",
    ("C", "C"): "ε4",
}

# SNPedia APOE page magnitude values. These are importance ratings, not a
# calibrated probability, odds ratio, or clinical risk score.
MAGNITUDE = {
    "ε1/ε1": 6.0,
    "ε1/ε2": 2.5,
    "ε1/ε3": 2.6,
    "ε1/ε4": 2.5,
    "ε2/ε2": 4.0,
    "ε2/ε3": 2.0,
    "ε2/ε4": 2.6,
    "ε3/ε3": 2.0,
    "ε3/ε4": 3.0,
    "ε4/ε4": 6.0,
    "ε1/ε3 OR ε2/ε4": 2.6,
}

EPSILON_ORDER = {"ε1": 1, "ε2": 2, "ε3": 3, "ε4": 4}


def open_text(path: str):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def normalize_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def locate_record(chrom: str, pos: int, ids: str) -> str | None:
    id_set = set(ids.split(";")) if ids not in {"", "."} else set()
    by_id = [marker for marker in MARKERS if marker in id_set]
    by_pos = [
        marker
        for marker, spec in MARKERS.items()
        if normalize_chrom(chrom) == "19" and spec["pos"] is not None and pos == spec["pos"]
    ]
    hits = set(by_id) | set(by_pos)
    if len(hits) > 1:
        raise ValueError(f"Record {chrom}:{pos} matches multiple target loci")
    return next(iter(hits), None)


def parse_gt(gt: str, alleles: list[str]) -> tuple[list[str] | None, bool]:
    phased = "|" in gt
    sep = "|" if phased else "/"
    fields = gt.split(sep)
    if len(fields) != 2 or any(x == "." for x in fields):
        return None, phased
    try:
        indices = [int(x) for x in fields]
        bases = [alleles[i].upper() for i in indices]
    except (ValueError, IndexError):
        return None, phased
    return bases, phased


def read_targets(vcf: str):
    samples: list[str] | None = None
    records: dict[str, dict] = {}

    with open_text(vcf) as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            if samples is None:
                raise ValueError("VCF is missing a #CHROM header")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"Malformed VCF record at line {line_number}")
            chrom, pos_text, ids, ref, alt_text = fields[:5]
            rsid = locate_record(chrom, int(pos_text), ids)
            if rsid is None:
                continue
            if rsid in records:
                raise ValueError(f"Duplicate target record for {rsid}")

            spec = LOCI.get(rsid)
            alts = alt_text.upper().split(",")
            if spec is not None and (ref.upper() != spec["ref"] or not spec["alts"].intersection(alts)):
                raise ValueError(
                    f"Unexpected alleles for {rsid}: {ref}>{alt_text}; "
                    f"expected GRCh38 {spec['ref']}>{'/'.join(sorted(spec['alts']))}"
                )
            format_keys = fields[8].split(":")
            if "GT" not in format_keys:
                raise ValueError(f"{rsid} has no GT FORMAT field")
            gt_index = format_keys.index("GT")
            ps_index = format_keys.index("PS") if "PS" in format_keys else None
            calls = {}
            for sample, sample_field in zip(samples, fields[9:]):
                values = sample_field.split(":")
                gt = values[gt_index] if gt_index < len(values) else "."
                bases, phased = parse_gt(gt, [ref.upper(), *alts])
                ps = values[ps_index] if ps_index is not None and ps_index < len(values) else None
                if ps in {None, "", "."}:
                    ps = None
                calls[sample] = {"bases": bases, "phased": phased, "ps": ps, "raw_gt": gt}
            records[rsid] = {
                "chrom": chrom,
                "pos": int(pos_text),
                "vcf_id": ids,
                "ref": ref.upper(),
                "alt": ",".join(alts),
                "qual": fields[5],
                "filter": fields[6],
                "info": fields[7],
                "calls": calls,
            }

    if samples is None:
        raise ValueError("VCF is missing a #CHROM header")
    missing = sorted(set(LOCI) - set(records))
    if missing:
        raise ValueError("VCF does not contain target record(s): " + ", ".join(missing))
    return samples, records


def canonical_pair(a: str, b: str) -> str:
    x, y = sorted((a, b), key=EPSILON_ORDER.get)
    return f"{x}/{y}"


def genotype_call(call_358: dict, call_7412: dict) -> dict:
    a = call_358["bases"]
    b = call_7412["bases"]
    base = {
        "rs429358_GT": call_358["raw_gt"],
        "rs429358_bases": "." if a is None else "/".join(a),
        "rs7412_GT": call_7412["raw_gt"],
        "rs7412_bases": "." if b is None else "/".join(b),
    }
    if a is None or b is None:
        return {**base, "phase_status": "missing", "APOE_genotype": "UNKNOWN",
                "call_status": "MISSING_GT", "magnitude": ".", "ε1_dose": ".",
                "ε2_dose": ".", "ε3_dose": ".", "ε4_dose": ".",
                "note": "One or both defining genotypes are missing or unsupported"}

    phase_usable = call_358["phased"] and call_7412["phased"]
    ps1, ps2 = call_358["ps"], call_7412["ps"]
    if phase_usable and ps1 is not None and ps2 is not None and ps1 != ps2:
        phase_usable = False

    if phase_usable:
        eps = [HAPLOTYPE_TO_EPSILON[(a[i], b[i])] for i in (0, 1)]
        genotype = canonical_pair(*eps)
        phase_status = "phased_same_PS" if ps1 is not None and ps2 is not None else "phased_no_PS"
        status = "CALLED"
        note = "" if phase_status == "phased_same_PS" else "Phased GT used; no shared PS value available"
    else:
        g358, g7412 = tuple(sorted(a)), tuple(sorted(b))
        if g358 == ("C", "T") and g7412 == ("C", "T"):
            genotype = "ε1/ε3 OR ε2/ε4"
            eps = []
            status = "AMBIGUOUS_PHASE"
            note = "Unphased double heterozygote; phase is required to distinguish the two diplotypes"
        else:
            lookup = {
                (("C", "C"), ("T", "T")): "ε1/ε1",
                (("C", "T"), ("T", "T")): "ε1/ε2",
                (("C", "C"), ("C", "T")): "ε1/ε4",
                (("T", "T"), ("T", "T")): "ε2/ε2",
                (("T", "T"), ("C", "T")): "ε2/ε3",
                (("T", "T"), ("C", "C")): "ε3/ε3",
                (("C", "T"), ("C", "C")): "ε3/ε4",
                (("C", "C"), ("C", "C")): "ε4/ε4",
            }
            genotype = lookup.get((g358, g7412), "UNKNOWN")
            eps = genotype.split("/") if genotype != "UNKNOWN" else []
            status = "CALLED" if genotype != "UNKNOWN" else "INVALID_COMBINATION"
            note = "" if status == "CALLED" else "Allele combination is not in the APOE lookup table"
        phase_status = "unphased" if not (call_358["phased"] and call_7412["phased"]) else "different_PS"

    doses = Counter(eps)
    return {
        **base,
        "phase_status": phase_status,
        "APOE_genotype": genotype,
        "call_status": status,
        "magnitude": MAGNITUDE.get(genotype, "."),
        "ε1_dose": doses.get("ε1", 0) if eps else ".",
        "ε2_dose": doses.get("ε2", 0) if eps else ".",
        "ε3_dose": doses.get("ε3", 0) if eps else ".",
        "ε4_dose": doses.get("ε4", 0) if eps else ".",
        "note": note,
    }


def gt_class_and_alt_dose(raw_gt: str) -> tuple[str, str]:
    fields = raw_gt.replace("|", "/").split("/")
    if len(fields) != 2 or any(x == "." for x in fields):
        return "MISSING", "."
    try:
        alleles = [int(x) for x in fields]
    except ValueError:
        return "INVALID", "."
    dose = sum(x != 0 for x in alleles)
    if alleles == [0, 0]:
        return "HOM_REF", "0"
    if alleles[0] == alleles[1]:
        return "HOM_ALT", str(dose)
    return "HET", str(dose)


def auxiliary_variant_summary(sample: str, records: dict) -> dict:
    """Score the retained non-epsilon markers by exact nucleotide genotype."""
    partial_score = 0.0
    contributors = []
    unresolved = []
    marker_results = {}

    for marker, metadata in GENOTYPE_SCORES.items():
        record = records.get(marker)
        if record is None:
            marker_results[marker] = {"GT": ".", "bases": ".", "magnitude": 0.0}
            unresolved.append(f"{marker}:RECORD_ABSENT_ASSUMED_ZERO")
            continue
        call = record["calls"][sample]
        if call["bases"] is None:
            marker_results[marker] = {"GT": call["raw_gt"], "bases": ".", "magnitude": 0.0}
            unresolved.append(f"{marker}:MISSING_GT_ASSUMED_ZERO")
            continue
        genotype = "".join(sorted(call["bases"]))
        magnitude = metadata["scores"].get(genotype)
        if magnitude is None:
            marker_results[marker] = {
                "GT": call["raw_gt"], "bases": genotype, "magnitude": 0.0
            }
            unresolved.append(f"{marker}:{genotype}:UNSCORED_GENOTYPE_ASSUMED_ZERO")
            continue
        marker_results[marker] = {
            "GT": call["raw_gt"], "bases": genotype, "magnitude": magnitude
        }
        partial_score += magnitude
        if magnitude > 0:
            contributors.append(f"{marker}:{genotype}:mag={magnitude:g}")

    ch = marker_results["rs121918393"]
    if "rs121918393" not in records:
        ch_carrier = "NO_VARIANT_RECORD"
        ch_a_dose = 0
    elif ch["bases"] == ".":
        ch_carrier = "UNKNOWN"
        ch_a_dose = "."
    else:
        ch_a_dose = ch["bases"].count("A")
        ch_carrier = "YES" if ch_a_dose > 0 else "NO"

    result = {
        "extra_magnitude_partial_score": round(partial_score, 6),
        "extra_magnitude_score": round(partial_score, 6),
        "magnitude_score_status": "COMPLETE" if not unresolved else "INCOMPLETE_ASSUMED_ZERO",
        "magnitude_score_unresolved": ";".join(unresolved) if unresolved else ".",
        "extra_magnitude_contributors": ";".join(contributors) if contributors else ".",
        "APOE_Christchurch_GT": ch["GT"],
        "APOE_Christchurch_genotype": ch["bases"],
        "APOE_Christchurch_A_dose": ch_a_dose,
        "APOE_Christchurch_carrier": ch_carrier,
        "APOE_Christchurch_magnitude": ch["magnitude"],
    }
    for marker, values in marker_results.items():
        result[f"{marker}_GT"] = values["GT"]
        result[f"{marker}_genotype"] = values["bases"]
        result[f"{marker}_magnitude"] = values["magnitude"]
    return result


def write_variant_table(path: str, samples: list[str], records: dict) -> None:
    fields = [
        "sample", "marker", "expected_GRCh38_pos", "genotype_score_map",
        "record_present", "CHROM", "POS", "VCF_ID", "REF", "ALT",
        "QUAL", "FILTER", "GT", "genotype_bases", "genotype_magnitude", "zygosity",
        "non_reference_dose", "INFO",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for marker, metadata in MARKERS.items():
            record = records.get(marker)
            for sample in samples:
                if record is None:
                    writer.writerow({
                        "sample": sample, "marker": marker,
                        "expected_GRCh38_pos": metadata["pos"] if metadata["pos"] is not None else ".",
                        "genotype_score_map": "." if marker in LOCI else json.dumps(GENOTYPE_SCORES[marker]["scores"], sort_keys=True),
                        "record_present": "NO", "CHROM": "chr19", "POS": ".", "VCF_ID": ".",
                        "REF": ".", "ALT": ".", "QUAL": ".", "FILTER": ".", "GT": ".",
                        "genotype_bases": ".", "genotype_magnitude": 0.0, "zygosity": "RECORD_ABSENT",
                        "non_reference_dose": ".", "INFO": ".",
                    })
                    continue
                call = record["calls"][sample]
                zygosity, dose = gt_class_and_alt_dose(call["raw_gt"])
                genotype = "." if call["bases"] is None else "".join(sorted(call["bases"]))
                genotype_magnitude = (
                    "." if marker in LOCI
                    else GENOTYPE_SCORES[marker]["scores"].get(genotype, 0.0)
                )
                writer.writerow({
                    "sample": sample, "marker": marker,
                    "expected_GRCh38_pos": metadata["pos"] if metadata["pos"] is not None else ".",
                    "genotype_score_map": "." if marker in LOCI else json.dumps(GENOTYPE_SCORES[marker]["scores"], sort_keys=True),
                    "record_present": "YES", "CHROM": record["chrom"], "POS": record["pos"],
                    "VCF_ID": record["vcf_id"], "REF": record["ref"], "ALT": record["alt"],
                    "QUAL": record["qual"], "FILTER": record["filter"], "GT": call["raw_gt"],
                    "genotype_bases": genotype, "genotype_magnitude": genotype_magnitude,
                    "zygosity": zygosity, "non_reference_dose": dose, "INFO": record["info"],
                })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="Plain-text or gzip/BGZF VCF containing the APOE region")
    parser.add_argument("-o", "--output", required=True, help="Per-sample output TSV")
    parser.add_argument("--summary", help="Optional JSON summary output")
    parser.add_argument(
        "--variants-output",
        help="Optional long-format TSV containing every listed APOE-region marker for every sample",
    )
    args = parser.parse_args()

    try:
        samples, records = read_targets(args.vcf)
        rows = []
        for sample in samples:
            epsilon = genotype_call(
                records["rs429358"]["calls"][sample],
                records["rs7412"]["calls"][sample],
            )
            auxiliary = auxiliary_variant_summary(sample, records)
            epsilon_magnitude = epsilon["magnitude"]
            overall_partial = (
                round(epsilon_magnitude + auxiliary["extra_magnitude_partial_score"], 6)
                if isinstance(epsilon_magnitude, (int, float)) else "."
            )
            overall = (
                round(epsilon_magnitude + auxiliary["extra_magnitude_score"], 6)
                if isinstance(epsilon_magnitude, (int, float))
                else "."
            )
            rows.append({
                "sample": sample,
                **epsilon,
                **auxiliary,
                "overall_magnitude_partial_score": overall_partial,
                "overall_magnitude_score": overall,
            })

        fieldnames = list(rows[0]) if rows else ["sample"]
        with open(args.output, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        if args.variants_output:
            write_variant_table(args.variants_output, samples, records)

        if args.summary:
            summary = {
                "input_vcf": str(Path(args.vcf)),
                "n_samples": len(rows),
                "call_status_counts": dict(Counter(row["call_status"] for row in rows)),
                "APOE_genotype_counts": dict(Counter(row["APOE_genotype"] for row in rows)),
                "requested_marker_count": len(MARKERS),
                "markers_present": [marker for marker in MARKERS if marker in records],
                "markers_absent": [marker for marker in MARKERS if marker not in records],
            }
            with open(args.summary, "w", encoding="utf-8") as handle:
                json.dump(summary, handle, indent=2, ensure_ascii=False)
                handle.write("\n")
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
