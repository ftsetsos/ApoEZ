# ApoEZ
ApoE allele genotyper and risk scoring from WGS VCF data

This standalone Python 3 script reads a plain or gzip/BGZF-compressed GRCh38
VCF slice and calls APOE epsilon diplotypes from:

- `rs429358`, chr19:44908684, T>C
- `rs7412`, chr19:44908822, C>T

No third-party Python packages are required.

## Usage

```bash
python3 apoe_genotyper.py input.apoe.vcf.gz \
  --output apoe_genotypes.tsv \
  --variants-output apoe_region_variants.tsv \
  --summary apoe_summary.json
```

The two records may be identified by rsID or by their exact GRCh38 positions.
Both defining variants must be present. The script validates their REF/ALT
alleles before making calls.

The optional `--variants-output` file contains the two epsilon-defining SNPs
and the seven retained scored markers. It is written in long format with one
row per sample/marker. Each row includes the nucleotide genotype, its exact
genotype-specific magnitude when defined, zygosity, non-reference dosage,
REF/ALT, FILTER, QUAL, INFO, and expected GRCh38 position. Markers absent from
the VCF are retained as `record_present=NO`, so absent records are not confused
with homozygous-reference calls.

## Calling logic

| rs429358 haplotype base | rs7412 haplotype base | APOE allele |
|---|---|---|
| C | T | ε1 |
| T | T | ε2 |
| T | C | ε3 |
| C | C | ε4 |

If both GT fields are phased, the corresponding haplotypes are paired. If
both variants are unphased heterozygotes (`C/T` and `C/T`), the result is
reported as `ε1/ε3 OR ε2/ε4` with `AMBIGUOUS_PHASE`; the script does not guess
ε2/ε4 merely because ε1 is rare. When both records have PS values, they must
match for the phase to be used.

`magnitude` reproduces the genotype magnitude shown on SNPedia's APOE page.
It is an editorial importance rating, not a probability, odds ratio, polygenic
risk score, or clinically validated risk measure. For association analyses,
the allele-dose columns—especially `ε4_dose`—are usually more interpretable.

## Quick test

```bash
python3 apoe_genotyper.py example.vcf \
  -o example.calls.tsv \
  --variants-output example.variants.tsv \
  --summary example.summary.json
```

Expected special cases include a phased ε2/ε4 call, an unphased ambiguous
double heterozygote, and a missing call.

## Output columns

- `sample`
- raw GT and decoded bases for both SNPs
- `phase_status`
- `APOE_genotype`
- `call_status`
- SNPedia `magnitude`
- `ε1_dose`, `ε2_dose`, `ε3_dose`, `ε4_dose`
- `note`

The optional variant table is deliberately separate from the epsilon genotype
table. Rare APOE-region variants do not redefine the conventional epsilon
diplotype.

## Per-person auxiliary-variant score

The main `apoe_genotypes.tsv` assigns exact genotype-specific magnitudes to:

| Marker | Genotype magnitudes |
|---|---|
| rs449647 | AA=2; AT=1.2; TT=1.5 |
| rs267606664 | AG=4; GG=0 |
| rs121918393 | AA=6; AC=5; CC=0 |
| rs387906567 | CC=0; CT=3 |
| rs121918394 | AA=0; AG=5 |
| rs199768005 | AT=2.1 |
| rs4420638 | AA=0; AG=2; GG=3 |

The conventional APOE genotype magnitude is calculated separately from
`rs429358` and `rs7412`. Their individual SNP magnitudes are not added again,
avoiding double-counting APOE.

- `extra_magnitude_score` is the sum of the available exact genotype
  magnitudes. Auxiliary records that are absent, have a missing GT, or have an
  unlisted genotype contribute zero.
- `overall_magnitude_score` is the APOE ε-genotype magnitude plus
  `extra_magnitude_score`.
- `extra_magnitude_partial_score` and `overall_magnitude_partial_score` are
  retained for compatibility and have the same numeric values under this
  zero-contribution policy.
- `magnitude_score_status` is `COMPLETE` when all seven extra genotypes are
  available and scored; otherwise it is `INCOMPLETE_ASSUMED_ZERO`.
- `magnitude_score_unresolved` lists every absent, missing, or unlisted call
  that was assigned zero.
- `extra_magnitude_contributors` lists all nonzero components.

Each marker also receives `_GT`, `_genotype`, and `_magnitude` columns in the
main output. These sums are descriptive prioritization indices, not validated
clinical or biological risk scores; protective and adverse variants are not
assigned opposite signs.

The Christchurch variant is reported explicitly as:

- `APOE_Christchurch_GT`
- `APOE_Christchurch_genotype`
- `APOE_Christchurch_A_dose`
- `APOE_Christchurch_carrier`
- `APOE_Christchurch_magnitude`

For GRCh38, Christchurch is chr19:44908756 C>A, APOE p.Arg154Ser (historically
called R136S). Its genotype magnitude is CC=0, AC=5, or AA=6. The carrier flag
checks the called nucleotide genotype for the A allele. If the Christchurch
record is absent from the VCF, it contributes zero and the carrier field is
`NO_VARIANT_RECORD`; it does not invalidate the overall magnitude.

This is intended for research workflows. Any clinical or participant-facing
interpretation needs an appropriately validated assay and review process.

## Cohort plots and summaries

Use `plot_apoe_results.py` on the per-sample genotype output:

```bash
python3 plot_apoe_results.py apoe_genotypes.tsv \
  --output-dir apoe_plots \
  --prefix genesis_apoe
```

Requirements are Python 3 with `pandas`, `numpy`, and `matplotlib`. The plotting
script uses a noninteractive backend and can run directly in an LSF job.

It generates eight figures as both PNG and SVG:

1. Four-panel cohort dashboard
2. APOE genotype counts and percentages
3. ε1–ε4 allele frequencies
4. Overall-magnitude histogram
5. Overall magnitude stratified by APOE genotype
6. Genotype composition at the seven auxiliary markers
7. Nonzero prevalence and cumulative magnitude by marker
8. Score completeness and assumed-zero QC

The same figures are collected into `<prefix>.cohort_overview.pdf`. Additional
TSVs report cohort metrics, APOE genotype counts, allele frequencies,
marker-genotype distributions, marker contributions, and score-status counts.
