# Octoploid examples

The synthetic assembly contains 16 chromosomes in two homolog groups of eight.
Run these commands from the repository root.

## Two original sources

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --sources examples/octoploid/two_sources.tsv \
  --homolog-groups examples/octoploid/groups.tsv --ploidy 8 \
  --source-copies hap1=4 --source-copies hap2=4 \
  --min-aln-len 1 --out-dir octoploid_two_source_output
```

Expected: `hap1.fa` and `hap2.fa`, each with eight chromosomes; four from each
homolog group. All group checks pass.

## Eight known sources

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --sources examples/octoploid/sources.tsv \
  --homolog-groups examples/octoploid/groups.tsv --ploidy 8 \
  --min-aln-len 1 --out-dir octoploid_output
```

Expected: `copy01.fa` through `copy08.fa`, each with two chromosomes; one from
each homolog group. All group checks pass.
