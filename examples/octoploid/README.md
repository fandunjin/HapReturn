# Two-source example

The synthetic assembly contains 16 chromosome records. Use the ordinary hap1/hap2
command from the repository root:

```bash
hapreturn \
  --final-fasta examples/octoploid/final.fa \
  --hap1-paf examples/octoploid/hap1.paf \
  --hap2-paf examples/octoploid/hap2.paf \
  --ploidy 8 --min-aln-len 1 --out-dir example_result
```

Expected: `hap1.fa` and `hap2.fa`, each with eight chromosomes; zero unresolved
records. Records are ordered as two consecutive groups of eight homologs.
Other even ploidies use the same command with their corresponding --ploidy.
