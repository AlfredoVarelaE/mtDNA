# FileName: human_SNVs

## Strategy for the identification of polymorphic single nucleotide variants in the human mitochondrial genome

This document is supplementary to the article:
Alfredo Varela-Echavarría, Kenya L. Contreras-Ramírez, Carlos Lozano-Flores, and Maribel Hernández-Rosales (2024) Detection of single nucleotide variants in the mitochondrial genome of healthy mice and humans.

(citation to be inserted here)

To map precisely the reads corresponding to the circular junction regions at the D-loop, tail-to-head junction sequences were made for the human mitochondrial genome by joining the 3kb of the end of the rCRS reference genome to 3kb of the start of the same sequence (THjnct)(Figure below). All reads from each sample were then aligned to the THjnct sequence and the full length linear sequence (FL).



 This code parses and filters mutect2 output of FL mapping, selects single nucleotide events, inverts alleles if necessary to have major allele (MaA) counts in field 5 and minor allele (MiA) counts in field 6 adding the “i” sufix to both alleles in REF and ALT if inversion is made, eliminates strand bias both in MaA and MiA keeping only those in which strand ratios fall within the range 0.66-1.5, selects events in which MiA depth is at least 10 with at least five events on each strand, eliminates 1kb from each tail and head ends, and calculates allele frequency (MiA/MaA). Oneliners are separate so that progress can be verified at each step during test runs before batch runs are made.

 ```bash
zcat seqfile_FL.vcf.gz | \
grep ^chrM | \
grep  -e '0/1:' -e '0|1\:'| \
awk  '{print $2"\t"$4"\t"$5"\t"$10}' | \
awk -F ':' '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}'|\
awk -F ',' '{print $1"\t"$2"\t"$3"\t"$4}' | \
awk '($7>=0.8) && ($1>1000) && ($1<15570){print $0}'|sort -k1n,1 > seqfile_FL.HAP.fltrd
```

Filtering mutect2 output of tail-to-head junction (THjnct) mapping:
This code parses and selects variants from mutect2 output of THjnct mapping as described for FL mapping. It also selects events that lie within the 1 kb at the tail and 1kb at the head of the genome also correcting the position to the full length genome.

```bash
zcat seqfile_TH.vcf.gz | \
grep ^chrM | \
grep  -e '0/1:' -e '0|1\:'| \
awk  '{print $2"\t"$4"\t"$5"\t"$10}' | \
awk -F ':' '{print $1"\t"$2"\t"$3"\t"$5"\t"$6}'|\
awk -F ',' '{print $1"\t"$2"\t"$3"\t"$4}' | \
awk '{if ($7>=0.8 && $1>2000 && $1<3001) print $1+13569"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11; else if ($7>=0.8 && $1>3000 && $1<4001) print $1-3000"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}'> seqfile_TH.HAP.fltrd.1
```
Both, FL and THjnct events are then merged into a single variant file for each sample.
