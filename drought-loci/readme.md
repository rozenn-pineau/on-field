Goal: to intersect the sync pooled sequencing file with the minimum filtering with the drought loci bed file. 

### Prepping the file for bedtools intersect
```
# --> Remove Scaffold_ from chromosome name (bedtools wants integer)
cat excluded_min_filtered.sync | \sed s/^\Scaffold_//g > excluded_min_filtered_integer.sync
# --> Prep bed format for syn file and intersect
awk -F'\t' -v OFS='\t' 'BEGIN{OFS=FS} {$3 = $2 OFS $3} 1' excluded_min_filtered_integer.sync > excluded_min_filtered_integer.sync.bed
#I had to run this step as a batch job on the cluster
# --> Split into Scaffolds (bedtools intersect command line out of memory if not)
awk -F'\t' '{print $0 > ($1 ".sync")}' ../excluded_min_filtered_integer.sync.bed
#also ran on the cluster as a batch job
```

