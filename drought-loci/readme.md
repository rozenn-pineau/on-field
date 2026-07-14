Goal: to intersect the sync pooled sequencing file with the minimum filtering with the drought loci bed file. 

### Copy file from storage on cds to scratch
```
cp /cds3/kreiner/2024_onfield_pools/sync/excluded_mid_filtered.sync ./
#takes a while because the file is big (>600GB for min filtered file)
```

### Prepping the file for bedtools intersect

(1) Remove Scaffold_ from chromosome name (bedtools wants integer)

(2) Make bed from sync file for bedtools intersect

(3) Split into scaffolds (memory intensive job)

```
# --> Remove Scaffold_ from chromosome name (bedtools wants integer)
cat excluded_mid_filtered.sync | \sed s/^\Scaffold_//g > excluded_mid_filtered_integer.sync
```
Following steps ran as a batch job:

```
#!/bin/bash
#SBATCH --job-name=split
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10G

# --> Prep bed format for syn file and intersect
cd /scratch/midway2/rozennpineau/on-field/sync
awk -F'\t' -v OFS='\t' 'BEGIN{OFS=FS} {$3 = $2 OFS $3} 1' excluded_mid_filtered_integer.sync > excluded_mid_filtered_integer.sync.bed

# --> Split into Scaffolds (bedtools intersect command line out of memory if not)
cd /scratch/midway2/rozennpineau/on-field/sync/scaffolds
awk -F'\t' '{print $0 > ($1 ".sync")}' ../excluded_mid_filtered_integer.sync.bed
```


### Run bedtools intersect

```
#ran as a batch job

#conda environment
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner/rpineau/bedtools

bed=/scratch/midway2/rozennpineau/on-field/drought_bed/ancestry_corrected_inflated_gemma_gwas_893_header_numeric.assoc.bed
out=/scratch/midway2/rozennpineau/on-field/intersect
cd $out

#intersect
#bedtools version v2.31.1 on March 31
cd /scratch/midway2/rozennpineau/on-field/sync/scaffolds

for sync in *.sync; do
        echo $sync
        bedtools intersect -a $bed -b $sync -wb -f 0.99 -r >> ../../intersect/intersect_sync_drought_stringent.sync
        #-r 1 -f 1 requires that there is at least 1 bp match, with 100% of the regions matching
        # -wb : keep the overlap from b file
        bedtools intersect -a $bed -b $sync -wb >> ../../intersect/intersect_sync_drought.sync
done

```

The more stringent bedtools command is the one to keep, as the less stringent one kept the loci before and after the focal SNP. 

### Clean intersection file
```
#keep sync file part only
cut -f 8- intersect_sync_drought_stringent.sync > intersect_sync_drought_stringent_noheader.sync
#add header
cat ../sync/header intersect_sync_drought_stringent_noheader.sync > intersect_drought_min.sync
```

### From sync file to AF matrix

Shortcuts/assumptions of this script: 
Only biallelic (not tri-allelic) SNP are being considered. 
Anything with less than 5 reads is set to NA. 

```
#Rscript
#upload sync file
sync <- read.table("/scratch/midway3/rozennpineau/on-field/drought-loci/intersect/intersect_sync_drought_stringent_noheader.sync", sep = "\t", header = F)

ncol <- dim(sync)[2]
nloci <- dim(sync)[1]
start_col <- 4
npool <- ncol - start_col

lookup_table <- c("A" = 1, "T" = 2, "C" = 3, "G" = 4, "N" = 5)

#prepare output file
AF_mat <- matrix(NA, nloci, npool) #keep track of allele frequencies for each pool

#Loop through loci
for (locus in 1:nloci) {
        #locus <- 1
        alt <- sync[locus,4]
        #Loop through pools
        for (pool in 1:npool) {
                #pool <- 1
                cur_pool <- as.numeric(unlist(strsplit(sync[locus, start_col + pool],":")))
                #calculate sum
                total_num_reads <- sum(cur_pool)
                #if total reads < 5 --> NA
                if (total_num_reads <= 5 | is.na(total_num_reads)) { AF_mat[locus, pool] <- NA}
                #calculate AF based on alternate allele
                else if (total_num_reads > 5) {AF_mat[locus, pool] <- cur_pool[lookup_table[alt]]/total_num_reads}
                }
        }

#add chrom, pos, pos, ALT to the AF matrix before writing out
AF_mat_export <- cbind(sync[,1:4], AF_mat)
colnames(AF_mat_export) <- c("chrom", "pos", "pos", "alt", 1:npool)
write.table(AF_mat_export, "drought_loci_AF_matrix.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
```




