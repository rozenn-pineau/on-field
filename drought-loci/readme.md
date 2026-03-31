Goal: to intersect the sync pooled sequencing file with the minimum filtering with the drought loci bed file. 

### Prepping the file for bedtoold intersect
--> Remove Scaffold_ from chromosome name (bedtools wants integer)
```
cat excluded_min_filtered.sync | \sed s/^\Scaffold_//g > excluded_min_filtered_integer.sync
```
### Prep bed format for syn file and intersect
```
#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --time=5:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-kreiner
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10GB

#conda environment
module load python/anaconda-2022.05
source /software/python-anaconda-2022.05-el8-x86_64/etc/profile.d/conda.sh
conda activate /project/kreiner

# Print pos1 pos2 and rest of file with tab delimitations in sync file
cd /scratch/midway3/rozennpineau/on-field/drought-loci/sync
awk -F'\t' -v OFS='\t' '{$3=$2} 1' excluded_min_filtered_integer.sync > excluded_min_filtered_integer.sync.bed

#load files
cd /scratch/midway3/rozennpineau/on-field/drought-loci/intersect
#before clumping set of drought loci (893 loci with FDR < 0.05)
bed=/scratch/midway3/rozennpineau/on-field/drought-loci/beds/ancestry_corrected_inflated_gemma_gwas_893_header.assoc.bed
#sync file with the minimim filtering approach as we are looking for specific loci
sync=/scratch/midway3/rozennpineau/on-field/drought-loci/sync/excluded_min_filtered_integer.sync.bed

#intersect
#bedtools version v2.31.1 on March 31
bedtools intersect -a $bed -b $sync -wb -f 0.99 -r > intersect_sync_drought_stringent.sync
#-r 1 -f 1 requires that there is at least 1 bp match, with 100% of the regions matching
# -wb : keep the overlap from b file
bedtools intersect -a $bed -b $sync -wb > intersect_sync_drought.bed
```
