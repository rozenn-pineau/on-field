#Rscript
#upload sync file
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/drought_loci")
sync <- read.table("intersect_sync_drought_stringent.sync", sep = "\t", header = F)

#Keep only the syn file part of the dataset
sync_clean <- sync[,c(8:100)]
ncol <- dim(sync_clean)[2]
nloci <- dim(sync_clean)[1]
start_col <- 4
npool <- ncol - start_col #88 pools

lookup_table <- c("A" = 1, "T" = 2, "C" = 3, "G" = 4, "N" = 5)

#prepare output file
AF_mat <- matrix(NA, nloci, npool) #keep track of allele frequencies for each pool

#Loop through loci
for (locus in 1:nloci) {
  #locus <- 1
  alt <- sync_clean[locus,4]
  #Loop through pools
  for (pool in 1:npool) {
    #pool <- 1
    cur_pool <- as.numeric(unlist(strsplit(sync_clean[locus, start_col + pool],":")))
    #calculate sum
    total_num_reads <- sum(cur_pool)
    #if total reads < 5 --> NA
    if (total_num_reads <= 5 | is.na(total_num_reads)) { AF_mat[locus, pool] <- NA}
    #calculate AF based on alternate allele
    else if (total_num_reads > 5) {AF_mat[locus, pool] <- cur_pool[lookup_table[alt]]/total_num_reads}
  }
}

#add chrom, pos, pos, ALT to the AF matrix before writing out
pool_name <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/3.OnField/3.Data/metadata/pool_order.txt", 
                        sep = "\t", header = F)
AF_mat_export <- cbind(sync_clean[,1:4], AF_mat)
colnames(AF_mat_export) <- pool_name
write.table(AF_mat_export, "drought_loci_AF_matrix_2024_89pools.csv", sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
