library(data.table)
library(dplyr)

dt <- unlist(fread("md5sum_cloud.tsv", sep=NULL, header=FALSE))
dt <- data.table(file_cloud=grep("Hashes", dt, value=TRUE), md5_cloud=grep("md5", dt, value=TRUE))

dt_cloud <- dt %>% mutate(
file_cloud=gsub("Hashes.*additive-tsvs/(.*):", "\\1", file_cloud),
md5_cloud=gsub(".*:\t\t", "", md5_cloud))

# Created the md5 sums on the cluster too (see additive_md5sum_cluster.sh in this folder for details).
# Copy them down
# system("scp dpalmer@login02:*md5sum_cluster.tsv .")

# Read in and munge
dt_md5_list <- list()
i <- 1
for(file in c("both_sexes_md5sum_cluster.tsv", "male_md5sum_cluster.tsv", "female_md5sum_cluster.tsv")) {
	dt_md5_list[[i]] <- fread(file, header=FALSE)
	i <- i+1
}
dt_cluster <- rbindlist(dt_md5_list) %>% rename(md5_cluster=V1, file_cluster=V2) %>% mutate(file_cluster = gsub(".*/", "", file_cluster))

setdiff(dt_cluster$file_cluster, dt_cloud$file_cloud)
setdiff(dt_cloud$file_cloud, dt_cluster$file_cluster)

dt <- merge(dt_cluster %>% rename(file=file_cluster), dt_cloud %>% rename(file=file_cloud), by="file")

# Check that they're all the same.
all(dt$md5_cloud == dt$md5_cluster)

# They are, great!
