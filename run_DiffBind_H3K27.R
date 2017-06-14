# nathan dot lazar at gmail dot com

# Runs both an occupancy analysis and differential binding analysis
# on samples listed in the input file.
# A number of plots are generated which are stored in the supplied 
# output directory

# Usage: run_DiffBind.R <sample file> <pooled sample file>
# <pseudo sample file> <pseudo pooled sample files> 
# <fragment length> <output directory>
# The sample sheet should be a comma separated file with a row for 
# each sample. Fields include Sample ID, Tissue, Factor,
# Condition, Treatment, Replicate, Read .bam file ,ControlID,
# Control .bam file, Peak .bed.gz file
# Note: RefSeq genes are read in using library GenomicFeatures

# Example: run_DiffBind.R H3K27me3_samples.csv H3K27me3_pool_samples.csv
#                         H3K27me3_pseudo_samples.csv 
#                         H3K27me3_ps_pool_samples.csv
#                         150  H3K27_analysis/


.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.2")
# library(ChIPQC)
library(dplyr)
library(DiffBind)
library(GenomicFeatures)
library(AnnotationDbi)
# library(org.Mm.eg.db)

args <- commandArgs(TRUE)
all.samples.csv <- args[1]
pool.samples.csv <- args[2]
pseudo.samples.csv <- args[3]
ps.pool.samples.csv <- args[4]
frag.size <- args[5]
out.dir <- args[6]
gene_file <- 'mm10_refseq.bed'

# all.samples.csv <- 'H3K27me3_samples.csv'
# pool.samples.csv <- 'H3K27me3_pool_samples.csv'
# pseudo.samples.csv <- 'H3K27me3_pseudo_samples.csv'
# ps.pool.samples.csv <- 'H3K27me3_ps_pool_samples.csv'
# frag.size <- 150
# out.dir <- 'H3K27_analysis/'

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

# Read in sample sheets (which were written by hand)
all.samples <- read.csv(all.samples.csv)
pool.samples <- read.csv(pool.samples.csv)
pseudo.samples <- read.csv(pseudo.samples.csv)
ps.pool.samples <- read.csv(ps.pool.samples.csv)

# Merge peak sets to get a master list of peaks and get a 
# binding matrix with a value for each peak and sample
all.samples.dba <- dba(sampleSheet=all.samples.csv, peakCaller='raw')
pool.samples.dba <- dba(sampleSheet=pool.samples.csv, peakCaller='raw')
pseudo.samples.dba <- dba(sampleSheet=pseudo.samples.csv, peakCaller='raw')
ps.pool.samples.dba <- dba(sampleSheet=ps.pool.samples.csv, peakCaller='raw')

# Close automatic plotting
dev.off()

print(all.samples.dba)
print(pool.samples.dba)
print(pseudo.samples.dba)
print(ps.pool.samples.dba)

# Function to plot the distribution of log(peak sizes)
plot_peak_hist <- function(peaks, sample, file) {
  png(file=paste0(out.dir, file))
  hist(log(peaks$V3 - peaks$V2), 
       main=paste(sample, 'peak size distribution'),
       xlab='log peak size', ylab='density', breaks=200)
  dev.off()
}

# Plot distributions of log peak sizes for each sample
for(i in 1:8)
  plot_peak_hist(all.samples.dba$peaks[[i]], all.samples.dba$samples$SampleID[[i]], 
                 paste0(all.samples.dba$samples$SampleID[[i]], '_peak_hist.png'))

# Plot the distribution of log peak sizes for pooled samples
plot_peak_hist(pool.samples.dba$peaks[[1]], pool.samples.dba$samples$SampleID[[1]], 
               paste0(pool.samples.dba$samples$SampleID[[1]], '_peak_hist.png'))
plot_peak_hist(pool.samples.dba$peaks[[2]], pool.samples.dba$samples$SampleID[[2]], 
               paste0(pool.samples.dba$samples$SampleID[[2]], '_peak_hist.png'))

# Plot distributions of log peak sizes for pseudo samples and pooled pseudo samples
for(i in 1:8)
  plot_peak_hist(pseudo.samples.dba$peaks[[i]], pseudo.samples.dba$samples$SampleID[[i]],
                 paste0(pseudo.samples.dba$samples$SampleID[[i]], '_peak_hist.png'))
plot_peak_hist(ps.pool.samples.dba$peaks[[1]], 'Alc_pseudo_pool',
               'Alc_pseudo_pool_peak_hist.png')
plot_peak_hist(ps.pool.samples.dba$peaks[[2]], 'Con_pseudo_pool',
               'Con_pseudo_pool_peak_hist.png')

# Plot chromosomal distribution of peaks for all samples
plot_chrom_dist <- function(peaks, sample, file) {
  df <- peaks %>% tbl_df %>% group_by(V1) %>% 
    summarise(length(V2)) %>% data.frame
  vec <- df[,2]
  names(vec) <- df[,1]
  png(file=paste0(out.dir, file))
  barplot(vec, main=paste(sample, 'Peak distribution'),
    ylab='Peaks', las=2)
  dev.off()
}

for(i in 1:8)
  plot_chrom_dist(all.samples.dba$peaks[[i]], all.samples.dba$samples$SampleID[[i]],
                  paste0(all.samples.dba$samples$SampleID[[i]], '_chrom_dist.png'))
plot_chrom_dist(pool.samples.dba$peaks[[1]], pool.samples.dba$samples$SampleID[[1]],
                paste0(pool.samples.dba$samples$SampleID[[1]], '_chrom_dist.png'))
plot_chrom_dist(pool.samples.dba$peaks[[2]], pool.samples.dba$samples$SampleID[[2]],
                paste0(pool.samples.dba$samples$SampleID[[2]], '_chrom_dist.png'))

############################################################
# 1) Occupancy analysis: Look at differences in called peaks
#    between samples and conditions

all.olap.rate <- dba.overlap(all.samples.dba, mode=DBA_OLAP_RATE)
pool.olap.rate <- dba.overlap(pool.samples.dba, mode=DBA_OLAP_RATE)
pseudo.olap.rate <- dba.overlap(pseudo.samples.dba, mode=DBA_OLAP_RATE)
ps.pool.olap.rate <- dba.overlap(ps.pool.samples.dba, mode=DBA_OLAP_RATE)

print("All overlap rates:"); print(all.olap.rate)
print("Pooled samples overlap rates:"); print(pool.olap.rate)
print("Pseudo replicate overlap rates:"); print(pseudo.olap.rate)
print("Pseudo replicate pooled overlap rates:"); print(ps.pool.olap.rate)

# Plot the overlap rates for different numbers of samples
png(paste0(out.dir, 'olap_rate.png'))
par(mfrow=c(2,2))
plot(all.olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='All samples overlap rate')
plot(pool.olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Pooled samples overlap rate')
plot(pseudo.olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Pseudo replicates overlap rate')
plot(ps.pool.olap.rate, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Pseudo pooled overlap rate')
dev.off()

# Plot the overlap rates for just alcohol samples
olap.alc <- dba.overlap(all.samples.dba, all.samples.dba$masks$Alcohol, mode=DBA_OLAP_RATE)
olap.ps.alc <- dba.overlap(pseudo.samples.dba, pseudo.samples.dba$masks$Alcohol, mode=DBA_OLAP_RATE)
png(paste0(out.dir, 'olap_rate_alc.png'), width=480*2)
par(mfrow=c(1,2))
plot(olap.alc, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Alcohol samples overlap rates')
plot(olap.ps.alc, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Pseudo rep. alcohol overlap rates')
dev.off()

# Plot the overlap rates for just control samples
olap.con <- dba.overlap(all.samples.dba, all.samples.dba$masks$Control, mode=DBA_OLAP_RATE)
olap.ps.con <- dba.overlap(pseudo.samples.dba, pseudo.samples.dba$masks$Control, mode=DBA_OLAP_RATE)
png(paste0(out.dir, 'olap_rate_con.png'), width=480*2)
par(mfrow=c(1,2))
plot(olap.con, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Control samples overlap rates')
plot(olap.ps.con, type='b', ylab='# peaks', xlab='Overlap at least this many peaksets',
     main='Pseudo rep. control overlap rates')
dev.off()

# Venn diagrams for the two conditions separately
png(paste0(out.dir, 'Alc_venn.png'))
dba.plotVenn(all.samples.dba, all.samples.dba$masks$Alcohol, main='Alcohol')
dev.off()
png(paste0(out.dir, 'Pseudo_alc_venn.png'))
dba.plotVenn(pseudo.samples.dba, pseudo.samples.dba$masks$Alcohol, main='Pseudo alcohol')
dev.off()
png(paste0(out.dir, 'Con_venn.png'))
dba.plotVenn(all.samples.dba, all.samples.dba$masks$Control, main='Control')
dev.off()
png(paste0(out.dir, 'Pseudo_con_venn.png'))
dba.plotVenn(pseudo.samples.dba, pseudo.samples.dba$masks$Control, main='Control')
dev.off()

# Generate a peakset for each condition where peaks are in all 4 samples
all.samples.dba <- dba.peakset(all.samples.dba, consensus=-DBA_REPLICATE, minOverlap=4)
pseudo.samples.dba <- dba.peakset(pseudo.samples.dba, consensus=-DBA_REPLICATE, minOverlap=4)

# Venn diagram between conditions
png(paste0(out.dir, 'Alc_vs_Con_Venns.png'), height=480*1.5, width=480*1.5)
par(mfrow=c(2,2))
dba.plotVenn(all.samples.dba, mask=all.samples.dba$masks$Consensus,
             main='All samples')
dba.plotVenn(pool.samples.dba, mask=c(1,2), main='Pooled samples')
dba.plotVenn(pseudo.samples.dba, mask=pseudo.samples.dba$masks$Consensus,
             main='Pseudo samples')
dba.plotVenn(ps.pool.samples.dba, mask=c(1,2),
             main='Pooled pseudo samples')
dev.off()

# Using only peak caller data make a correlation heatmap.
# This gives a clustering of the samples based on 
# "occupancy scores" using the cross-correlations of each 
# row of the binding matrix:
pdf(paste0(out.dir, 'peak_cc_heatmap.pdf'))
plot(all.samples.dba)
dev.off()

# Count reads in each peak for each sample and get FRiP
all.samples.dba <- dba.count(all.samples.dba, minOverlap=1,
                             fragmentSize=frag.size)
pool.samples.dba <- dba.count(pool.samples.dba, minOverlap=1,
                              fragmentSize=frag.size)

print(all.samples.dba)
print(pool.samples.dba)
print(pseudo.samples.dba)
print(ps.pool.samples.dba)

# There's a problem with chromsome names below. 
# Reload pool.samples.dba
pool.samples.dba <- dba(sampleSheet=pool.samples.csv, peakCaller='raw')

# Retrieve the peaks only found in alcohol or control GRanges objects
diff.gr <- dba.overlap(pool.samples.dba, pool.samples.dba$masks$All)
alc.gr <- diff.gr$onlyA
con.gr <- diff.gr$onlyB
both.gr <- diff.gr$inAll

# Save the granges objects
save(alc.gr, file=paste0(out.dir, 'alc_gr.Rdata'))
save(con.gr, file=paste0(out.dir, 'con_gr.Rdata'))
save(both.gr, file=paste0(out.dir, 'both_gr.Rdata'))

# There's a problem with chromsome names below.
# Reload all.samples.dba
re.all.samples.dba <- dba(sampleSheet=pool.samples.csv, peakCaller='raw')

all.alc <- dba.overlap(re.all.samples.dba, re.all.samples.dba$masks$Alcohol)
all.con <- dba.overlap(re.all.samples.dba, re.all.samples.dba$masks$Control)

# Get a list of overlaps of peaks in all alcohol samples with each object 
# in the list of control peak overlaps and vice vers
alc.laps <- lapply(all.con, function(x) findOverlaps(all.alc$inAll, x))
con.laps <- lapply(all.alc, function(x) findOverlaps(all.con$inAll, x))

# Make a GRanges object of peaks that are in all the alcohol samples
# and 0,1,2,3 or 4 of the control samples
all.alc.list <- list()
n <- length(all.alc$inAll)
which.lap.0 <- (1:n)[-unique(unlist(sapply(alc.laps, function(x) x@queryHits)))]
which.lap.1 <- (1:n)[-unique(c(which.lap.0,
  alc.laps[['AandB']]@queryHits, alc.laps[['AandC']]@queryHits, alc.laps[['AandD']]@queryHits,
  alc.laps[['BandC']]@queryHits, alc.laps[['BandD']]@queryHits, alc.laps[['CandD']]@queryHits,
  alc.laps[['notA']]@queryHits, alc.laps[['notB']]@queryHits, alc.laps[['notC']]@queryHits,
  alc.laps[['inAll']]@queryHits))]
which.lap.2 <- (1:n)[-unique(c(which.lap.0, which.lap.1,
  alc.laps[['notA']]@queryHits, alc.laps[['notB']]@queryHits, alc.laps[['notC']]@queryHits,
  alc.laps[['inAll']]@queryHits))]
which.lap.3 <- (1:n)[-unique(c(which.lap.0, which.lap.1, which.lap.2,
  alc.laps[['inAll']]@queryHits))]

all.alc.list[['con.0']] <- all.alc$inAll[which.lap.0]
all.alc.list[['con.1']] <- all.alc$inAll[which.lap.1]
all.alc.list[['con.2']] <- all.alc$inAll[which.lap.2]
all.alc.list[['con.3']] <- all.alc$inAll[which.lap.3]
all.alc.list[['con.4']] <- subsetByOverlaps(all.alc$inAll, all.con$inAll)

# Add the full set of peaks in alcohol to the list
all.alc.list[['all.peaks']] <- all.alc$inAll

# Make Granges object of peaks in all control and 0,1,2,3,4 alcohol samples
all.con.list <- list()
n <- length(all.con$inAll)
which.lap.0 <- (1:n)[-unique(unlist(sapply(con.laps, function(x) x@queryHits)))]
which.lap.1 <- (1:n)[-unique(c(which.lap.0,
  con.laps[['AandB']]@queryHits, con.laps[['AandC']]@queryHits, con.laps[['AandD']]@queryHits,
  con.laps[['BandC']]@queryHits, con.laps[['BandD']]@queryHits, con.laps[['CandD']]@queryHits,
  con.laps[['notA']]@queryHits, con.laps[['notB']]@queryHits, con.laps[['notC']]@queryHits,
  con.laps[['inAll']]@queryHits))]
which.lap.2 <- (1:n)[-unique(c(which.lap.0, which.lap.1,
  con.laps[['notA']]@queryHits, con.laps[['notB']]@queryHits, con.laps[['notC']]@queryHits,
  con.laps[['inAll']]@queryHits))]
which.lap.3 <- (1:n)[-unique(c(which.lap.0, which.lap.1, which.lap.2,
  con.laps[['inAll']]@queryHits))]

all.con.list[['alc.0']] <- all.con$inAll[which.lap.0]
all.con.list[['alc.1']] <- all.con$inAll[which.lap.1]
all.con.list[['alc.2']] <- all.con$inAll[which.lap.2]
all.con.list[['alc.3']] <- all.con$inAll[which.lap.3]
all.con.list[['alc.4']] <- subsetByOverlaps(all.con$inAll, all.alc$inAll)

# Add the full set of peaks in alcohol to the list
all.con.list[['all.peaks']] <- all.con$inAll

# Find genes overlapping these peak sets.
#########################################
# Download genes and make GRanges object
# genes.txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")
# genes.gr <- genes(genes.txdb)
# genes.gr$symbol <- mapIds(org.Mm.eg.db, keys=genes.gr$gene_id,
#                          column='SYMBOL', keytype='ENTREZID', multiVals='first')

genes <- read.table(gene_file)
genes.gr <- makeGRangesFromDataFrame(genes, keep.extra.columns=T,
  seqnames.field='V1', start.field='V2', end.field='V3')
mcols(genes.gr)[[1]] <- as.vector(mcols(genes.gr)[[1]])
names(mcols(genes.gr))[1] <- 'gene_id'

alc.peaks.lap.genes.list <- lapply(all.alc.list, subsetByOverlaps, genes.gr)
con.peaks.lap.genes.list <- lapply(all.con.list, subsetByOverlaps, genes.gr)

# Function to get gene ids for overlaps
get_lap_ids <- function(query, subject) {
  laps <- findOverlaps(query, subject)
  lap.df <- tbl_df(data.frame(q_idx=laps@queryHits, 
                              s_idx=laps@subjectHits,
                              s_id=subject$gene_id[laps@subjectHits]))
  id.df <- group_by(lap.df, q_idx) %>%
    summarise(ids=paste(s_id, sep=',', collapse=','))
  return(id.df$ids)
}

# # Function to get gene symbols for overlaps
# get_symbols <- function(query, subject) {
#   laps <- findOverlaps(query, subject)
#   lap.df <- tbl_df(data.frame(q_idx=laps@queryHits, 
#                               s_idx=laps@subjectHits,
#                               sym=subject$symbol[laps@subjectHits]))
#   symbol.df <- group_by(lap.df, q_idx) %>%
#     summarise(symbols=paste(sym, sep=',', collapse=','))
#   return(symbol.df$symbols)
# }

# Add gene ids & symbols to these GRanges as metadata

for(i in 1:length(alc.peaks.lap.genes.list)) {
  alc.peaks.lap.genes.list[[i]]$gene_ids <- 
    get_lap_ids(all.alc.list[[i]], genes.gr)
#  alc.peaks.lap.genes.list[[i]]$gene.symbol <- 
#    get_symbols(all.alc.list[[i]], genes.gr)
}

for(i in 1:length(con.peaks.lap.genes.list)) {
  con.peaks.lap.genes.list[[i]]$gene_ids <- 
    get_lap_ids(all.con.list[[i]], genes.gr)
#  con.peaks.lap.genes.list[[i]]$gene.symbol <- 
#    get_symbols(all.con.list[[i]], genes.gr)
}

alc.genes.lap.peaks.list <- 
  lapply(all.alc.list, function(x) subsetByOverlaps(genes.gr, x))

con.genes.lap.peaks.list <- 
  lapply(all.con.list, function(x) subsetByOverlaps(genes.gr, x))

# Sort the above GRanges objects
tmp <- sapply(alc.genes.lap.peaks.list, sort)
alc.genes.lap.peaks.list <- tmp
tmp <- sapply(con.genes.lap.peaks.list, sort)
con.genes.lap.peaks.list <- tmp

# Find genes that have peaks unique to alcohol and peaks unique to contols
alc.and.con.list <- intersect(c(alc.genes.lap.peaks.list[['con.0']]$gene_id,
                                alc.genes.lap.peaks.list[['con.1']]$gene_id,
                                alc.genes.lap.peaks.list[['con.2']]$gene_id,
                                alc.genes.lap.peaks.list[['con.3']]$gene_id),
                              con.genes.lap.peaks.list[['alc.0']]$gene_id)
con.and.alc.list <- intersect(c(con.genes.lap.peaks.list[['alc.0']]$gene_id,
                                con.genes.lap.peaks.list[['alc.1']]$gene_id,
                                con.genes.lap.peaks.list[['alc.2']]$gene_id,
                                con.genes.lap.peaks.list[['alc.3']]$gene_id),
                              alc.genes.lap.peaks.list[['con.0']]$gene_id)

# Make a data frame of peaks found in all 4 alcohol samples
# sorted by the number of control samples that also have these peaks 
alc.peaks.df <- data.frame(chr=seqnames(all.alc.list[['con.0']]),
                           start=start(all.alc.list[['con.0']])-1,
                           end=end(all.alc.list[['con.0']]),
                           lap.con=0)
alc.peaks.df <- rbind(alc.peaks.df, 
  data.frame(chr=seqnames(all.alc.list[['con.1']]),
             start=start(all.alc.list[['con.1']])-1,
             end=end(all.alc.list[['con.1']]),
             lap.con=1))
alc.peaks.df <- rbind(alc.peaks.df, 
  data.frame(chr=seqnames(all.alc.list[['con.2']]),
             start=start(all.alc.list[['con.2']])-1,
             end=end(all.alc.list[['con.2']]),
             lap.con=2))
alc.peaks.df <- rbind(alc.peaks.df, 
  data.frame(chr=seqnames(all.alc.list[['con.3']]),
             start=start(all.alc.list[['con.3']])-1,
             end=end(all.alc.list[['con.3']]),
             lap.con=3))

# Make a data frame of peaks found in all 4 control samples
# sorted by the number of alcohol samples that also have these peaks
con.peaks.df <- data.frame(chr=seqnames(all.con.list[['alc.0']]),
                           start=start(all.con.list[['alc.0']])-1,
                           end=end(all.con.list[['alc.0']]),
                           lap.alc=0)
con.peaks.df <- rbind(con.peaks.df, 
  data.frame(chr=seqnames(all.con.list[['alc.1']]),
             start=start(all.con.list[['alc.1']])-1,
             end=end(all.con.list[['alc.1']]),
             lap.alc=1))
con.peaks.df <- rbind(con.peaks.df, 
  data.frame(chr=seqnames(all.con.list[['alc.2']]),
             start=start(all.con.list[['alc.2']])-1,
             end=end(all.con.list[['alc.2']]),
             lap.alc=2))
con.peaks.df <- rbind(con.peaks.df, 
  data.frame(chr=seqnames(all.con.list[['alc.3']]),
             start=start(all.con.list[['alc.3']])-1,
             end=end(all.con.list[['alc.3']]),
             lap.alc=3))

# Make a data frame of genes with peaks in all 4 alcohol samples
# sorted by the number of control samples that also have these peaks
# and with a column to indicate whether the gene also has a peak that
# is in all 4 control samples and no alcohol samples

alc.genes.df <- data.frame(chr=seqnames(alc.genes.lap.peaks.list[['con.0']]),
                           start=start(alc.genes.lap.peaks.list[['con.0']]),
                           end=end(alc.genes.lap.peaks.list[['con.0']]),
                           gene_id=alc.genes.lap.peaks.list[['con.0']]$gene_id,
#                           symbol=alc.genes.lap.peaks.list[['con.0']]$symbol,
                           lap.con=0)
alc.genes.df <- rbind(alc.genes.df,
  data.frame(chr=seqnames(alc.genes.lap.peaks.list[['con.1']]),
             start=start(alc.genes.lap.peaks.list[['con.1']]),
             end=end(alc.genes.lap.peaks.list[['con.1']]),
             gene_id=alc.genes.lap.peaks.list[['con.1']]$gene_id,
#             symbol=alc.genes.lap.peaks.list[['con.1']]$symbol,
             lap.con=1))
alc.genes.df <- rbind(alc.genes.df,
  data.frame(chr=seqnames(alc.genes.lap.peaks.list[['con.2']]),
             start=start(alc.genes.lap.peaks.list[['con.2']]),
             end=end(alc.genes.lap.peaks.list[['con.2']]),
             gene_id=alc.genes.lap.peaks.list[['con.2']]$gene_id,
#             symbol=alc.genes.lap.peaks.list[['con.2']]$symbol,
             lap.con=2))
alc.genes.df <- rbind(alc.genes.df,
  data.frame(chr=seqnames(alc.genes.lap.peaks.list[['con.3']]),
             start=start(alc.genes.lap.peaks.list[['con.3']]),
             end=end(alc.genes.lap.peaks.list[['con.3']]),
             gene_id=alc.genes.lap.peaks.list[['con.3']]$gene_id,
#             symbol=alc.genes.lap.peaks.list[['con.3']]$symbol,
             lap.con=3))

# Add a column indicating whether the gene also has peaks that are unique 
# to control samples 
alc.genes.df$con.peaks <- ' '
alc.genes.df$con.peaks[alc.genes.df$gene_id %in% alc.and.con.list] <- 'Y'

# Make a data frame of genes with peaks in all 4 control samples
# sorted by the number of alcohol samples that also have these peaks
# and with a column to indicate whether the gene also has a peak that
# is in all 4 alcohol samples and no control samples

con.genes.df <- data.frame(chr=seqnames(con.genes.lap.peaks.list[['alc.0']]),
                           start=start(con.genes.lap.peaks.list[['alc.0']]),
                           end=end(con.genes.lap.peaks.list[['alc.0']]),
                           gene_id=con.genes.lap.peaks.list[['alc.0']]$gene_id,
#                           symbol=con.genes.lap.peaks.list[['alc.0']]$symbol,
                           lap.alc=0)
con.genes.df <- rbind(con.genes.df,
  data.frame(chr=seqnames(con.genes.lap.peaks.list[['alc.1']]),
             start=start(con.genes.lap.peaks.list[['alc.1']]),
             end=end(con.genes.lap.peaks.list[['alc.1']]),
             gene_id=con.genes.lap.peaks.list[['alc.1']]$gene_id,
#             symbol=con.genes.lap.peaks.list[['alc.1']]$symbol,
             lap.alc=1))
con.genes.df <- rbind(con.genes.df,
  data.frame(chr=seqnames(con.genes.lap.peaks.list[['alc.2']]),
             start=start(con.genes.lap.peaks.list[['alc.2']]),
             end=end(con.genes.lap.peaks.list[['alc.2']]),
             gene_id=con.genes.lap.peaks.list[['alc.2']]$gene_id,
#             symbol=con.genes.lap.peaks.list[['alc.2']]$symbol,
             lap.alc=2))
con.genes.df <- rbind(con.genes.df,
  data.frame(chr=seqnames(con.genes.lap.peaks.list[['alc.3']]),
             start=start(con.genes.lap.peaks.list[['alc.3']]),
             end=end(con.genes.lap.peaks.list[['alc.3']]),
             gene_id=con.genes.lap.peaks.list[['alc.3']]$gene_id,
#             symbol=con.genes.lap.peaks.list[['alc.3']]$symbol,
             lap.alc=3))

# Add a column indicating whether the gene also has peaks that are unique 
# to control samples 
con.genes.df$alc.peaks <- ' '
con.genes.df$alc.peaks[con.genes.df$gene_id %in% con.and.alc.list] <- 'Y'

# Function to write GRanges object to a tab delimited file
write_gr <- function(gr, file) {
  df <- data.frame(chr=seqnames(gr),
                   start=start(gr)-1,
                   end=end(gr))
  df <- cbind(df, mcols(gr))
  write.table(df, file=file, quote=F, sep='\t', row.names=F, col.names=T)
}

# Write differentially bound peak data frames to tab delmited files
write.table(alc.peaks.df, file=paste0(out.dir, "alc_specific_peaks.txt"),
            quote=F, row.names=F, col.names=T, sep='\t')
write.table(con.peaks.df, file=paste0(out.dir, "con_specific_peaks.txt"),
            quote=F, row.names=F, col.names=T, sep='\t')

# Write genes w/ differentially bound peaks to files
write.table(alc.genes.df, file=paste0(out.dir, "alc_genes_lap_peaks.txt"),
            quote=F, row.names=F, col.names=T, sep='\t')
write.table(con.genes.df, file=paste0(out.dir, "con_genes_lap_peaks.txt"),
            quote=F, row.names=F, col.names=T, sep='\t')

save.image(paste0(out.dir, 'image.Rdata'))

# Differential binding analysis:
# This isn't recommended for broad peaks and the regions found did not 
# seem to be truely differentially bound when looking on IGV
#################################################################################

# Depreciated
#############
# Get report for all ChIP samples using ChIPQC
# all.qc <- ChIPQC(all.samples, 'mm10', fragmentLength=frag.size)
# ChIPQCreport(all.qc, reportName='H3K17', facetBy='Factor',
#             colourBy='Replicate',
#              reportFolder=paste0(out.dir, 'all_samples'))
# 
# Get a DBA object of just the consensus peaks 
# in all samples for each group: Alc vs. Con
# consensus.dba <- dba(all.samples.dba, mask=all.samples.dba$masks$Consensus,
#                      minOverlap=4)
# con2.dba <- dba(all.samples.dba, minOverlap=4,
#                 mask=!all.samples.dba$masks$Consensus)
#
# Get a list of GRanges objects of peaks unique to each group
# all.samples.GRlist = dba.overlap(all.samples.dba, all.samples.dba$masks$Consensus)
# pool.samples.GRlist = dba.overlap(pool.samples.dba, mask=c(1,2))
# all.samples.GRlist$onlyA
# all.samples.GRlist$onlyB
# pool.samples.GRlist$onlyA
# pool.samples.GRlist$onlyB

# Subset the alcohol and control specific peaks that overlap genes
# alc.peaks.lap.genes.gr <- subsetByOverlaps(alc.gr, genes.gr)
# con.peaks.lap.genes.gr <- subsetByOverlaps(con.gr, genes.gr)
# both.peaks.lap.genes.gr <- subsetByOverlaps(both.gr, genes.gr)

# Add gene ids & symbols to these GRanges as metadata

# alc.peaks.lap.genes.gr$gene_ids <- get_lap_ids(alc.gr, genes.gr)
# con.peaks.lap.genes.gr$gene_ids <- get_lap_ids(con.gr, genes.gr)
# both.peaks.lap.genes.gr$gene_ids <- get_lap_ids(both.gr, genes.gr)

# alc.peaks.lap.genes.gr$symbol <- get_symbols(alc.gr, genes.gr)
# con.peaks.lap.genes.gr$symbol <- get_symbols(con.gr, genes.gr)
# both.peaks.lap.genes.gr$symbol <- get_symbols(both.gr, genes.gr)

# Subset the alcohol and control specific genes that overlap peaks
# alc.genes.lap.peaks.gr <- subsetByOverlaps(genes.gr, alc.gr)
# con.genes.lap.peaks.gr <- subsetByOverlaps(genes.gr, con.gr)
# both.genes.lap.peaks.gr <- subsetByOverlaps(genes.gr, both.gr)

# Subset the genes uniquely overlapping alcohol or control peaks
# alc.only.genes.lap.peaks.gr <- alc.genes.lap.peaks.gr[!(alc.genes.lap.peaks.gr$gene_id %in%
#   con.genes.lap.peaks.gr$gene_id)]
# con.only.genes.lap.peaks.gr <- con.genes.lap.peaks.gr[!(con.genes.lap.peaks.gr$gene_id %in%
#   alc.genes.lap.peaks.gr$gene_id)]

# Sort the above GRanges objects
# alc.genes.lap.peaks.gr <- sort(alc.genes.lap.peaks.gr)
# con.genes.lap.peaks.gr <- sort(con.genes.lap.peaks.gr)
# alc.only.genes.lap.peaks.gr <- sort(alc.only.genes.lap.peaks.gr)
# con.only.genes.lap.peaks.gr <- sort(con.only.genes.lap.peaks.gr)

# Function to write GRanges object to a tab delimited file
# write_gr <- function(gr, file) {
#   df <- data.frame(chr=seqnames(gr),
#                    start=start(gr)-1,
#                    end=end(gr))
#   df <- cbind(df, mcols(gr))
#   write.table(df, file=file, quote=F, sep='\t', row.names=F, col.names=T)
# }
# 
# # Write differentially bound GRanges objects to a tab delmited file
# write_gr(alc.gr, paste0(out.dir, "alc_specific_peaks.txt"))
# write_gr(con.gr, paste0(out.dir, "con_specific_peaks.txt"))
# 
# # Write differentially bound peaks overlapping genes
# write_gr(alc.peaks.lap.genes.gr, paste0(out.dir, "alc_peaks_lap_genes.txt"))
# write_gr(con.peaks.lap.genes.gr, paste0(out.dir, "con_peaks_lap_genes.txt"))
# 
# # Write genes w/ differentially bound peaks to files
# write_gr(alc.genes.lap.peaks.gr, paste0(out.dir, "alc_genes_lap_peaks.txt"))
# write_gr(con.genes.lap.peaks.gr, paste0(out.dir, "con_genes_lap_peaks.txt"))
# 
# # Write genes w/ differentially bound peaks unique to one condition to files
# write_gr(alc.only.genes.lap.peaks.gr, paste0(out.dir, "alc_only_genes_lap_peaks.txt"))
# write_gr(con.only.genes.lap.peaks.gr, paste0(out.dir, "con_only_genes_lap_peaks.txt"))

