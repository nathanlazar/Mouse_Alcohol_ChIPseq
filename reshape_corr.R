# nathan dot lazar at gmail dot com

# Reads in correlations.txt with all pairwise correlations between 
# coverage.gz files. ENCODE "blacklisted" reagions have been removed
# before calculating corrlations

library(reshape2)
args <- commandArgs(TRUE)

corrs <- read.table(args[[1]], header=F,
                    stringsAsFactors=F, sep=" ",
                    col.names=c('samp1', 'vs', 'samp2', 'cor'))

corrs <- corrs[,-2]
corrs$samp2 <- sub(':', '', corrs$samp2)

cor.mat <- acast(corrs, samp1~samp2, value.var='cor')

write.table(cor.mat, file=args[[2]], quote=F,
            sep='\t')
