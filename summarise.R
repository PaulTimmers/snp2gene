#!/usr/bin/Rscrpit

#----
# Setup environment
#----

library(data.table)
args <- commandArgs(T)

genes <- fread(args[1])
snps <- fread(args[2], header=F, col.names="rsid", sep="\t")
a <- genes[snps, on="rsid"]

#----
# Summarise
#----

# SNPs with no matches at all
answ <- a[is.na(cyto),.(rsid,gene=NA)]
a <- a[!rsid %in% answ$rsid,]


# SNPs without any gene names
answ <- rbind(answ, a[is.na(dist1),.(rsid,gene=cyto)])
a <- a[!rsid %in% answ$rsid,]


# SNPs with only one gene name and SNPs inside genes
answ <- rbind(answ, a[is.na(gene2) | dist1==0 | sign(dist1) == sign(dist2),.(rsid,gene=gene1)])
a <- a[!rsid %in% answ$rsid,]


# SNPs with two gene names
answ <- rbind(answ, a[, .(rsid, gene=ifelse(dist1 <= dist2, paste(gene1,gene2,sep="/"), paste(gene2,gene1,sep="/")))])
a <- a[!rsid %in% answ$rsid,]


# Shorten
answ[,gene:=gsub("/LOC.*$|/MIR|/LINC.*$","",gene)]
answ[,gene:=gsub("^.*LOC.*/|^.*MIR.*/|^.*LINC.*/","",gene)]


genes$summary <- answ[match(genes$rsid,answ$rsid),gene]

# Export
write.table(genes, args[1], quote=FALSE, sep="\t", row.names=FALSE)