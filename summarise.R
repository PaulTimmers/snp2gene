#!/usr/bin/Rscrpit

#----
# Setup environment
#----

library(data.table)
args <- commandArgs(T)

genes <- fread(args[1])
snps <- fread(args[2], header=F, col.names="rsid")
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
answ <- rbind(answ, a[is.na(gene2) | dist1==0 | dist1 == dist2 ,.(rsid,gene=gene1)])
a <- a[!rsid %in% answ$rsid,]


# SNPs with two gene names
answ <- rbind(answ, a[is.na(gene3), .(rsid, gene=ifelse(dist1 <= dist2, paste(gene1,gene2,sep="/"), paste(gene2,gene1,sep="/")))])
a <- a[!rsid %in% answ$rsid,]


# SNPs with three gene names
for (i in 1:nrow(a)) {
	b <- a[i,.(gene=c(gene1,gene2,gene3),dist=c(dist1,dist2,dist3))]
	b <- paste(b[order(dist),gene][c(1,3)],collapse="/")
	answ <- rbind(answ, data.table(rsid=a[i,rsid],gene=b))
}
a <- a[!rsid %in% answ$rsid,]


# Shorten
answ[,gene:=gsub("/LOC.*$|/MIR","",gene)]
answ[,gene:=gsub("^.*LOC.*/|^.*MIR.*/","",gene)]


# Export
write.table(answ, args[1], quote=FALSE, sep="\t", row.names=FALSE)