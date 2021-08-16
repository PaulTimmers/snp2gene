#!/usr/bin/env Rscript

#----
# Setup environment
#----

library(data.table)

# Functions

comb <- function(gene1, gene2) {

l1 <- length(gene1)
summary <- array(NA, l1)

for(i in 1:l1) {
	
	g1 <- gene1[i]
	g2 <- gene2[i]
	
	for(j in 1:nchar(g1)) {
		overlap <- grepl(paste0("^", substr(g1,0,j)), g2)
		if(j == nchar(g1) && overlap) {j <- j+2; break}
		if(!overlap) break
	}
	summary[i] <- paste(g1,substr(g2,j-1,nchar(g2)),sep="/")

	}
return(summary)
}



# Variable

args <- commandArgs(T)
genes <- fread(args[1])
snps <- fread(args[2], header=F, col.names="rsid", sep="\t")
a <- genes[snps, on="rsid"]


#----
# Adjust small elements
#----

# Remove small elements from list, as long as not all genes nearby are small elements
a[!is.na(gene2) & grepl("LOC|MIR|LINC",gene1) & !(grepl("LOC|MIR|LINC",gene2) & grepl("LOC|MIR|LINC",gene3)) ,c("gene1","dist1","gene2","dist2","gene3","dist3"):=.(gene2,dist2,gene3,dist3,NA,NA)]

# Repeat in case second element is also small
a[!is.na(gene2) & grepl("LOC|MIR|LINC",gene1) & !(grepl("LOC|MIR|LINC",gene2) & grepl("LOC|MIR|LINC",gene3)) ,c("gene1","dist1","gene2","dist2","gene3","dist3"):=.(gene2,dist2,gene3,dist3,NA,NA)]

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
answ <- rbind(answ, a[is.na(gene2) | (dist1==0 & dist2 !=0) | sign(dist1) == sign(dist2) & sign(dist2) != 0 ,.(rsid,gene=gene1)])
a <- a[!rsid %in% answ$rsid,]


# SNPs with two gene names
answ <- rbind(answ, a[, .(rsid, gene=ifelse(dist1 <= dist2, comb(gene1,gene2), comb(gene2,gene1)))])
a <- a[!rsid %in% answ$rsid,]


# Shorten
answ[,gene:=gsub("/LOC.*$|/MIR.*$|/LINC.*$","",gene)]
answ[,gene:=gsub("^.*LOC.*/|^.*MIR.*/|^.*LINC.*/","",gene)]


genes$summary <- answ[match(genes$rsid,answ$rsid),gene]

# Export
write.table(genes, args[1], quote=FALSE, sep="\t", row.names=FALSE)
