
setwd("C:/Users/jamsh/OneDrive/Desktop/Data_Analysis/Project_1/Melica/Siri_data/")

library(stringr)
library(DESeq2)
library(dplyr)
library(tidyr)
library(edgeR)
library(ggplot2)
library(DEFormats)

######
# Computing values needed in time-profiles
######

countsTable <- read.table("countMatrix.rsem.txt", header= TRUE, row.names=1)
countsTableRound <- round(countsTable)
cnames <- names(countsTableRound)
for (i in 1:length(names(countsTableRound))){
  cnames[i] = str_replace(cnames[i], "expressions.", "")
  cnames[i] = str_replace(cnames[i], ".genes.results", "")
}
names(countsTableRound) <- cnames

filteredCountsTableRound <- countsTableRound[rowSums(countsTableRound)>10*ncol(countsTableRound),]


# ----------------------------------
## FR
# ----------------------------------
conds <- read.table("samples_FRaFRb_T0T3T4.txt", header=TRUE ,stringsAsFactors=TRUE, row.names=1)
t0t3t4 <- c("4.T0.FRa.4.T0.FRa",
            "5.T0.FRa.5.T0.FRa", 
            "6.T0.FRa.6.T0.FRa",
            "7.T0.FRa.7.T0.FRa",
            "8.T0.FRa.8.T0.FRa",
            "T3.FRa2.4.T3.FRa2.4",
            "T3.FRa5.2.T3.FRa5.2",        
            "T3.FRa7.3.T3.FRa7.3",
            "T3.FRa.T3.FRa",   
            "T4.FRa11.4.T4.FRa11.4",
            "T4.FRa14.4.T4.FRa14.4",
            "T4.FRa.T4.FRa",
            "T4.FRb11.3.T4.FRb11.3" 
)
sub_countsTableRound <- filteredCountsTableRound[, t0t3t4]


dds <- DESeqDataSetFromMatrix(countData = sub_countsTableRound, 
                              colData = conds, 
                              design = ~ Genotype)
normFactors <- calcNormFactors(sub_countsTableRound, method="TMM")
N <- colSums(sub_countsTableRound) #vector of library size
tmm.counts <- N*normFactors/exp(mean(log(N*normFactors)))
sizeFactors(dds) <- tmm.counts
dds <- DESeq(dds)
#
vsd <- vst(dds, blind=FALSE, nsub=100)
mat <- assay(vsd)
mat <- mat - rowMeans(mat)
genotypes <- as.data.frame(colData(dds))["Genotype"]
#
FR.profiles = data.frame(T0 = rowMeans(mat[,genotypes=="T0"]),
                         T3 = rowMeans(mat[,genotypes=="T3"]),
                         T4 = rowMeans(mat[,genotypes=="T4"]))
# normalizing the scale of time profiles
prof.max <- max(FR.profiles)
prof.min <- min(FR.profiles)
prof.rng = prof.max - prof.min
FR.profiles <- (FR.profiles - prof.min) / prof.rng



conds <- read.table("samples_BIdBIe_T0T3T4.txt", header=TRUE ,stringsAsFactors=TRUE, row.names=1)
t0t3t4 <- c("T0.10.BIdR.3.T0.10.BIdR.3",
            "T0.11.BIdR.4.T0.11.BIdR.4",
            "T0.12.BIdR.5.T0.12.BIdR.5",  
            "T0.9.BIdR.2.T0.9.BIdR.2",    
            "T3.45.BId2.T3.45.BId2",
            "T3.46.BId9.2.T3.46.BId9.2", 
            "T3.47.BId14.3.T3.47.BId14.3",
            "T4.48.BId9.3.T4.48.BId9.3",
            "T4.49.BIe15.3.T4.49.BIe15.3",
            "T4.50.BIe6.3.T4.50.BIe6.3")
sub_countsTableRound <- filteredCountsTableRound[, t0t3t4]

dds <- DESeqDataSetFromMatrix(countData = sub_countsTableRound, 
                              colData = conds, 
                              design = ~ Genotype)
normFactors <- calcNormFactors(sub_countsTableRound, method="TMM")
N <- colSums(sub_countsTableRound) #vector of library size
tmm.counts <- N*normFactors/exp(mean(log(N*normFactors)))
sizeFactors(dds) <- tmm.counts
dds <- DESeq(dds)
#
vsd <- vst(dds, blind=FALSE, nsub=100)
mat <- assay(vsd)
mat <- mat - rowMeans(mat)
genotypes <- as.data.frame(colData(dds))["Genotype"]
#
BI.profiles = data.frame(T0 = rowMeans(mat[,genotypes=="T0"]),
                         T3 = rowMeans(mat[,genotypes=="T3"]),
                         T4 = rowMeans(mat[,genotypes=="T4"]))
# normalizing the scale of time profiles
prof.max <- max(BI.profiles)
prof.min <- min(BI.profiles)
prof.rng = prof.max - prof.min
BI.profiles <- (BI.profiles - prof.min) / prof.rng


### concatenating the profiles
names(FR.profiles) <- c("FR.T0", "FR.T3", "FR.T4")
names(BI.profiles) <- c("BI.T0", "BI.T3", "BI.T4")
profiles <- cbind(FR.profiles, BI.profiles)

write.table(profiles, "KMeans/rev4_results/profiles.tsv", sep = '\t', quote = FALSE)



###############
# Clustering 
###############

gc()
km.res <- kmeans(profiles, 6, nstart = 10)

write(km.res$cluster, "KMeans/rev4_results/cluster_indices.txt", ncolumns=1)
km.res.inds <- read.table("KMeans/rev4_results/cluster_indices.txt", header=FALSE)


cluster.ind <- 2
cluster <- profiles[km.res$cluster==cluster.ind,]
cluster <- profiles[km.res.inds==cluster.ind,]
cluster.size <- dim(cluster)[1]
print(cluster.size)

FR.T0.mean <- mean(cluster$FR.T0)
FR.T3.mean <- mean(cluster$FR.T3)
FR.T4.mean <- mean(cluster$FR.T4)
BI.T0.mean <- mean(cluster$BI.T0)
BI.T3.mean <- mean(cluster$BI.T3)
BI.T4.mean <- mean(cluster$BI.T4)
mean.df <- data.frame("LFC"=as.numeric(c(FR.T0.mean, 
                                         FR.T3.mean, 
                                         FR.T4.mean,
                                         BI.T0.mean, 
                                         BI.T3.mean, 
                                         BI.T4.mean)),
                      "times"=c(1,2,3,5,6,7))


colors = c("darkcyan", "darkolivegreen3", "deeppink", "darkorange")


p <- ggplot()
for (ind in 1:200) { 
  dummy_df <- data.frame("LFC"=as.numeric(c(cluster[ind,1], 
                                            cluster[ind,2], 
                                            cluster[ind,3],
                                            cluster[ind,4],
                                            cluster[ind,5],
                                            cluster[ind,6])),
                         "times"=c(1,2,3,5,6,7))
  
  p <- p + geom_line(data=dummy_df, aes(x=times, y=LFC), color=colors[1])
}
p <- p + geom_line(data=mean.df, aes(x=times, y=LFC), color="red", size=1.5)

dev.new()
p + scale_x_continuous(name ="Cold Length (W)", 
                       breaks=c(1,2,3,5,6,7),
                       labels=c("FR-T0", "FR-T3", "FR-T4",
                                "BI-T0", "BI-T3", "BI-T4"),
                       limits=c(1,7)) +
  scale_y_continuous(name="Normalized LFC",
                     limits=c(0.4,.75)) +
  labs(title=paste("n = ", cluster.size))+
  theme(plot.title = element_text(size=18, hjust=0.5),
        text = element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA))



profiles["MNUTANS.r1.scaffold_52G00608510_MNUTANS.r1.scaffold_52G00608510",]



ind = which(rownames(profiles)=="MNUTANS.r1.scaffold_52G00608510_MNUTANS.r1.scaffold_52G00608510")
dummy_df <- data.frame("LFC"=as.numeric(c(profiles[ind,1], 
                                          profiles[ind,2], 
                                          profiles[ind,3],
                                          profiles[ind,4],
                                          profiles[ind,5],
                                          profiles[ind,6])),
                       "times"=c(1,2,3,5,6,7))

p <- ggplot() + geom_line(data=dummy_df, aes(x=times, y=LFC), color=colors[1], size=1.5)

dev.new()
p + scale_x_continuous(name ="Samples", 
                       breaks=c(1,2,3,5,6,7),
                       labels=c("FR-T0", "FR-T3", "FR-T4",
                                "BI-T0", "BI-T3", "BI-T4"),
                       limits=c(1,7)) +
  scale_y_continuous(name="Normalized Expression Level",
                     limits=c(0.45,.65)) +
  labs(title=paste("VRN1"))+ 
  theme(plot.title = element_text(size=18, hjust=0.5),
        text = element_text(size=14),
        panel.border = element_rect(colour = "black", fill=NA))
