# Arielle Johnson
# 4/16/24
# Comparing cyathium with flower

library(tximport)
library(rhdf5)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggpubr)
library(stringr)

#setting working directory to where files are 
directory <- "~/Library/CloudStorage/Box-Box/Box Storage/Labwork/Cyathium paper/RNAseq analysis/Working directory"
setwd(directory)

#load the results from Euphorbia peplus 
Ep_stam <- read_csv("m_stam_res.csv") %>% 
  filter(log2FoldChange>0)

Ep_pistil <- read_csv("m_pistil_res.csv") %>% 
  filter(log2FoldChange>0)

Ep_filiform <- read_csv("filiform_res.csv") %>% 
  filter(log2FoldChange>0)

Ep_involucre <- read_csv("involucre_res.csv") %>% 
  filter(log2FoldChange>0)


#making the tx2gene table for arabidopsis
tx2gene <- read.table("Arabidopsis_gene_names.txt", sep = "\t", header = F)
geneid <- list()
for (i in 1:length(tx2gene$V1)){
  geneid[i] <- unlist(strsplit(tx2gene$V1[i], split = "[.]"))[1]
}
tx2gene$geneid <- geneid
names(tx2gene)[names(tx2gene) == 'V1'] <- "transcriptid"
class(tx2gene$geneid) <- "character"
#head(tx2gene)
write.table(tx2gene, file='Arabidopsis_tx2gene.tsv', quote=F, sep='\t', col.names = T, row.names=F)
  

#do DESeq2 on Arabidopsis data 

#import sample information file (file, sample, condition)
my_samples <- read.table(file = 'arabidopsis_samples.csv', sep = ',', header = TRUE)
#head(my_samples)

#make the "file" column into a vector called "files"
files <- my_samples[,"file"]

#set the names of the files to be the sample names
names(files) <- my_samples[,c("sample")]

#make the "sample" and "condition" columns into a new dataframe, "sampleTable"
sampleTable <- my_samples[,c("sample","condition")]

#importing counts with tximport
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")
#head(txi$counts)

#Create a DESeq dataset from the txi dataset
ddsTxi <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
#it converts characters to factors, that is fine

# get genes with summed counts greater than 20
sumcounts <- rowSums(counts(ddsTxi))
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsTxi <- ddsTxi[keep,]

#run the DESeq pipeline for genes with counts greater than 20 only
dds <- DESeq(ddsTxi)

#variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

#plot the PCA
plotPCA(vsd, intgroup="condition")

# arabidopsis carpel DE
carpel_res <- results(dds, alpha=0.05,  contrast=c("condition","carpels","leaf"))
summary(carpel_res)
carpel_res <- carpel_res[order(carpel_res$pvalue),]
At_carpel <- subset(carpel_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID") %>% 
  filter(log2FoldChange>0)

# arabidopsis petal DE
petal_res <- results(dds, alpha=0.05,  contrast=c("condition","petals","leaf"))
summary(petal_res)
petal_res <- petal_res[order(petal_res$pvalue),]
At_petal <- subset(petal_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID") %>% 
  filter(log2FoldChange>0)

# arabidopsis filament DE
filament_res <- results(dds, alpha=0.05,  contrast=c("condition","filaments","leaf"))
summary(filament_res)
filament_res <- filament_res[order(filament_res$pvalue),]
At_filament <- subset(filament_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID") %>% 
  filter(log2FoldChange>0)


# arabidopsis anther DE
anther_res <- results(dds, alpha=0.05,  contrast=c("condition","anthers","leaf"))
summary(anther_res)
anther_res <- anther_res[order(anther_res$pvalue),]
At_anther <- subset(anther_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID") %>% 
  filter(log2FoldChange>0)


# arabidopsis sepal DE
sepal_res <- results(dds, alpha=0.05,  contrast=c("condition","sepals","leaf"))
summary(sepal_res)
sepal_res <- sepal_res[order(sepal_res$pvalue),]
At_sepal <- subset(sepal_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID") %>% 
  filter(log2FoldChange>0)


#load the orthogroups
orthogroups <- read.table("Orthogroups_only_At_and_Ep.txt", sep = "\t", header = T)
#tail(orthogroups)

#get rid of the empty orthogroups
orthogroups <- orthogroups[!(orthogroups$E_peplus=="" & orthogroups$A_thaliana==""), ]
#tail(orthogroups)

orthogroups %>% 
  filter(str_detect(E_peplus,"Ep_chr6_g16520"))

#make new dataframe to collect results
mydata <- data.frame("organ"=character(), "At_upreg"=numeric(), "At_not_upreg"=numeric(), "No_At_gene"=numeric())

#E. peplus pistillate flower vs. Arabidopsis carpel
Ep_pistil_orthos <- orthogroups %>% 
  filter(str_detect(E_peplus, paste(Ep_pistil$gene_ID, collapse = "|")))
#total # of orthogroups with Ep genes upreg in pistil
Ep_orthogroup_total <- nrow(Ep_pistil_orthos)
#the # of those that are upreg in Arabidopsis thaliana carpel
At_and_Ep_upreg <- Ep_pistil_orthos %>% 
  filter(str_detect(A_thaliana, paste(At_carpel$gene_ID, collapse = "|"))) %>% 
  nrow()
#the ones where there's not an At gene in the orthogroup with the Ep gene
No_At_gene <- nrow(Ep_pistil_orthos[(Ep_pistil_orthos$A_thaliana==""),])
#the ones where At gene exists but is not upreg
At_not_upreg <- Ep_orthogroup_total-At_and_Ep_upreg-No_At_gene
newrow <-  list("carpel", At_and_Ep_upreg, At_not_upreg, No_At_gene)
mydata <- rbind(mydata, newrow)


#E. peplus filiform structure vs. Arabidopsis petal
Ep_filiform_orthos <- orthogroups %>% 
  filter(str_detect(E_peplus, paste(Ep_filiform$gene_ID, collapse = "|")))
#total # of orthogroups with Ep genes upreg in filiform str.
Ep_orthogroup_total <- nrow(Ep_filiform_orthos)
#the # of those that are upreg in Arabidopsis thaliana petal
At_and_Ep_upreg <- Ep_filiform_orthos %>% 
  filter(str_detect(A_thaliana, paste(At_petal$gene_ID, collapse = "|"))) %>% 
  nrow()
#the ones where there's not an At gene in the orthogroup with the Ep gene
No_At_gene <- nrow(Ep_filiform_orthos[(Ep_filiform_orthos$A_thaliana==""),])
#the ones where At gene exists but is not upreg
At_not_upreg <- Ep_orthogroup_total-At_and_Ep_upreg-No_At_gene
newrow <-  list("petal", At_and_Ep_upreg, At_not_upreg, No_At_gene)
mydata <- rbind(mydata, newrow)

#E. peplus involucre vs. Arabidopsis sepal
Ep_involucre_orthos <- orthogroups %>% 
  filter(str_detect(E_peplus, paste(Ep_involucre$gene_ID, collapse = "|")))
#total # of orthogroups with Ep genes upreg in involucre
Ep_orthogroup_total <- nrow(Ep_involucre_orthos)
#the # of those that are upreg in Arabidopsis thaliana sepal
At_and_Ep_upreg <- Ep_involucre_orthos %>% 
  filter(str_detect(A_thaliana, paste(At_sepal$gene_ID, collapse = "|"))) %>% 
  nrow()
#the ones where there's not an At gene in the orthogroup with the Ep gene
No_At_gene <- nrow(Ep_involucre_orthos[(Ep_involucre_orthos$A_thaliana==""),])
#the ones where At gene exists but is not upreg
At_not_upreg <- Ep_orthogroup_total-At_and_Ep_upreg-No_At_gene
newrow <-  list("sepal", At_and_Ep_upreg, At_not_upreg, No_At_gene)
mydata <- rbind(mydata, newrow)

#E. peplus staminate flower vs. Arabidopsis filament
Ep_stam_orthos <- orthogroups %>% 
  filter(str_detect(E_peplus, paste(Ep_stam$gene_ID, collapse = "|")))
#total # of orthogroups with Ep genes upreg in staminate flower
Ep_orthogroup_total <- nrow(Ep_stam_orthos)
#the # of those that are upreg in Arabidopsis thaliana filament or stamen
stamen_and_filament <- c(At_filament$gene_ID, At_anther$gene_ID)
At_and_Ep_upreg <- Ep_stam_orthos %>% 
  filter(str_detect(A_thaliana, paste(stamen_and_filament, collapse = "|"))) %>% 
  nrow()
#the ones where there's not an At gene in the orthogroup with the Ep gene
No_At_gene <- nrow(Ep_stam_orthos[(Ep_stam_orthos$A_thaliana==""),])
#the ones where At gene exists but is not upreg
At_not_upreg <- Ep_orthogroup_total-At_and_Ep_upreg-No_At_gene
newrow <-  list("stamen", At_and_Ep_upreg, At_not_upreg, No_At_gene)
mydata <- rbind(mydata, newrow)



colnames(mydata)= c("organ", "At_upreg", "At_not_upreg", "No_At_gene")

mydata

pdata <- pivot_longer(mydata, c("At_upreg","At_not_upreg", "No_At_gene"), 
                      names_to = "arabidopsis", values_to = "num_orthogroups") %>%
  mutate(organ = factor(organ, levels=c("sepal", "petal", "stamen", "carpel"))) %>%
  mutate(arabidopsis = factor(arabidopsis, levels=c("No_At_gene", "At_not_upreg", "At_upreg")))

ggplot(pdata, aes(fill=arabidopsis, y=num_orthogroups, x=organ)) + 
  geom_bar(position="stack", stat="identity")

ggplot(pdata, aes(fill=arabidopsis, y=num_orthogroups, x=organ)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  scale_fill_manual(values=c("#e69f00", "#f0e442", "#009e73"), 
                    labels = c("no arabidopsis gene", 
                               "arabidopsis gene not upregulated", 
                               "arabidopsis gene upregulated")) +
  ylab("Proportion of orthogroups") +
  theme(axis.title.x=element_blank(), legend.title=element_blank()) 
ggsave("stacked_barplot_whorl_orthogroups.jpg", width=7, height=3, units = "in")

