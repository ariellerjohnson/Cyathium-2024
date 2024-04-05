#RNAseq analysis for cyathium paper 
#Arielle Johnson, 4/1/23

#load packages 
library(tximport)
library(rhdf5)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ggpubr)

#setting working directory to where files are 
#files: kallisto output folders, tx2gene file, sample info file
directory <- "~/Library/CloudStorage/Box-Box/Box Storage/Labwork/Cyathium paper/RNAseq analysis/Working directory"
setwd(directory)

#checking that the right files & directories are present
list.files(directory)

#import tx2gene file (shows which transcripts correspond to which genes)
tx2gene <- read.table(file = 'E_peplus_tx2gene.tsv', sep = '\t', header = TRUE)
#head(tx2gene)

#import sample information file (file, sample, condition)
my_samples <- read.table(file = 'my_samples.tsv', sep = '\t', header = TRUE)
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

#variance stabilizing transformation for visualization
vsd <- vst(dds, blind=TRUE)

#plot the PCA
#plotPCA(vsd, intgroup="condition")
dat <- plotPCA(vsd, intgroup="condition", returnData=TRUE) %>% 
  mutate(group = factor(group, levels=c("cyathophyll", "nectary", "involucre", "filiform", "m. stam.", "m. pistil.", "y. pistil.","y. stam.", "lg. primord.", "sm. primord.")))

#shape version
p <- ggplot(dat,aes(x=PC1,y=PC2,
                    shape=group, 
                    col=group))
p <- p + geom_point(size=3) + scale_shape_manual(values=1:nlevels(dat$group)) +
  theme_minimal()
p

#decided to add PC3 and PC4 inspired by the following vignette:
##https://www.biostars.org/p/333436/
rv <- rowVars(assay(vsd))
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
#summary(pca)
PC1=pca$x[,1]
PC2=pca$x[,2]
#yes these are the same data as in "dat", good
PC3=pca$x[,3]
PC4=pca$x[,4]

plot(PC3, PC4)

moar_dat <- PC3 %>%
  enframe(value="PC3") %>% 
  left_join(dat, by="name")

super_dat <- PC4 %>%
  enframe(value="PC4") %>% 
  left_join(moar_dat, by="name")

p2 <- ggplot(super_dat,aes(x=PC3,y=PC4,
                           shape=group, 
                           col=group))
p2 <- p2 + geom_point(size=3) + scale_shape_manual(values=1:nlevels(dat$group)) +
  theme_minimal()
p2

#plotting PC1/PC2 and PC3/PC4 together
ggarrange(p, p2, common.legend=TRUE, legend = "left")
ggsave("PC_1_thru_4.tiff", width=6, height=3, units="in")
ggsave("PC_1_thru_4.png", width=6, height=3, units="in")

#get genes with top variation across all conditions for heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=T),100)

#plot the heatmap
pheatmap(
  assay(vsd)[topVarGenes,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE, 
  scale = "row"
)

#check the annotations of the top var genes
ahrd_annos <- read.table(file = 'ahrd_output_sorted.txt', sep = '\t', header = TRUE, skip=1)
ahrd_annos$gene_ID <- gsub("\\.t\\d$", "", ahrd_annos$Protein.Accession)
my_annos <- ahrd_annos[!duplicated(ahrd_annos$gene_ID), c("Human.Readable.Description", "gene_ID")]
my_topVar <- as.data.frame(assay(vsd)[topVarGenes,])
topVar_annos <- my_annos[my_annos$gene_ID %in% rownames(my_topVar),]

#looking at individual gene
gene_counts <- plotCounts(dds, gene = "Ep_chr1_g00024",
                          intgroup="condition", returnData = T)

#plot dots connected by lines for individual gene expression
ggplot(gene_counts, aes(x = condition, y = count, color = condition, group=1)) +
  scale_y_log10() +  geom_point() + 
  stat_summary(geom="line", fun="mean", size=0.5, color="black", linewidth = 1.2)

##based on BLAST results, Ep_chr1_g00024 appears to be homologous
##to ABORTED MICROSPORES from Arabidopsis which controls
##pollen wall formation... good gut check 

#heatmap of floral identity genes of interest
genes_of_interest <- read.table(file = 'Floral_genes_of_interest.txt', sep = '\t', header = TRUE)

#getting averaged values across organs for heatmap
my_vsd <- as.data.frame(assay(vsd))

vsd_averaged <- data.frame(row.names=rownames(my_vsd))
vsd_averaged$"cyathophyll"<- rowMeans(my_vsd[, grep("_B\\d$", colnames(my_vsd))])
vsd_averaged$"nectary"<- rowMeans(my_vsd[, grep("_N\\d$", colnames(my_vsd))])
vsd_averaged$"involucre"<- rowMeans(my_vsd[, grep("_S\\d$", colnames(my_vsd))])
vsd_averaged$"filiform"<- rowMeans(my_vsd[, grep("_P\\d$", colnames(my_vsd))])
vsd_averaged$"m. stam."<- rowMeans(my_vsd[, grep("_A\\d$", colnames(my_vsd))])
vsd_averaged$"m. pistil."<- rowMeans(my_vsd[, grep("_G\\d$", colnames(my_vsd))])
vsd_averaged$"y. pistil."<- rowMeans(my_vsd[, grep("_TG\\d$", colnames(my_vsd))])
vsd_averaged$"y. stam."<- rowMeans(my_vsd[, grep("_TA\\d$", colnames(my_vsd))])
vsd_averaged$"lg. primord."<- rowMeans(my_vsd[, grep("_LP\\d$", colnames(my_vsd))])
vsd_averaged$"sm. primord."<- rowMeans(my_vsd[, grep("_SP\\d$", colnames(my_vsd))])


#Plotting the floral genes of interest heatmap
myvsd <- vsd_averaged %>% 
  rownames_to_column(var="gene_ID")

floral_genes <- read.table(file = 'Floral_genes_of_interest.txt', sep = '\t', header = TRUE)

floral_vsd <- inner_join(myvsd, floral_genes, by="gene_ID") %>% 
  column_to_rownames(var="gene_name") %>% 
  arrange(class)

##the good one for the figure, scaled by row/gene
pheatmap(
  select(floral_vsd,-c("gene_ID", "class")), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE, 
  scale = "row", 
  annotation_row = select(floral_vsd, "class"), 
  annotation_colors = list(class=c(A="#CC79A7", B="#E69F00", C="#F0E442", D="#009E73", E="#0072B2", other="#000000")), 
  angle_col = 45
)
#saved as tiff: 550 x 500

#not scaled by row (can somewhat see relative levels of expression between genes)
pheatmap(
  floral_vsd, 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE, 
  #scale = "row"
)

##The floral genes not included because the counts were too low
floral_genes[!(floral_genes$gene_ID %in% rownames(vsd_averaged)),]



##COMPARING CONDITIONS: each organ vs. cyathophyll

#nectary vs. cyathophyll
nectary_res <- results(dds, alpha=0.05,  contrast=c("condition","nectary","cyathophyll"))
summary(nectary_res)
nectary_res <- nectary_res[order(nectary_res$pvalue),]
nectary_res_df <- subset(nectary_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
nectary_res_df_annos <- left_join(nectary_res_df, my_annos, by= "gene_ID")
write_csv(nectary_res_df_annos, file = "nectary_res.csv")

#involucre vs. cyathophyll
involucre_res <- results(dds, alpha=0.05,  contrast=c("condition","involucre","cyathophyll"))
summary(involucre_res)
involucre_res <- involucre_res[order(involucre_res$pvalue),]
involucre_res_df <- subset(involucre_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
involucre_res_df_annos <- left_join(involucre_res_df, my_annos, by= "gene_ID")
write_csv(involucre_res_df_annos, file = "involucre_res.csv")


#filiform vs. cyathophyll
filiform_res <- results(dds, alpha=0.05,  contrast=c("condition","filiform","cyathophyll"))
summary(filiform_res)
filiform_res <- filiform_res[order(filiform_res$pvalue),]
filiform_res_df <- subset(filiform_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
filiform_res_df_annos <- left_join(filiform_res_df, my_annos, by= "gene_ID")
write_csv(filiform_res_df_annos, file = "filiform_res.csv")


#mature staminate vs. cyathophyll
m_stam_res <- results(dds, alpha=0.05,  contrast=c("condition","m. stam.","cyathophyll"))
summary(m_stam_res)
m_stam_res <- m_stam_res[order(m_stam_res$pvalue),]
m_stam_res_df <- subset(m_stam_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
m_stam_res_df_annos <- left_join(m_stam_res_df, my_annos, by= "gene_ID")
write_csv(m_stam_res_df_annos, file = "m_stam_res.csv")

#mature pistillate vs. cyathophyll
m_pistil_res <- results(dds, alpha=0.05,  contrast=c("condition","m. pistil.","cyathophyll"))
summary(m_pistil_res)
m_pistil_res <- m_pistil_res[order(m_pistil_res$pvalue),]
m_pistil_res_df <- subset(m_pistil_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
m_pistil_res_df_annos <- left_join(m_pistil_res_df, my_annos, by= "gene_ID")
write_csv(m_pistil_res_df_annos, file = "m_pistil_res.csv")


#young pistillate vs. cyathophyll
y_pistil_res <- results(dds, alpha=0.05,  contrast=c("condition","y. pistil.","cyathophyll"))
summary(y_pistil_res)
y_pistil_res <- y_pistil_res[order(y_pistil_res$pvalue),]
y_pistil_res_df <- subset(y_pistil_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
y_pistil_res_df_annos <- left_join(y_pistil_res_df, my_annos, by= "gene_ID")
write_csv(y_pistil_res_df_annos, file = "y_pistil_res.csv")


#young staminate vs. cyathophyll
y_stam_res <- results(dds, alpha=0.05,  contrast=c("condition","y. stam.","cyathophyll"))
summary(y_stam_res)
y_stam_res <- y_stam_res[order(y_stam_res$pvalue),]
y_stam_res_df <- subset(y_stam_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
y_stam_res_df_annos <- left_join(y_stam_res_df, my_annos, by= "gene_ID")
write_csv(y_stam_res_df_annos, file = "y_stam_res.csv")

#large primordium vs. cyathophyll
lg_primord_res <- results(dds, alpha=0.05,  contrast=c("condition","lg. primord.","cyathophyll"))
summary(lg_primord_res)
lg_primord_res <- lg_primord_res[order(lg_primord_res$pvalue),]
lg_primord_res_df <- subset(lg_primord_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
lg_primord_res_df_annos <- left_join(lg_primord_res_df, my_annos, by= "gene_ID")
write_csv(lg_primord_res_df_annos, file = "lg_primord_res.csv")

#small primordium vs. cyathophyll
sm_primord_res <- results(dds, alpha=0.05,  contrast=c("condition","sm. primord.","cyathophyll"))
summary(sm_primord_res)
sm_primord_res <- sm_primord_res[order(sm_primord_res$pvalue),]
sm_primord_res_df <- subset(sm_primord_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
sm_primord_res_df_annos <- left_join(sm_primord_res_df, my_annos, by= "gene_ID")
write_csv(sm_primord_res_df_annos, file = "sm_primord_res.csv")


#mature pistillate vs young pistillate
p_age_res <- results(dds, alpha=0.05,  contrast=c("condition","m. pistil.","y. pistil."))
summary(p_age_res)
p_age_res <- p_age_res[order(p_age_res$pvalue),]
p_age_res_df <- subset(p_age_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
p_age_res_df_annos <- left_join(p_age_res_df, my_annos, by= "gene_ID")
write_csv(p_age_res_df_annos, file = "p_age_res.csv")

#mature staminate vs young staminate
s_age_res <- results(dds, alpha=0.05,  contrast=c("condition","m. stam.","y. stam."))
summary(s_age_res)
s_age_res <- s_age_res[order(s_age_res$pvalue),]
s_age_res_df <- subset(s_age_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
s_age_res_df_annos <- left_join(s_age_res_df, my_annos, by= "gene_ID")
write_csv(s_age_res_df_annos, file = "s_age_res.csv")

#mature staminate vs mature pistillate
mature_res <- results(dds, alpha=0.05,  contrast=c("condition","m. stam.","m. pistil."))
summary(mature_res)
mature_res <- mature_res[order(mature_res$pvalue),]
mature_res_df <- subset(mature_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
mature_res_df_annos <- left_join(mature_res_df, my_annos, by= "gene_ID")
write_csv(mature_res_df_annos, file = "m_stam_vs_pistil_res.csv")

#young staminate vs young pistillate
young_res <- results(dds, alpha=0.05,  contrast=c("condition","y. stam.","y. pistil."))
summary(young_res)
young_res <- young_res[order(young_res$pvalue),]
young_res_df <- subset(young_res, padj < 0.01) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="gene_ID")
young_res_df_annos <- left_join(young_res_df, my_annos, by= "gene_ID")
write_csv(young_res_df_annos, file = "y_stam_vs_pistil_res.csv")

#copied and pasted all of these csvs into the same Excel file for supplementary info
