# Arielle Johnson
# 4/16/24
# GO term enrichment for cyathium paper

library(topGO)
library(Rgraphviz)

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

Ep_nectary <- read_csv("nectary_res.csv") %>% 
  filter(log2FoldChange>0)

Ep_y_pistil <- read_csv("y_pistil_res.csv") %>% 
  filter(log2FoldChange>0)

Ep_y_stam <- read_csv("y_stam_res.csv") %>% 
  filter(log2FoldChange>0)

#make a list that is the GO terms and which genes they correspond to

my_GO <- read_tsv("Euphorbia_peplus_annot.txt", col_names=c("gene_ID", "GO_term"))
my_list <- my_GO$GO_term
names(my_list) <- my_GO$gene_ID

#get the genes of interest
stam_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_stam$gene_ID]
pistil_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_pistil$gene_ID]
filiform_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_filiform$gene_ID]
involucre_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_involucre$gene_ID]
nectary_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_nectary$gene_ID]
y_pistil_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_y_pistil$gene_ID]
y_stam_InterestingGenes <- my_GO$gene_ID[my_GO$gene_ID %in% Ep_y_stam$gene_ID]

#make a list of all the genes and show which ones are the genes of interest
#using 0's (not of interest) and 1's (of interest)
geneNames <- my_GO$gene_ID

stam_geneList<-factor(as.integer(geneNames%in%stam_InterestingGenes))
names(stam_geneList)<-geneNames
pistil_geneList<-factor(as.integer(geneNames%in%pistil_InterestingGenes))
names(pistil_geneList)<-geneNames
filiform_geneList<-factor(as.integer(geneNames%in%filiform_InterestingGenes))
names(filiform_geneList)<-geneNames
involucre_geneList<-factor(as.integer(geneNames%in%involucre_InterestingGenes))
names(involucre_geneList)<-geneNames
nectary_geneList<-factor(as.integer(geneNames%in%nectary_InterestingGenes))
names(nectary_geneList)<-geneNames
y_stam_geneList<-factor(as.integer(geneNames%in%y_stam_InterestingGenes))
names(y_stam_geneList)<-geneNames
y_pistil_geneList<-factor(as.integer(geneNames%in%y_pistil_InterestingGenes))
names(y_pistil_geneList)<-geneNames

#make a topGOdata dataset and specify subontology
#Biological Process (BP), Cellular Component (CC), Molecular Function (MF)
stam_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=stam_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
stam_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=stam_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
pistil_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=pistil_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
pistil_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=pistil_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
filiform_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=filiform_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
filiform_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=filiform_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
involucre_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=involucre_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
involucre_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=involucre_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
nectary_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=nectary_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
nectary_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=nectary_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
y_stam_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=y_stam_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
y_stam_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=y_stam_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
y_pistil_MF_GOdata<-new("topGOdata",ontology="MF",allGenes=y_pistil_geneList, annot=annFUN.gene2GO,gene2GO=my_list)
y_pistil_BP_GOdata<-new("topGOdata",ontology="BP",allGenes=y_pistil_geneList, annot=annFUN.gene2GO,gene2GO=my_list)


#run a Fisher's test 
test.stat<-new("classicCount",testStatistic=GOFisherTest,name="Fishertest")

stam_MF_resultFisher<-getSigGroups(stam_MF_GOdata,test.stat)
stam_BP_resultFisher<-getSigGroups(stam_BP_GOdata,test.stat)
pistil_MF_resultFisher<-getSigGroups(pistil_MF_GOdata,test.stat)
pistil_BP_resultFisher<-getSigGroups(pistil_BP_GOdata,test.stat)
filiform_MF_resultFisher<-getSigGroups(filiform_MF_GOdata,test.stat)
filiform_BP_resultFisher<-getSigGroups(filiform_BP_GOdata,test.stat)
involucre_MF_resultFisher<-getSigGroups(involucre_MF_GOdata,test.stat)
involucre_BP_resultFisher<-getSigGroups(involucre_BP_GOdata,test.stat)
nectary_MF_resultFisher<-getSigGroups(nectary_MF_GOdata,test.stat)
nectary_BP_resultFisher<-getSigGroups(nectary_BP_GOdata,test.stat)
y_stam_MF_resultFisher<-getSigGroups(y_stam_MF_GOdata,test.stat)
y_stam_BP_resultFisher<-getSigGroups(y_stam_BP_GOdata,test.stat)
y_pistil_MF_resultFisher<-getSigGroups(y_pistil_MF_GOdata,test.stat)
y_pistil_BP_resultFisher<-getSigGroups(y_pistil_BP_GOdata,test.stat)

#export it to view in Excel
#(The number of topNodes has to be less than the number of GO terms in the analysis)
stam_MF_GO<-GenTable(stam_MF_GOdata,classic=stam_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(stam_MF_GO, file="stam_MF_GO.csv")
stam_BP_GO<-GenTable(stam_BP_GOdata,classic=stam_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(stam_BP_GO, file="stam_BP_GO.csv")
pistil_MF_GO<-GenTable(pistil_MF_GOdata,classic=pistil_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(pistil_MF_GO, file="pistil_MF_GO.csv")
pistil_BP_GO<-GenTable(pistil_BP_GOdata,classic=pistil_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(pistil_BP_GO, file="pistil_BP_GO.csv")
filiform_MF_GO<-GenTable(filiform_MF_GOdata,classic=filiform_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(filiform_MF_GO, file="filiform_MF_GO.csv")
filiform_BP_GO<-GenTable(filiform_BP_GOdata,classic=filiform_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(filiform_BP_GO, file="filiform_BP_GO.csv")
involucre_MF_GO<-GenTable(involucre_MF_GOdata,classic=involucre_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(involucre_MF_GO, file="involucre_MF_GO.csv")
involucre_BP_GO<-GenTable(involucre_BP_GOdata,classic=involucre_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(involucre_BP_GO, file="involucre_BP_GO.csv")
nectary_MF_GO<-GenTable(nectary_MF_GOdata,classic=nectary_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(nectary_MF_GO, file="nectary_MF_GO.csv")
nectary_BP_GO<-GenTable(nectary_BP_GOdata,classic=nectary_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(nectary_BP_GO, file="nectary_BP_GO.csv")
y_stam_MF_GO<-GenTable(y_stam_MF_GOdata,classic=y_stam_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(y_stam_MF_GO, file="y_stam_MF_GO.csv")
y_stam_BP_GO<-GenTable(y_stam_BP_GOdata,classic=y_stam_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(y_stam_BP_GO, file="y_stam_BP_GO.csv")
y_pistil_MF_GO<-GenTable(y_pistil_MF_GOdata,classic=y_pistil_MF_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(y_pistil_MF_GO, file="y_pistil_MF_GO.csv")
y_pistil_BP_GO<-GenTable(y_pistil_BP_GOdata,classic=y_pistil_BP_resultFisher,topNodes=200) %>% 
  filter(classic<0.05)
write.csv(y_pistil_BP_GO, file="y_pistil_BP_GO.csv")

