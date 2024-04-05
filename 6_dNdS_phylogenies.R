#Plotting PAML outputs for cyathium paper
#Arielle Johnson 4/1/24
#modified from: https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/ggtree/inst/doc/treeImport.html#parsing-paml-output
#with additional advice and edits from Jacob Landis

library(ggtree)
library(treeio)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(makeunique)

#set working directory to where files are
setwd("~/Library/CloudStorage/Box-Box/Box Storage/Labwork/Cyathium paper/R plotting of PAML results")

#### LFY
mlc <- read.codeml_mlc("LFY_OG0007405_model1_NSSites0.txt")
mlc

raxml <- read.newick("LFY_OG0007405_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"LFY")
  }
}
species_label=make_unique(species_label, sep =" LFY_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "LFY_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

LFY <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 


#### AP3
mlc <- read.codeml_mlc("AP3_OG0001408_model1_NSSites0.txt")
mlc

raxml <- read.newick("AP3_OG0001408_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"AP3")
  }
}
species_label=make_unique(species_label, sep =" AP3_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "AP3_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

AP3 <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### PHYA

mlc <- read.codeml_mlc("PHYA_OG0007018_model1_NSSites0.txt")
mlc

raxml <- read.newick("PHYA_OG0007018_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"PHYA")
  }
}
species_label=make_unique(species_label, sep =" PHYA_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "PHYA_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

PHYA <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### AG

mlc <- read.codeml_mlc("AG_OG0008174_model1_NSSites0.txt")
mlc

raxml <- read.newick("AG_OG0008174_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"AG")
  }
}
species_label=make_unique(species_label, sep =" AG_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "AG_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

AG <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### PI

mlc <- read.codeml_mlc("PI_OG0004848_model1_NSSites0.txt")
mlc

raxml <- read.newick("PI_OG0004848_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"PI")
  }
}
species_label=make_unique(species_label, sep =" PI_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "PI_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

PI <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### UFO

mlc <- read.codeml_mlc("UFO_OG0008319_model1_NSSites0.txt")
mlc

raxml <- read.newick("UFO_OG0008319_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"UFO")
  }
}
species_label=make_unique(species_label, sep =" UFO_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "UFO_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

UFO <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### AP1-CAL

mlc <- read.codeml_mlc("AP1-CAL_OG0000496_model1_NSSites0.txt")
mlc

raxml <- read.newick("AP1-CAL_OG0000496_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT1G26310", species_label), "A. thaliana CAL")
species_label=replace(species_label, grep("AT1G69120", species_label), "A. thaliana AP1")
species_label=replace(species_label, grep("AT3G30260", species_label), "A. thaliana AGL79")
species_label=replace(species_label, grep("AT5G60910", species_label), "A. thaliana AGL8")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(grepl("A. thaliana",species_label[i], fixed = T)==T){
    next
  }
  else {
    if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"AP1")
    }
  }
}
species_label=make_unique(species_label, sep =" AP1_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "AP1-CAL_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

AP1_CAL <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 

#### SEP3

mlc <- read.codeml_mlc("SEP3_OG0003695_model1_NSSites0.txt")
mlc

raxml <- read.newick("SEP3_OG0003695_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"SEP3")
  }
}
species_label=make_unique(species_label, sep =" SEP3_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "SEP3_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

SEP3 <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 


#### SEP1

mlc <- read.codeml_mlc("SEP1_OG0000458_model1_NSSites0.txt")
mlc

raxml <- read.newick("SEP1_OG0000458_cds.new.gene_tree.raxml.support", node.label = "support")

my_merge <- merge_tree(mlc,raxml)

#label genes by species
gene_names <-as_tibble(raxml) %>% pull(label)
gene_names <-gene_names[!is.na(gene_names)]
species_label=replace(gene_names, grep("Casp", gene_names),"E. lathyris")
species_label=replace(species_label, grep("Ep", species_label), "E. peplus")
species_label=replace(species_label, grep("Manes", species_label), "M. esculenta")
species_label=replace(species_label, grep("Glyma", species_label), "G. max")
species_label=replace(species_label, grep("XP_", species_label), "R. communis")
species_label=replace(species_label, grep("AT", species_label), "A. thaliana")
species_label=replace(species_label, grep("Potri", species_label), "P. trichocarpa")
species_label=replace(species_label, grep("KAF", species_label), "H. brasiliensis")
for(i in 1:length(species_label)){
  if(!species_label[i] %in% species_label[duplicated(species_label)]){
    species_label[i] <- paste(species_label[i],"SEP1")
  }
}
species_label=make_unique(species_label, sep =" SEP1_", wrap_in_brackets = F)
d <-tibble(label=gene_names,species=species_label)
write_csv(d, file = "SEP1_orthologs.csv")

tree2 <- full_join(my_merge, d, by = "label")

SEP1 <- ggtree(tree2, branch.length = "t", aes(color=dN_vs_dS),size=1.5) + scale_color_continuous(name='dN/dS', limits=c(0, 0.9), oob=scales::squish, low="blue", high="red") +
  theme_tree2(legend.position=c(0.9, 0.8)) + geom_tiplab(aes(label=species), hjust=-0.03, size=2.7) + ggplot2::xlim(0, 3.9) + theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm")) + 
  geom_text(aes(label=support), color="black", vjust=-0.5, hjust=1.25, size=2.5) 


##Arrange trees into larger figures

ggarrange(PHYA, AG, AP3, PI, LFY, UFO,
          ncol = 2,
          nrow = 3, 
          common.legend = T, 
          labels = c("PHYA", "AG", "AP3", "PI", "LFY", "UFO")
) %>% ggsave(filename = "PAML_results.png", width=6, height=8.5,  units = "in")


ggarrange(SEP1, SEP3,
          ncol = 2,
          nrow = 1, 
          common.legend = T, 
          labels = c("SEP1", "SEP3")
) %>% ggsave(filename = "SEP_PAML_results.png", width=8.5, height=5,  units = "in")

ggsave(AP1_CAL, filename = "AP1_PAML_results.png", width=6, height=5,  units = "in")


