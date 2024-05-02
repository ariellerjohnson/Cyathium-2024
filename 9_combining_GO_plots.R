#Arielle Johnson
#4/16/24
#Combining GO plots for cyathium paper

#ran this after running the four files for the individual figures

library(ggpubr)

ggarrange(nectary_BP, nectary_MF,involucre_BP, filiform_BP,
          ncol = 1,
          nrow = 4, 
          common.legend = T, 
          labels = c("Nectary Biological Process Enrichment", 
                     "Nectary Molecular Function Enrichment",
                     "Involucre Biological Process Enrichment", 
                     "Filiform Structure Biological Process Enrichment"), 
          vjust = 1.8, 
          hjust = -0.15
) %>% ggsave(filename = "Figure XX-- GO analysis.png", width=9, height=11,  units = "in")
