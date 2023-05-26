# Author : Jane Schadtler-Law
# Contact : jane.schadtler-law@etu.u-paris.fr
# Date : April 2023
#
# This program allows to calculates proportions of isoforms and creates stacked 
# barplot. The input data is list of isoforms fin each samples. This list has 
# been filtered, if one RNA isoform has at least 10 copies in at least one 
# sample, the isoform is kept. For the barplot representation, the isoforms 
# with a proportion of less than 0.5% are filtered.

# Libraries
library(RColorBrewer)
library(ggplot2)

# Barplot
res = read.table("/Users/janelaw/Desktop/Uni/M2/Projet_Long/CACNB1_matrice_key10_VS_BC.csv", header = T, sep=",", numerals ="no.loss", check.names = F)
row.names(res) <- as.character(res[,1])
res <- res[,-1]

# Get proportion of transcripts
res_prop <- apply(res,2,function(x){x/sum(x)})

# Separe by amplicon
res_prop_amorces <-res_prop[,grepl("U1_E1_Cacnb1E_fwrd_U1_mProm_int_Cacnb1D_fw_U2_Ex14_Cacnb1E_rv_U2_Ex13end_Cacnb1E_rv_", colnames(res))]  # Toutes les amorces 
res_prop_E1_E14 <-res_prop[,grepl("U1_E1_Cacnb1E_fwrd_U2_Ex14_Cacnb1E_rv_BC", colnames(res))] # E1 - E14
res_prop_E1_E13 <-res_prop[,grepl("U1_E1_Cacnb1E_fwrd_U2_Ex13end_Cacnb1E_rv", colnames(res))] # E1 - Ex 13
res_prop_PromInt_E14 <-res_prop[,grepl("U1_mProm_int_Cacnb1D_fw_U2_Ex14_Cacnb1E_rv_BC", colnames(res))] # Prom_int - Ex 14
res_prop_PromInt_E13 <-res_prop[,grepl("U1_mProm_int_Cacnb1D_fw_U2_Ex13end_Cacnb1E_rv_BC", colnames(res))]  # Prom_int - Ex 13

# get rid of transcripts which aren't recovered by sequences
res_prop_amorces_fil <- res_prop_amorces[rowSums(res_prop_amorces[])>0, ]
res_prop_E1_E14_fil <- res_prop_E1_E14[rowSums(res_prop_E1_E14[])>0, ]
res_prop_E1_E13_fil <- res_prop_E1_E13[rowSums(res_prop_E1_E13[])>0, ]
res_prop_PromInt_E14_fil <- res_prop_PromInt_E14[rowSums(res_prop_PromInt_E14[])>0, ]
res_prop_PromInt_E13_fil <- res_prop_PromInt_E13[rowSums(res_prop_PromInt_E13[])>0, ]

# rank 
prop_amorces_ranked <- res_prop_amorces_fil[order(apply(res_prop_amorces_fil,1,sum), decreasing = T), ]
prop_E1_E14_ranked <- res_prop_E1_E14_fil[order(apply(res_prop_E1_E14_fil,1,sum), decreasing = T), ]
prop_E1_E13_ranked <- res_prop_E1_E13_fil[order(apply(res_prop_E1_E13_fil,1,sum), decreasing = T), ]
prop_PromInt_E14_ranked <- res_prop_PromInt_E14_fil[order(apply(res_prop_PromInt_E14_fil,1,sum), decreasing = T), ]
prop_PromInt_E13_ranked <- res_prop_PromInt_E13_fil[order(apply(res_prop_PromInt_E13_fil,1,sum), decreasing = T), ]

# Stacked Bar plot
titles <- c("Les 4 amorces", "Amorces E1-E14", "Amorces E1-E13", "Amorces PromInt-E14", "Amorces PromInt-E13")
stades <- c("Embryon stade précoce", "Embryon stade avancée", "Nourisson", "Adulte dose moyenne 2M", "Adulte dose elevée 2M", "Adulte dose elevée 3M")
cont <- 1

# Loop start
for(dataset in list(prop_amorces_ranked, prop_E1_E14_ranked, prop_E1_E13_ranked, prop_PromInt_E14_ranked,prop_PromInt_E13_ranked)){
  dataset_fil <- dataset[rowSums(dataset[])>0.005, ]
  nb_isofrom <- dim(dataset_fil)[1]
  ind <- rep(stades, each = nb_isofrom)
  transcripts <- rep(row.names(dataset_fil) , 6)
  val <- NULL
  for(i in 1:length(colnames(dataset_fil))){
    val <- c(val,dataset_fil[,i] )
  }
  dataset_formated <- data.frame(ind,transcripts,val)
  
  mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb_isofrom)

  ggplot(dataset_formated, aes(fill=transcripts, y=val, x=factor(ind, stades))) + 
    geom_bar(position="stack", stat="identity") + 
    theme_bw()+
    scale_fill_manual(values = mycolors) + 
    theme(axis.text.x = element_text(angle = 60, vjust = 0.8, hjust=0.8))+
    labs(title=titles[cont],
         x = "Échantillons",
         y = "Proportions des isoformes",
         fill = "Détail des isoformes" )+
    theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
          legend.key.height = unit(0.3, 'cm'), #change legend key height
          legend.key.width = unit(0.3, 'cm'), #change legend key width
          legend.title = element_text(size=8), #change legend title font size
          legend.text = element_text(size=5)) #change legend text font size
 
 ggsave(
   filename = paste0("/Users/janelaw/Desktop/Uni/M2/Projet_Long/Graph/Barplot_",cont,".pdf"),
   plot = last_plot(),
   scale = 1,
   width = 8,
   height = 6,
   units = "in",
 )
 
 cont <- cont +1
} # Loop end
