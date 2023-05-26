# Author : Jane Schadtler-Law
# Contact : jane.schadtler-law@etu.u-paris.fr
# Date : April 2023
#
# This program allow to extract indexes of inserted segment of sequences and 
# indexes of anomalies. Anomalies are sequences which have a segment on the 
# other strand (of opposite sign). The results of the alignment in .maf format 
# have been previously converted to .tab format.
#
# Usage Rscript index_insertion_anomaly.R file_with_indexes.tab
# Example : Rscript index_insertion_anomaly.R res_alignmentind_6_mProm_E13.tab


# Retrieve arguments 
Args <- commandArgs(trailingOnly = TRUE)
if(length(Args)<1){
  stop("Missing file")
}

tab <- read.table(Args[1], sep="\t") # results of the alignment in .tab format

# Initializing column names for the future tab
to_keep <- c("Sequence", "Start", "End")
anomalie <- c("Sequence", "Start", "Length", "Normal strand", "Actual strand")

# Loop to categorize each sequence : anomaly, insertion or normal
for(seq in unique(tab[,7])){
  strand <- unique(tab[which(tab[,7] == seq),10]) # get sign of the strand 
  # Anomalies
  if(length(strand)>1){ 
    if(sum(tab[which(tab[,7] == seq),10] == "+") > sum(tab[which(tab[,7] == seq),10] == "-")){
      strand <- "+"
      ano <- "-"
      } else {
      strand <- "-"
      ano <- "+"
      }
    for(a in  which(tab[which(tab[,7] == seq),10] == ano)){
      anomalie <- rbind(anomalie, c(seq, tab[a,8], tab[a,9], strand, ano))
      }
    }
  # Insertions
  if(strand == "+"){ 
    for(i in head(which(tab[,7] == seq),-1)){
      if(tab[i+1,8] - (tab[i,8]+tab[i,9]) > 10){ 
        to_add <- c(seq, tab[i,8]+tab[i,9], tab[i+1,8]-1)
        to_keep <- rbind(to_keep,to_add)
      }
    }
  }
  else if(strand == "-"){
    for(i in head(rev(which(tab[,7] == seq)),-1)){
      if(tab[i-1,8] - (tab[i,8]+tab[i,9]) > 10){ 
        to_add <- c(seq, tab[i,8]+tab[i,9], tab[i-1,8]-1)
        to_keep <- rbind(to_keep,to_add)
      }
    }
  }
}

# Write in .txt the indexes of insertions and anomalies in separate files
write.table(to_keep[-1,], file = "./insertion_list.txt",row.names = F, col.names = F, quote=F)
write.table(anomalie, file = "./anomalie_list.txt",row.names = F, col.names = F, quote=F)
