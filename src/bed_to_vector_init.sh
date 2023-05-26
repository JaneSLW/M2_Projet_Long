#!/bin/bash

# Script name 	: bed_to_vector_init.sh
# Description	: This program allows to build a binary vector for each read
# from the alignment results. It divides the bed file based on the length of
# the exons. Each file will be filtered by a threshold alignement score which
# defines the presence or absence of the exon. Then only the reads that start
# with either the first or third exon and end with the 20th or last exon are 
# kept.
#
# Author	: Jane Schadtler-Law
# Date	: 2023-03-27


REFNAME=${1}	#$1 reference
WDIR=${2}	# working directory
BEDold=${3}	# bed file
STARTS=${4} # example  (1 3)
STOPS=${5}	# example (20 21)

echo $REFNAME $WDIR $BEDold
FNAME=`basename ${REFNAME} | cut -d"." -f1`
mkdir -p ${WDIR}
BED=`basename ${BEDold}`

cp ${BEDold} ${WDIR}/${BED}

cd ${WDIR}

# Dividing bed file based on the length of the exon
# Example : If the exon is less than 25 bp then the match threshold is 60 %
awk '{if ( ($3-$2+1) <=25) {printf("%s\n",$0);}}' ${BED} > BED_600.bed 
awk '{if ((($3-$2+1) <=50) && (($3-$2+1)>25)) {printf("%s\n",$0);}}' ${BED} > BED_185.bed
awk '{if ((($3-$2+1) <=80) && (($3-$2+1)>50)) {printf("%s\n",$0);}}' ${BED} > BED_140.bed
awk '{if ((($3-$2+1) <=100) && (($3-$2+1)>80)) {printf("%s\n",$0);}}' ${BED} > BED_125.bed
awk '{if (($3-$2+1) >100) {printf("%s\n",$0);}}' ${BED} > BED_100.bed

for L in 600 185 140 125 100 ; do
bedtools intersect -f 0.${L} -wo -a BED_${L}.bed -b ../${REFNAME}
done > ${FNAME}.txt
sort -k8,8 -k2,2n ${FNAME}.txt > ${FNAME}_sorted.txt

ALL_Exons=`awk '{if (NR==1) printf("%s",$4); else printf(",%s",$4);}END{printf("\n");}' ${BED}`

cut -f4,8,10 ${FNAME}_sorted.txt | uniq | awk -v allexons=${ALL_Exons} 'BEGIN{n=split(allexons,a,","); }{  if ($2 != old){if (NR>1) {while (i<=n){printf("0"); i=i+1;} printf("\t%s\t",old);} i=1; for (k=1;k<j; k++){printf("%s",b[k]);} printf("\n"); old=$2;j=1;} while (i<=n){if (a[i]!=$1){printf("0");i=i+1;}else{ i=i+1; printf("1"); break;}  } b[j]=$3; j=j+1;}' | sort -k1,1nr > ${FNAME}_vector.txt
 
cat ${FNAME}_vector.txt | grep -v -E "\+\-|\-\+" | awk -v allexons=${ALL_Exons} 'BEGIN{n=split(allexons,a,","); }{nn=split($1,b,"");for(i=1; i<=n; i++){ if (b[i]=="1"){ if (a[i] ~ /ASS5/){ if (b[i-1]=="1"){printf("1");}else{printf("0");} }else{ if (a[i] ~ /ASS3/){ if (b[i+1]=="1"){printf("1"); }else{printf("0");} }else{ printf("1"); } } }else{printf("0");} } printf("\t%s\t%s\n",$2,$3); }' > ${FNAME}_vector_corrected.txt

cat ${FNAME}_vector_corrected.txt | cut -f1 | sort | uniq -c | sort -k1,1nr| awk '{printf("%s\t%s\n",$1,$2);}' > ${FNAME}_abundance_filter.txt 

for i in "${STARTS[@]}"
do
	# Start pattern
	NAME="Exon_${i}$"
	IND=`grep -n ${NAME} ${BED} | cut -d : -f1`
	pat_s=`pat_s=$(printf "%${IND}s");echo ${pat_s// /.}`
	
	if [ ${i} -gt 1 ];
	then
		for past_exon in `seq ${i}`
		do
			echo ${past_exon}
			IND_past_exon=`grep -n Exon_${past_exon}$ ${BED} | cut -d : -f1`
			pat_s=`echo ${pat_s} | sed s/./0/${IND_past_exon}`
		done
	fi

	pat_s=`echo ${pat_s} | sed s/./1/${IND}`

	# End pattern
	nb_l=`wc -l ${BED} | awk '{print $1}'`
	Ind=`grep -n Exon_${STOPS[0]}$ ${BED} | cut -d : -f1`
	len_pat=`expr ${nb_l} - ${Ind}`
	pat_end=`pat_end=$(printf "%${len_pat}s");echo ${pat_end// /0}`
	for i_last in "${STOPS[@]}"
	do
		echo ${i_last}
		ind=`grep -n Exon_${i_last}$ ${BED} | cut -d : -f1`
		ind=`expr ${ind} - ${Ind}`
		if [[ ${ind} -eq 0 ]]; then
			ind=1
			pat_end=`echo ${pat_end} | sed s/./1/${ind}`
		else
			pat_end=`echo ${pat_end} | sed s/./1/${ind}`
			before_ind=`expr ${ind} - 1`
			pat_end=`pat_deb=$(printf "%${before_ind}s");echo ${pat_deb// /.}`${pat_end[@]:${before_ind}} 
		fi
		
		if [[ `grep -n Exon_${i_last}_ASS5p ${BED} | cut -d : -f1` ]]; then
			ind_5p=`expr ${ind} + 1`
			pat_end=`echo ${pat_end} | sed s/./1/${ind_5p}`
		fi

		grep "^${pat_s}"  ${FNAME}_vector_corrected.txt | grep -v -E "\+\-|\-\+" | awk '{printf("%s\t%s\n",$2,$1)}' | grep "${pat_end}$" | cut -f2  | sort | uniq -c | sort -k1,1nr| awk -v allexons=${ALL_Exons} 'BEGIN{nexons=split(allexons,b,",");} {n=split($2,a,""); if (n==nexons){printf("%s\t%s\n",$1,$2);}}' > ${FNAME}_${i}_${i_last}_abundance_filter_11.txt 
	done
done
