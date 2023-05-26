"""This program allows to get the complete sequence from the indexes 
written in a text file. It will look into a fastq file which contains
all sequences and write the sequences corresponding to the indexes in
a fasta file. The output file can then be used for alignement.

Usage:
======
python sequence_insertion.py
"""

__authors__ = ("Jane Schadtler-Law")
__contact__ = ("jane.schadtler-law@etu.u-paris.fr")
__date__ = "2023-04-06"

with open("/data-isilon/FRED_ALL/FRED/PROJECT/TEST_NANOPORE/PROJECT12/BC_FASTQ/58_individu_6_220523_D4_U1_mProm_int_Cacnb1D_fw_U2_Ex13end_Cacnb1E_rv_BC30-U1_BC30-U2.fastq", "r") as fast, open("/home/jschadtler/Desktop/Proj_long/PROJET_cacnb1/insertion_list.txt", "r") as anomalies, open("/home/jschadtler/Desktop/Proj_long/PROJET_cacnb1/seq_insertion_list.fasta", "w") as output:
	lines_fast = fast.readlines()
	lines_ano = anomalies.readlines() 
	for ano in lines_ano:
		seq_ano = ano.split(" ")
		for i in range(0, len(lines_fast)):
			line = lines_fast[i]
			if line.startswith("@"+seq_ano[0]):
				output.write(">"+seq_ano[0]+'\n')
				output.write(lines_fast[i+1][int(seq_ano[1]):int(seq_ano[1])+int(seq_ano[2])]+'\n')
