# This program:
# 1 - reads in a batch_1.sumstats.tsv file for a  run of populatiosn with a given set of paramaters
# 2 - takes the first SNP in a contig and writes it to an output file to be used as a Whitelist
# 3 - prints out the contig that is associated with each used SNP to a separate file (as a check that I've done things correctly)
# 4 - adds up all of the loci that are NOT used, which, again, serves as a check that I've done things properly

# ****NOTE**** MUST check the number of lines in the header of the sumstats file before running this
# ****NOTE**** MUST to update the names of the output files approproiately before running this program

#The column values for the batch_1.sumstats.tsv file:
Batch_ID_Sumstats = 0
Locus_ID_Sumstats = 1
Chromosome_Sumstats = 2
Basepair_Sumstats = 3
Column_Sumstats = 4
Population_ID_Sumstats = 5
P_Nucleotide_Sumstats = 6
Q_Nucleotide_Sumstats = 7
Number_of_Individuals_Sumstats = 8
P_Sumstats = 9
Observed_Heterozygosity_Sumstats = 10
Observed_Homozygosity_Sumstats = 11
Expected_Heterozygosity_Sumstats = 12
Expected_Homozygosity_Sumstats = 13
Pi_Sumstats = 14
Smoothed_Pi_Sumstats = 15
Smoothed_Pi_P_value_Sumstats = 16
FIS_Sumstats = 17
Smoothed_FIS_Sumstats = 18
Smoothed_FIS_P_value_Sumstats = 19
Private_allele_Sumstats = 20



#Open STACKS sumstats file with all of the F_is information:
Stats_file =  open('/PATH_TO_DATA/batch_1.sumstats.tsv')
#Open an optout file for the F_is info:
output_Whitelist =  open('/PATH_TO_DATA/Whitelist_SNPsFromTrueSamples.txt', 'wt')



#Set counter to skip lines:
i=0

for data in Stats_file:	
	if i<3:
		i=i+1
		dataLinePrev = data.rstrip('\n').split('\t')
	else:
		dataLine = data.rstrip('\n').split('\t')
		if dataLine[Locus_ID_Sumstats] != dataLinePrev[Locus_ID_Sumstats] or dataLine[Column_Sumstats] != dataLinePrev[Column_Sumstats]:
			output_Whitelist.write(str(dataLine[Locus_ID_Sumstats]) + '\t' + str(dataLine[Column_Sumstats]) + '\n')
		dataLinePrev = dataLine

Stats_file.close()
output_Whitelist.close()




