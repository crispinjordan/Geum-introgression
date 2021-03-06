#This is a very simple file that reads through a batch_1.sumstats.tsv file and simply lists every catalog locus that is present to create a white-list.
#This is useful if I want to use the same set of loci in different analyses - i.e., take the loci that were used in 1 analysis with this script, and pass them
#in a whitelist for another STACKS analysis.

#Line counter for reading file:
i=0

#Number the columns:
BatchID = 0
LocusID = 1
Chr = 2
BP = 3
Col = 4
PopID = 5
PNuc = 6
QNuc = 7
N = 8
P = 9

contig_list = []

population_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

#open an output file for the whitelist:
whitefile = open('/PATH_TO_DATA/Whitelist_FILTERED_STRUCTURE.txt', 'wt')


for line in population_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<7:
		i=i+1
	else:
		current_line = line.rstrip('\n').split('\t')
		if current_line[Chr] not in contig_list:
			contig_list.append(current_line[Chr])
			whitefile.write(str(current_line[LocusID]) + '\t' + str(current_line[Col]) + '\n')

population_file.close()
whitefile.close()

