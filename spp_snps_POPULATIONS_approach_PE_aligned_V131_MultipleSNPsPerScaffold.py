#Make a counter and have it only start doing anything when we've passed the first 3 lines.

#Line counter for reading file:
i=0
#variable to remember what the previously printed LocusID was, so that we only print data from unique LocusIDs - initialize as 0, as I think this will not occur in STACKS:
uniqLocusID=0

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

#list to contain present line's info:
previous_line = []
present_line = []

#population_file = open("/home/crispinjordan/AllWork/RAD_data/berwickshire_ddRAD/alignedSTACKS/stacksM6CorrAlt/batch_1.sumstatsH_obsLessThanHalf.tsv")
population_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

#open an output file for the whitelist:
whitefile = open('/PATH_TO_DATA/Whitelist_M6_PE_CorrAltV135_Filtered_NoM_MultipleSNPsPerScaffold_April28.txt', 'wt')
#Open an output file that will contain all the information about a locus, to be used in a 'dict' in a following analysis:
dictfile = open('/PATH_TO_DATA/DictInfo_M6_PE_CorrAltV135_Filtered_NoM_MultipleSNPsPerScaffold_April28.txt', 'wt')

for line in population_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<3:
		i=i+1
	elif i==3:
		previous_line = line.rstrip('\n').split('\t')
		i = i+1
	else:
		present_line = line.rstrip('\n').split('\t')
		#print str(i)
		i = i +1
		#Check whether we are examining the same snp (LocusID, Chr, BP, Col) but that (1) species differ (PopID), (2) different alleles predominate in each species, noting that STACKS v1.31 places alternately fixed alleles in a diagonal in the batch_1.sumstats.tsv file:
		if previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] != present_line[PopID] and previous_line[PNuc] != present_line[QNuc] and previous_line[QNuc] =='-' and present_line[PNuc] =='-':
			print str(previous_line)
			print str(present_line)
			#print the output to the whitelist:
			whitefile.write(str(present_line[LocusID])+'\t'+str(present_line[Col])+'\n')
			#print information for the dict:
			dictfile.write(str(previous_line[BatchID])+'\t'+str(present_line[BatchID])+'\t'+str(previous_line[LocusID])+'\t'+str(present_line[LocusID])+'\t'+str(previous_line[Chr])+'\t'+str(present_line[Chr])+'\t'+str(previous_line[BP])+'\t'+str(present_line[BP])+'\t'+ str(previous_line[Col])+'\t'+str(present_line[Col])+'\t'+str(previous_line[PopID])+'\t'+str(present_line[PopID])+'\t'+str(previous_line[PNuc])+'\t'+str(present_line[PNuc])+'\t'+str(previous_line[QNuc])+'\t'+str(present_line[QNuc])+'\n')
		#re-assign the current line as the previous line
		previous_line=present_line

population_file.close()
whitefile.close()
dictfile.close()
