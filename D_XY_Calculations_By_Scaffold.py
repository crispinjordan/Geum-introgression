# This program calculates Dxy between two samples - calculates dxy separately per scaffold.

#Specify the length of the loci:
locus_length = 117

#Make a counter and have it only start doing anything when we've passed the first 3 lines.

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
HetObs = 10
Fis = 17


#list to contain present line's info:
previous_line = []
present_line = []

d_xy_current = 0.0
d_xy_total = 0.0


#This file contains the data to be analyzed:
population_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

dxy_file = open("/PATH_TO_DATA/Dxy_pure_FullFilter.txt",'wt')

dxy_file.write("Scaffold" + "\t" + "no_RAD_loci" + "\t" + "d_xy_total" + "\t" + "Mean_Dxy" + "\n")

no_RAD_loci = 0
d_xy_total = 0

for line in population_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<3:
		i=i+1
	elif i==3:
		previous_line = line.rstrip('\n').split('\t')
		i = i+1
	else:
		present_line = line.rstrip('\n').split('\t')
		i = i +1
		#Assign locusID flags for current and previous line:
		prev_LocusID = previous_line[LocusID]
		pres_LocusID = present_line[LocusID]
		if prev_LocusID != pres_LocusID:
			no_RAD_loci = no_RAD_loci + 1
		#Assign chromosome (scaffold) flags for current and previous line:
		prev_Chr = previous_line[Chr]
		pres_Chr = present_line[Chr]
		#Check whether switch scaffolds:
		if prev_Chr != pres_Chr:
			dxy_file.write(str(prev_Chr) + "\t" + str(no_RAD_loci) + "\t" + str(d_xy_total) + "\t" + str(d_xy_total / (no_RAD_loci*locus_length) ) + "\n")
			no_RAD_loci = 0
			d_xy_total = 0
			
		#Check whether we are examining the same snp (LocusID, Chr, BP, Col) but that (1) species differ (PopID), (2) different alleles predominate in each species, noting that STACKS v1.31 places alternately fixed alleles in a diagonal in the batch_1.sumstats.tsv file
#**********	If analyzing 'true' populations:
		if previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'riv_true' and present_line[PopID]== 'urb_true':
#**********	If analyzing 'ber' populations:
#		if previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'riv_ber' and present_line[PopID]== 'urb_ber':
			d_xy_current = float(previous_line[P]) * (1.0-float(present_line[P])) + (1.0-float(previous_line[P]))*float(present_line[P])
			d_xy_total = d_xy_total + d_xy_current

		#re-assign the current line as the previous line
		previous_line=present_line

no_RAD_loci = no_RAD_loci + 1
dxy_file.write(str(prev_Chr) + "\t" + str(no_RAD_loci) + "\t" + str(d_xy_total) + "\t" + str(d_xy_total / (no_RAD_loci*locus_length) ) + "\n")

population_file.close()
dxy_file.close()
