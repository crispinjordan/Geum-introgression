# This program will count the different types of polymorphism outcomes for snps in our pure population samples.

# This program ALSO looks for contigs where there is excess heterozygosity (i.e., Het_obs > 0.5, based on STACKS output in the batch_1.sumstats.tsv file), and produced a long series of 'grep' commands that can be used to remove these contigs from the dataset used for snp analyses.

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

#Define the counters for the number of polymorphism cases considered (they are listed here in the same order as being examined, below):
altfixed = 0
altfixedOpp = 0
rivPoly1 = 0
urbPoly1 = 0
rivPoly2 = 0
urbPoly2 = 0
sharedPoly = 0

totalCases = 0

candidateHomeolog = 0

#list to contain present line's info:
previous_line = []
present_line = []

#We will create a list to hold all of the contigs with troublesome Het_observed (i.e., Het_obs > 0.5 - in all cases we also find negative Fis - We also find many other cases of negative Fis, But I place less faith in these values.  For example, in some cases, we obviously have somehting wrong with the Het_obs, because the allele frequency at the site can be 0.5 but heterozygosity = 1.0, so that every individual is heterozygous - this smells of paralogy problems.)

chr_list = []



#This file contains the SNP data from STACKS to be analyzed:
population_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

#open an output file for count data:
cntfile = open('/PATH_TO_DATA/TrueVsBer.txt', 'wt')

#Print cases where we have either high heterozygosity, negaive Fis, or both - these are not necssarily useful, but what the hell?:
HetAndFisfile = open('PATH_TO_DATA/HetAndFis.txt', 'wt')
Hetfile = open('PATH_TO_DATA/Het.txt', 'wt')
Fisfile = open('/PATH_TO_DATA/Fis.txt', 'wt')

#This is also not useful, but I never deleted it from an earlier version - it includes a list of grep commands to remove troublesome loci.
#But, the current version of the program crates a whitelist of loci to be used by STACKS - this is done at the very bottom of the file.
grepfile = open('/PATH_TO_DATA/GrepCommands.txt', 'wt')


#Add a counter to count the number of contigs - start at 1 due to how I do the counting, below:
contigTotal = 1

for line in population_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<3:
		i=i+1
	elif i==3:
		previous_line = line.rstrip('\n').split('\t')
		i = i+1
		#Here, I simply list possibly troublesome contigs with respect to heterozygosity (I do not do anything further with them in this python program)
		if float(previous_line[HetObs]) > 0.5 and float(previous_line[Fis]) < 0:
			HetAndFisfile.write(str(previous_line)+'\n')
			candidateHomeolog = candidateHomeolog + 1
			chr_list.append(previous_line[Chr])
		if float(previous_line[HetObs]) > 0.5 and float(previous_line[Fis]) > 0:
			Hetfile.write(str(previous_line)+'\n')
			candidateHomeolog = candidateHomeolog + 1
			chr_list.append(previous_line[Chr])
#		if float(previous_line[HetObs]) < 0.5 and float(previous_line[Fis]) < 0:
		if float(previous_line[Fis]) < 0:
			Fisfile.write(str(previous_line)+'\n')
			candidateHomeolog = candidateHomeolog + 1
			chr_list.append(previous_line[Chr])
	else:
		present_line = line.rstrip('\n').split('\t')
		#print str(i)
		i = i +1
		#Here, I simply list possibly troublesome contigs with respect to heterozygosity (I do not do anything further with them in this python program)
		if float(present_line[HetObs]) > 0.5 and float(present_line[Fis]) < 0:
			HetAndFisfile.write(str(present_line)+'\n')
			if present_line[Chr] not in chr_list:
				chr_list.append(present_line[Chr])
				candidateHomeolog = candidateHomeolog + 1
		if float(present_line[HetObs]) > 0.5 and float(present_line[Fis]) > 0:
			Hetfile.write(str(present_line)+'\n')
			if present_line[Chr] not in chr_list:
				chr_list.append(present_line[Chr])
				candidateHomeolog = candidateHomeolog + 1
#		if float(present_line[HetObs]) < 0.5 and float(present_line[Fis]) < 0:
		if float(present_line[Fis]) < 0:
			Fisfile.write(str(present_line)+'\n')
			if present_line[Chr] not in chr_list:
				chr_list.append(present_line[Chr])
				candidateHomeolog = candidateHomeolog + 1

		#Count the total number of contigs in the batch_1.sumstats.tsv file:
		if present_line[Chr] != previous_line[Chr]:
			contigTotal = contigTotal + 1


#NOTE:  Depending on the content of STACKS' batch1.sumstats.tsv file, it may be necessary to alter the PopID's being compared, below.  In the current implementation, the function looks for SNPs that differ in varions ways between two classes of urbanum samples, urbanum from Berwickshire (urb_ber) vs from 'allopatric' sites (urb_true)


		#Check whether we are examining the same snp (LocusID, Chr, BP, Col) but that (1) species differ (PopID), (2) different alleles predominate in each species, noting that STACKS v1.31 places alternately fixed alleles in a diagonal in the batch_1.sumstats.tsv file
		# Check whether we have alternately fixed snps:
		if previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[PNuc] != present_line[QNuc] and previous_line[QNuc] =='-' and present_line[PNuc] =='-':
			altfixed = altfixed + 1
			#print str(previous_line)
			#print str(present_line)
		
		# Check whether STACKS ever lists alternately fixed alleles in the opposite order / orientation (this should NOT occur):
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[QNuc] != present_line[PNuc] and previous_line[PNuc] =='-' and present_line[QNuc] =='-':
			altfixedOpp = altfixedOpp + 1
			#print str(previous_line)
			#print str(present_line)

		# Count instances of polymorphism in rivale where urbanum has one of the two segregating alleles:
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[PNuc] == present_line[PNuc] and previous_line[QNuc] !='-' and present_line[QNuc] =='-':
			rivPoly1 = rivPoly1 + 1
			#print str(previous_line)
			#print str(present_line)

		# Count instances of polymorphism in urbanum where rivale has one of the two segregating alleles:
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[PNuc] == present_line[PNuc] and previous_line[QNuc] =='-' and present_line[QNuc] !='-':
			urbPoly1 = urbPoly1 + 1
			#print str(previous_line)
			#print str(present_line)

		# Count instances of polymorphism in rivale where urbanum has one of the two segregating alleles, but considering the case where STACKS presents the data in the opposite orientation as for rivPoly1:
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[QNuc] == present_line[QNuc] and previous_line[PNuc] !='-' and present_line[PNuc] =='-':
			rivPoly2 = rivPoly2 + 1
			#print str(previous_line)
			#print str(present_line)

		# Count instances of polymorphism in urbanum where rivale has one of the two segregating alleles but considering the case where STACKS presents the data in the opposite orientation as for urbPoly1:
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[QNuc] == present_line[QNuc] and previous_line[PNuc] =='-' and present_line[PNuc] !='-':
			urbPoly2 = urbPoly2 + 1
			#print str(previous_line)
			#print str(present_line)


		# Count instances of polymorphism present in both species:
		elif previous_line[BatchID]==present_line[BatchID] and previous_line[LocusID]==present_line[LocusID] and previous_line[Chr]==present_line[Chr] and previous_line[BP]==present_line[BP] and previous_line[Col]==present_line[Col] and previous_line[PopID] == 'urb_ber' and present_line[PopID] == 'urb_true' and previous_line[PNuc] == present_line[PNuc] and previous_line[QNuc] == present_line[QNuc] and previous_line[PNuc] !='-' and present_line[PNuc] !='-' and previous_line[QNuc] !='-' and present_line[QNuc] !='-':
			sharedPoly = sharedPoly + 1
			#print str(previous_line)
			#print str(present_line)


		#re-assign the current line as the previous line
		previous_line=present_line

#Print to screen the number of possible homeologs:
print('The number of possible homoelogs based on heterozygosity and Fis equals:  ' + str(candidateHomeolog))
#print(str(chr_list)+'\n')
print('The number of troublesome contigs is:  '+str(len(chr_list)))
#If there are more than 2 troublesome loci, then create a list of grep commands that will allow me to easily remove these loci from the original file:
if len(chr_list)>2:
	grepfile.write('grep -v ' + str(chr_list[0]) + ' batch_1.sumstats.tsv | ')
	for i in range(1,len(chr_list)-1):
		grepfile.write('grep -v ' + str(chr_list[i]) + ' | ')
	grepfile.write('grep -v ' + str(chr_list[len(chr_list)-1]) + ' > batch_1.sumstatsH_obsLessThanHalf.tsv')
#Add all of the counters together:
totalCases = altfixed + altfixedOpp + rivPoly1 + urbPoly1 + rivPoly2 + urbPoly2 + sharedPoly
#Print output to file:
cntfile.write('The number of alternately fixed alleles (orientation 1) equals:' + str(altfixed)+'\n')
cntfile.write('The number of alternately fixed alleles (orientation 2) equals:' + str(altfixedOpp)+'\n'+'\n')
cntfile.write('The number of cases where rivale is polymorphic and urbanum is fixed for 1 rivale allele (orientation 1) equals:' + str(rivPoly1)+'\n')
cntfile.write('The number of cases where rivale is polymorphic and urbanum is fixed for 1 rivale allele (orientation 2) equals:' + str(rivPoly2)+'\n'+'\n')
cntfile.write('The number of cases where urbanum is polymorphic and rivale is fixed for 1 urbanum allele (orientation 1) equals:' + str(urbPoly1)+'\n')
cntfile.write('The number of cases where urbanum is polymorphic and rivale is fixed for 1 urbanum allele (orientation 2) equals:' + str(urbPoly2)+'\n'+'\n')
cntfile.write('The number of cases where both species are polymorphic equals:' + str(sharedPoly)+'\n'+'\n')
cntfile.write('The number of cases considered:' + str(totalCases)+'\n')

population_file.close()
cntfile.close()
HetAndFisfile.close()
Hetfile.close()
Fisfile.close()
grepfile.close()

########################################
#Now, create a white-list of markers, using the chr_list, above:

#Open the batch_1.sumstats.tsv file again:

populationAgain_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

whitelist_loci = []

i=0

for line in populationAgain_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<3:
		i=i+1
	else:
		current_line = line.rstrip('\n').split('\t')
		#Here, I check whether the current line involves a flagged chromosome and whether we have white-listed the locus yet - if neither is true, we add this locus to the whitelist
		if (current_line[Chr] not in chr_list) and (current_line[LocusID] not in whitelist_loci):
			whitelist_loci.append(current_line[LocusID])

whitelist_file = open('/home/crispinjordan/AllWork/RAD_data/berwickshire_ddRAD/alignedSTACKS/stacksM6CorrAlt/W_whitelist_HetAndFis_UK_SamplesFILTERED_ANALYSES_NoM.txt', 'wt')

print('the number of loci on the whitelist equals:  ' + str(len(whitelist_loci)))
print('the number of contigs equals:  ' + str(contigTotal))


for white in whitelist_loci:
	whitelist_file.write(white + '\n')

populationAgain_file.close()
whitelist_file.close()
