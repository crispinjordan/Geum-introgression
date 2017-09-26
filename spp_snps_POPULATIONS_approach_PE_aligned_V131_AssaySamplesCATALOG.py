# This program takes information from the file, DictInfo_M6_PE_CorrAltV131.txt, in order to assay the population samples and determine the alleles present within them.
# It which uses:
#	- the dict information from DictInfo_M6_PE_CorrAltV131.txt
#	- the information from xxx.matches.tsv in order to link the locus ID from the populations output (in DictInfo_M6_PE_CorrAltV131.txt) to the locus ID in the xxxx.snps.tsv file (i.e., the Locus ID in the Populations output refers, confusingly, to the catalog locus, whereas the locus ID in the snps.tsv file is the locus number IN A FOCAL SAMPLE).  We need to match the two different meanings of Locus ID
#	- xxxx.snps.tsv file for each individual, which will have the genotype information for each locus.



#NOTE - this file must be run using the correct:
#					- DictInfo file
#					- xxxx.matches.tsv file


import gzip
from Bio import SeqIO


#Open the output file to contain all the genotyping data:
output_snps_file =  open('/PATH_TO_DATA/GenotypeProportionsPerIndividual_Filtered_NoM.txt', 'wt')
output_snps_file.write('Sample' + '\t' + 'sample_type' + '\t' + 'loci_genotyped' + '\t' + 'rivale_allele_count' + '\t' + 'urbanum_allele_count' + '\t' + 'Third_allele_at_a_locus' + '\t' + 'No_Heterozygous_loci' + '\n')

#Open an additional file, where each line has the genotype for a sample at a focal locus:
output_indiv_genoptypes = open('/PATH_TO_DATA/GenotypePerLocusPerIndividual_Filtered_NoM.txt', 'wt')
output_indiv_genoptypes.write('sample' + '\t' + 'sample_type' + '\t' + 'contig' + '\t' + 'catalogLocus' + '\t' + 'column' + '\t' + 'rivale_allele' + '\t' + 'urbanum_allele' + '\t' + '\n')


#Number the columns within the imported data strings for each spp-specific snp, from the file DictInfo_M6_PE_CorrAltV131.txt:
previous_line_BatchID = 0
present_line_BatchID = 1
previous_line_LocusID = 2
present_line_LocusID = 3
previous_line_Chr = 4
present_line_Chr = 5
previous_line_BP = 6
present_line_BP = 7
previous_line_Col = 8
present_line_Col = 9
previous_line_PopID_Rivale = 10
present_line_PopID_urbanum = 11
previous_line_PNuc_RivaleAllele = 12
present_line_PNuc = 13
previous_line_QNuc = 14
present_line_QNuc_UrbanumAllele = 15



#list to accept the current info from DictInfo:
current_line = []
#list to accept the current info from sumstats.tsv:
#sumstats_line = []

#Create a dict to hold the information about each species-specific snp:
dict_spp_snps = {}

############################################################
#Notation for gaining information from the matches.tsv file:
############################################################
#Create a dict for the xxxx.matches.tsv file that matches the Catalog (Locus) ID to the Locus ID for the focal sample:
dict_matches_ID_info = {}

#List of order of entries for each line in the xxxx.matches.tsv file:
SQL_ID_matches = 0
Batch_ID_matches = 1
Catalog_ID_matches = 2
Sample_ID_matches = 3
Locus_ID_matches = 4
Haplotype_matches = 5
Stack_Depth_matches = 6
Log_likelihood_matches = 7

#List of columns in sample_snps_file that come from xxx.snps.tsv:
Sql_ID_snps = 0
Sample_ID_snps = 1
Locus_ID_snps = 2
Column_snps = 3
Type_snps = 4
Likelihood_ratio_snps = 5
Rank_1_snps = 6
Rank_2_snps = 7


#Use a Population Map from STACKS to both list the samples to be analyzed and to create a dict that will enumerate each sample's 'type', to be included in the output - this also includes information on the  hypothesized hybrid class:
sample_list_file = open("/PATH_TO_DATA/popmap_BerwickshireSampleListByIndiviudalTypePlusTruePure_4CatPlusTrue.txt")



#Use a dict to store what type of population each sample derives from:
sample_type_dict = {}
#Use a list to iterate through all the samples:
sample_list = []

for smp in sample_list_file:
	current_sample = smp.rstrip('\n').split('\t')
	#Add sample to list to be processed:
	sample_list.append(current_sample[0])
	#Add sample and it's population type to the dict:
	sample_type_dict[current_sample[0]] = str(current_sample[1])


#This file contains the data to be analyzed (i.e, the species-specific snps to be genotyped):
snp_info_file = open("/PATH_TO_DATA/DictInfo_M6_PE_CorrAltV135_Filtered_NoM.txt")

for line in snp_info_file:
	current_line = line.rstrip('\n').split('\t')
	#
	dictkey = str(current_line[previous_line_LocusID])
	dict_spp_snps[dictkey] = current_line


##########################################################################################################################################
#Now, import the matches.tsv file to be used to link the two types of locus ID
#In order to process each sample, we need to look at the matches.tsv file and snp file for each sample; so, create a list to work through
##########################################################################################################################################


#Iterate the following over each sample: **********************************************************************************************

for sample in sample_list:
	#Now, create the names of the files to be opened:


	sample_matches_file = gzip.open('/PATH_TO_DATA/' + sample + '.matches.tsv.gz')
	sample_snps_file = gzip.open('/PATH_TO_DATA/' + sample + '.snps.tsv.gz')
	
	#Make a Dict that can link the Catalog (Locus) ID to the Locus ID in the focal sample:
	#Start by emptying the Dict:
	dict_matches_ID_info.clear()
	#Set a flag that can be used to tell whether a locus should be added to a dict
	LocusPresent = 0
	
	#Counter to skip line:
	i = 0

	#This 'for' statement creates the dict that matches the two forms of Locus ID for the Catalog vs snps files:
	for matches_line in sample_matches_file:
		#Skip the first line of the file:
		if i==0:
			i=i+1
		else:
			# read a line of the matches file:
			current_line_matches = matches_line.rstrip('\n').split('\t')

			#print(current_line_matches)

			# check whether the catalog locus for that line is in the dict for the sepecies-specific snps:
			temp_matches_CHECK = str(current_line_matches[Catalog_ID_matches])

			#The Catalog_ID can appear more than once, so we check whether this locus is already stored in dict_matches_ID_info
			if temp_matches_CHECK in dict_spp_snps and str(current_line_matches[Locus_ID_matches]) not in dict_matches_ID_info:
				LocusPresent = 1
			#NOTE: we include the LocusPresent == 1 flag because the locus may be listed more than once, and this approach avoids adding it multiple times
			if LocusPresent == 1:
				# Assign the **Locus ID** as the key and the Catalog ID as the value (this allows one match the Locus ID in the snps.tsv file directly to the Catalog ID via this Dict):
				dict_matches_ID_info[current_line_matches[Locus_ID_matches]] = current_line_matches[Catalog_ID_matches]
				#Need to re-set the LocusPresent marker so that we only add a locus to dict_matches_ID_info if it is correct to do so (leaving it as ==1 will add all loci)
				LocusPresent = 0






	#Before processing the focal sample, set the counters of the data to zero:
	#Counter of the number of loci genotyped:
	loci_genotyped = 0
	#Counter of the number of rivale alleles found in a sample:
	rivale_allele_count = 0
	#Counter of the number of urbanum alleles found in a sample:
	urbanum_allele_count = 0
	#number of loci that have third alleles are found in a sample:
	Third_allele_at_a_locus = 0
	#A counter to count errors:
	ErrorCount = 0
	#A counter of the number of heterozygous loci:
	het_loci = 0


	#Now that we have a way of translating 'Locus ID' between the file types, we can now genoptype the focal sample:	
	#For each line, check whether the Locus ID for this sample is present in dict_matches_ID_info
	for calling_snp_line in sample_snps_file:
		#Read line of data for a given nucleitide position:
		focal_snp = calling_snp_line.rstrip('\n').split('\t')
		#Check whether the Locus for this nucleotide is present in dict, dict_matches_ID_info (which matched the Locus ID of the focal sample (Key) to the Locus ID in the catalog (Value)):
		if focal_snp[Locus_ID_snps] in dict_matches_ID_info:
			#When the focal Locus ID is one of the spp-specific loci, check whether this is the correct nucleotide postion:
			#Note that 'dict_matches_ID_info[focal_snp[Locus_ID_snps]]' calls for the Catalog ID for this locus, so calling dict_spp_snps[dict_matches_ID_info[focal_snp[Locus_ID_snps]]] will obtain the line of data is associated with this focal locus in STACKS' catalog.
			temp_line_snps = dict_spp_snps[dict_matches_ID_info[focal_snp[Locus_ID_snps]]]  #temp_line_snps holds the infor for this focal Catalog Locus:
				#Determine whether (a) we have the correct Column of data, and (b) whether the genptype of the nucleotide is known:
			if temp_line_snps[previous_line_Col] == focal_snp[Column_snps] and str(focal_snp[Type_snps]) != 'U':
				

				if str(focal_snp[Type_snps]) == 'O':
					#Compare the majority allele to the dictionary of snps:
					if focal_snp[Rank_1_snps] == temp_line_snps[previous_line_PNuc_RivaleAllele]:
						rivale_allele_count = rivale_allele_count + 2
						loci_genotyped = loci_genotyped + 1
						output_indiv_genoptypes.write(sample + '\t' + str(sample_type_dict[sample]) + '\t' + temp_line_snps[previous_line_Chr] + '\t' + temp_line_snps[previous_line_LocusID] + '\t' + temp_line_snps[previous_line_Col] + '\t' + str(temp_line_snps[previous_line_PNuc_RivaleAllele]) + '\t' + '-' + '\t' + '\n')
					elif focal_snp[Rank_1_snps] == temp_line_snps[present_line_QNuc_UrbanumAllele]:
						urbanum_allele_count = urbanum_allele_count + 2
						loci_genotyped = loci_genotyped + 1
						output_indiv_genoptypes.write(sample + '\t' + str(sample_type_dict[sample]) + '\t' + temp_line_snps[previous_line_Chr] + '\t' + temp_line_snps[previous_line_LocusID] + '\t' + temp_line_snps[previous_line_Col] + '\t' + '-' + '\t' + str(temp_line_snps[present_line_QNuc_UrbanumAllele]) + '\t' + '\n')
					else:
						Third_allele_at_a_locus = Third_allele_at_a_locus + 1
						#loci_genotyped = loci_genotyped + 1
						print('3rd allele present for this focal snp of sample:' + sample + '\n')
						print(str(focal_snp))
						print('...and for this focal spp-specific snp info:')
						print(str(temp_line_snps))
				elif str(focal_snp[Type_snps]) == 'E':
					#Check the 2 possible orientations for having a sample be heterozygous for the correct species-specific alleles: 
					if focal_snp[Rank_1_snps] == temp_line_snps[previous_line_PNuc_RivaleAllele] and focal_snp[Rank_2_snps] == temp_line_snps[present_line_QNuc_UrbanumAllele]:
						#Increment each allele type count by 1, since heterozygous:
						rivale_allele_count = rivale_allele_count + 1
						urbanum_allele_count = urbanum_allele_count + 1
						#Increment the locus counter:
						loci_genotyped = loci_genotyped + 1
						het_loci = het_loci + 1
						output_indiv_genoptypes.write(sample + '\t' + str(sample_type_dict[sample]) + '\t' + temp_line_snps[previous_line_Chr] + '\t' + temp_line_snps[previous_line_LocusID] + '\t' + temp_line_snps[previous_line_Col] + '\t' + str(temp_line_snps[previous_line_PNuc_RivaleAllele]) + '\t' + str(temp_line_snps[present_line_QNuc_UrbanumAllele]) + '\t' + '\n')
					elif focal_snp[Rank_1_snps] == temp_line_snps[present_line_QNuc_UrbanumAllele] and focal_snp[Rank_2_snps] == temp_line_snps[previous_line_PNuc_RivaleAllele]:
						#Increment each allele type count by 1, since heterozygous:
						rivale_allele_count = rivale_allele_count + 1
						urbanum_allele_count = urbanum_allele_count + 1
						#Increment the locus counter:
						loci_genotyped = loci_genotyped + 1
						het_loci = het_loci + 1
						output_indiv_genoptypes.write(sample + '\t' + str(sample_type_dict[sample]) + '\t' + temp_line_snps[previous_line_Chr] + '\t' + temp_line_snps[previous_line_LocusID] + '\t' + temp_line_snps[previous_line_Col] + '\t' + str(temp_line_snps[previous_line_PNuc_RivaleAllele]) + '\t' + str(temp_line_snps[present_line_QNuc_UrbanumAllele]) + '\t' + '\n')
					#If neither situation applies, presumably there is a 3rd allele present:
					else:
						Third_allele_at_a_locus = Third_allele_at_a_locus + 1
						#loci_genotyped = loci_genotyped + 1
						print('3rd allele present for this focal snp of sample:' + sample + '\n')
						print(str(focal_snp))
						print('...and for this focal spp-specific snp info:')
						print(str(temp_line_snps))
				else:
					print('Big error in ' + sample + ' because unknown Type of genotype' + '\n')
					ErrorCount = ErrorCount + 1

	#For DE-BUGGING only:
	print('Number of big errors: ' + str(ErrorCount) + '\n')



	#Write the results from genotyping this sample to the output file:
	output_snps_file.write(sample + '\t' + str(sample_type_dict[sample]) + '\t' + str(loci_genotyped) + '\t' + str(rivale_allele_count) + '\t' + str(urbanum_allele_count) + '\t' + str(Third_allele_at_a_locus) + '\t' + str(het_loci) + '\n')

	#Close the sample files (to be opened, again, for the next sample)
	sample_matches_file.close()
	sample_snps_file.close()



#De-bugging purposes:
#	for keys, values in dict_matches_ID_info.items():
#		print('key:' + keys)
#		print(values)


output_snps_file.close()
snp_info_file.close()
sample_list_file.close()
output_indiv_genoptypes.close()
