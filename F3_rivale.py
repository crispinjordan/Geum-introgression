# This program calculates f3 for the Berwickshire rivale samples.

#There are many cases where we have polymorphism within the Berwickshire samples but the allopatric samples are fixed for the same allele.  I'm not sure what this tells us about introgression.  So, I also output data to a 'clean' file, where polymorphism must occur either within or between the allopatric populations.

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


#This file contains the data to be analyzed:
population_file = open("/PATH_TO_DATA/batch_1.sumstats.tsv")

f3_file = open("/PATH_TO_DATA/f3_rivale.txt",'wt')
f3_file_clean = open("/PATH_TO_DATA/f3_rivale_clean.txt",'wt')

f3_file.write("scaffold" + "\t" + "locus" + "\t" + "Col" + "\t" + "riv_allo_freq" + "\t" + "urb_allo_freq" + "\t" + "riv_Ber_freq" + "\n")
f3_file_clean.write("scaffold" + "\t" + "locus" + "\t" + "Col" + "\t" + "riv_allo_freq" + "\t" + "urb_allo_freq" + "\t" + "riv_Ber_freq" + "\n")

for line in population_file:
#Skip the first 3 lines (lines 0,1,2) of the file:
	if i<4:
		i=i+1
	elif i==4:
		previous_line = line.rstrip('\n').split('\t')
		i = i+1
		if previous_line[PopID] == 'riv':
			riv_freq = previous_line[P]
		elif previous_line[PopID] == 'urb':
			urb_freq = previous_line[P]
		elif previous_line[PopID] == 'riv_ber':
			riv_ber_freq = previous_line[P]
		else:
			print("oops")
	else:
		present_line = line.rstrip('\n').split('\t')
		i = i +1
		#Assign locusID flags for current and previous line:
		prev_LocusID = previous_line[LocusID]
		pres_LocusID = present_line[LocusID]
		prev_Col = previous_line[Col]
		pres_Col = present_line[Col]

		if prev_LocusID == pres_LocusID and prev_Col == pres_Col:		
			if present_line[PopID] == 'riv':
				riv_freq = present_line[P]
			elif present_line[PopID] == 'urb':
				urb_freq = present_line[P]
			elif present_line[PopID] == 'riv_ber':
				riv_ber_freq = present_line[P]
		if prev_LocusID != pres_LocusID or prev_Col != pres_Col:
			f3_file.write(str(str(previous_line[Chr]) + "\t" + previous_line[LocusID]) + "\t" + str(previous_line[Col]) + "\t" + str(riv_freq) + "\t" + str(urb_freq) + "\t" + str(riv_ber_freq) + "\n")
			if (float(riv_freq) == 1.0 and float(urb_freq) != 1) or (float(riv_freq) !=1 and float(urb_freq) == 1) or (float(riv_freq) !=1 and float(urb_freq) != 1):
				f3_file_clean.write(str(previous_line[Chr]) + "\t" + str(previous_line[LocusID]) + "\t" + str(previous_line[Col]) + "\t" + str(riv_freq) + "\t" + str(urb_freq) + "\t" + str(riv_ber_freq) + "\n")
			if present_line[PopID] == 'riv':
				riv_freq = present_line[P]
			elif present_line[PopID] == 'urb':
				urb_freq = present_line[P]
			elif present_line[PopID] == 'riv_ber':
				riv_ber_freq = present_line[P]
		#re-assign the current line as the previous line
		previous_line=present_line

f3_file.write(str(previous_line[Chr]) + "\t" + str(previous_line[LocusID]) + "\t" + str(previous_line[Col]) + "\t" + str(riv_freq) + "\t" + str(urb_freq) + "\t" + str(riv_ber_freq) + "\n")
if (float(riv_freq) == 1.0 and float(urb_freq) != 1) or (float(riv_freq) !=1 and float(urb_freq) == 1) or (float(riv_freq) !=1 and float(urb_freq) != 1):
	f3_file_clean.write(str(previous_line[Chr]) + "\t" + str(previous_line[LocusID]) + "\t" + str(previous_line[Col]) + "\t" + str(riv_freq) + "\t" + str(urb_freq) + "\t" + str(riv_ber_freq) + "\n")

population_file.close()
f3_file.close()
f3_file_clean.close()
