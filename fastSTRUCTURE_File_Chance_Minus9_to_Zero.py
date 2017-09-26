#Changes 0 to -9 in STRUCTURE file 

o = open('/PATH_TO_DATA/batch_1.structure2.edited.str', 'w')
for l in open('/PATH_TO_DATA/batch_1.structure2.str'):
	info = l.strip().split('\t')
	for i in range(len(info[2:])):
		if info[i+2] == '0': info[i+2] = '-9'
	o.write('\t'.join(info) + '\n')
