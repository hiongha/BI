import re

must_list = []
with open('../mustgiven_CRC.txt','r') as F1:
	F1.readline()
	for line in F1:
		gene,aa,exon,chr,start,end,cds_mut = line.strip().split('\t')
		id = '_'.join([gene,chr,start,end,cds_mut])
		must_list.append(id)

brief = {}
with open('crc_brief_from_lsx','r') as F2:
	F2.readline()
	for line in F2:
		info = line.strip().split('\t')
		gene = info[1]	
		chr = info[4]
		start = info[5]
		end = info[6]
		cds_mut = info[7]
		type = info[20]
		id = '_'.join([gene,chr,start,end,cds_mut])
		brief[id] = type

for i in must_list:
#	print i
	if i in brief.keys():
		print brief[i]
	else:
		print "meiyou"
		
	

