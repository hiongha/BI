#!/usr/bin/python
#coding=utf-8
import sys
import re
from copy import deepcopy
file = sys.argv[1]

F = open(file,'r')
for line in F:
	info = line.strip().split('\t')
	ID,chr,start,cigar,seq,XF = info[0:6] #cigar 150M90S20M
	end = 0
	check = 1
	XF = re.sub(r'XF:i:(\d+)',r'\1',XF)
	cigar_tmp = re.sub(r'^(\d+)S',r'',cigar)
	cigar_tmp = re.sub(r'(\d+)S',r'',cigar)
	chai1 = re.split(r'[A-Za-z]',cigar_tmp) #150 90  20
	chai2 = re.split(r'(\d+)',cigar_tmp) #M S M

	for m in range(0,len(chai2)+1):
		if chai2[m] == 'M':
			if check == 1:
				end = start + chai1[m] -1
				check+=1
			else:
				start = end+1
				end = start + chai1[m]-1
		for p in range(start,end+1):
			key = '\t'.join([ID,chr,p,XF])
			if hash.has_key(key):
				pass
			else:
				hash[key] = 1
		elif chai2[m] == 'S' or chai2[m] == 'D':
			end = end + chai1[m]	

for each in hash.keys():
	print(each)







	












		
















