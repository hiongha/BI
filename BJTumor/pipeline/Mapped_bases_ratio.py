#!/usr/bin/python
#xinghe <xingh3223@berryoncology.com>

from __future__ import division
import os
import re
import subprocess
import sys

ReadId_MatchLen = {}

FC = sys.argv[1] #181214_MN00311_0020_A000H2HWGH_0
sample = sys.argv[2] #cfDNA2018a-1
bam = FC+'/'+sample + '.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam'
sam = os.path.basename(bam).split('.bam')[0]+'.sam'
samF = open(sam,'w')
subprocess.call(['samtools view %s'%bam],shell = True, stdout = samF)
samF.close()

f=open(sam,'r')
for line in f:

	reads_id = line.strip().split()[0]
	ReadLen = len(line.strip().split('\t')[9])
	cigar = line.strip().split()[5]
	XF = re.search(r'XF:i:(\d+)\t',line.strip()).group(1)

	if re.findall(r'(\d+?M)',cigar) !=[]:
		match_num_list = re.findall(r'(\d+?M)',cigar)
		match_num_int_list = [int(each.strip('M')) for each in match_num_list]
		match_num_each_read = sum(match_num_int_list)
		if not ReadId_MatchLen.has_key(reads_id):
			ReadId_MatchLen[reads_id] = [match_num_each_read,ReadLen,int(XF)]
		elif match_num_each_read > ReadId_MatchLen[reads_id]:
			ReadId_MatchLen[reads_id] = [match_num_each_read,ReadLen,int(XF)]
		else:
			pass

total_base = 0
total_matched_base = 0

for each in ReadId_MatchLen:
	total_base += ReadId_MatchLen[each][1]*ReadId_MatchLen[each][2]
	total_matched_base += ReadId_MatchLen[each][0]*ReadId_MatchLen[each][2]

Mapped_bases_ratio = total_matched_base/total_base
print "total_matched_bp_to_human:"+str(total_base)
print "matched_bases:"+str(total_matched_base)
resultF = open("result.txt",'a')
resultF.write("Mapped_bases_ratio:" + str(Mapped_bases_ratio)+'\n')
resultF.close()

