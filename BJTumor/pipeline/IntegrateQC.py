#!/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>

'''
将BZ666结果整合到csmart自动化分析得到的excel表格中.
'''
import argparse
import xlrd
import xlwt
from xlutils.copy import copy
import collections
import re
import os
from copy import deepcopy as deepcopy_from_copy
from Bio import SeqIO
import subprocess

parser=argparse.ArgumentParser(prog="",description="将BZ666QC结果整合到csmart结果中")
parser.add_argument("-b","--bz-qc",help='BZ666流程QC结果',default='result.txt')
parser.add_argument("-c","--cs-qc",help="csmart流程QC结果",default='181207_TPNB500312_0040_AH75JJBGX7_20.xls')
parser.add_argument("--brief-file",help="*.brief.xls文件",default='Z18L05895.brief.xls')
parser.add_argument("--type",help="样本类型",default='LC/CRC/...')
args=parser.parse_args()

#--------------------------初始化变量---------------------------------------

sample_type_brief = args.type
pipeline = '/share/work1/xinghe/proc/BZQC/pipeline'
root_dir = os.path.abspath('./')
must_given_LC_add = '/share/work1/xinghe/proc/BZQC/pipeline/config/mustgiven_%s_add.txt'%(sample_type_brief)
BZ_QC = root_dir +'/' + args.bz_qc
CS_QC = root_dir +'/' + args.cs_qc
CS_brief = root_dir + '/' + args.brief_file
OUT = re.match(r'(.+?_.+?_.+?_.+?)_.*.xls',CS_QC).group(1)+'_cSMART_%s_bz.xls'%(sample_type_brief)

#-------------------------对于在mustgiven和brief文件中都找不到的位点需要annovar注释--------
def TrRevSeq(seq):
	dic = {'A':'T','T':'A','G':'C','C':'G'}
	seq = seq[::-1]
	tr_seq = []
	for i in seq:
		tr_seq.append(dic[i.upper()])
	return ''.join(tr_seq)

def gene2nm(file = ''):
	dic = {}
	dic2 = {}
	with open(file,'r') as F:
		for line in F:
			gene,NM,strand = line.strip().split('\t')
			dic[gene] = NM
			dic2[gene] = strand
	return dic,dic2

gene2NM,gene2strand = gene2nm(file = '/share/Oncology/production/cSMART/cSMART/analysis/xinghe/DCE/test2/all_gene_strand.txt')

fa = '/share/public/database/genome/Homo_sapiens/UCSC-hg19/download/hg19.chr.fa'
def findseq(fa,chr,start,end = '' ):
	with open(fa,'r') as fa:
		for record in SeqIO.parse(fa,'fasta'):
			if record.id == chr:
				seq = record.seq
				onebpbeforestart = seq[int(start)-1-1]
				if end != '':
					xulie = ''.join(seq[int(start)-1:int(end)-1+1])
	if end == '':
		return onebpbeforestart
	elif end != '':
		return xulie,onebpbeforestart

def getrefalt(cds_mut,chr,start,end,gene):
	strand = gene2strand[gene]
	start = int(start)
	end = int(end)
	if re.search(r'>',cds_mut):
		onebpbeforestart = findseq(fa = fa,chr = chr,start = start)
		ref,alt = re.search(r'c.\d+_?\d+(?:\+?|\-?)(?:\d+)?([ATGC]*)>([ATGC]*)',cds_mut).groups()
		if strand == '-':
			ref = TrRevSeq(ref)
			alt = TrRevSeq(alt)
	elif re.search(r'delins[ATGC]+',cds_mut):
		ref,onebpbeforestart = findseq(fa = fa,chr = chr,start = start, end = end)
		alt = re.search(r'delins([ATGC]+)',cds_mut).group(1)
	elif re.search(r'delins\d+',cds_mut):
		ref,onebpbeforestart = findseq(fa = fa,chr = chr,start = start, end = end)
		alt = 'wujie'
	elif re.search(r'del([ATGC]+|\d+|$)',cds_mut):
#		start = start +1 
		ref,onebpbeforestart = findseq(fa = fa,chr = chr,start = start, end = end)
		alt = '-'
	elif re.search(r'\d+ins[ATGC]+',cds_mut):
		start = start + 1
		onebpbeforestart = findseq(fa = fa,chr = chr,start = start)
		ref = '-'
		alt = re.search(r'\d+ins([ATGC]+)',cds_mut).group(1)
		if strand == '-':
			ref = TrRevSeq(ref)
			alt = TrRevSeq(alt)
	elif re.search(r'\d+ins\d+',cds_mut):
		start = start +1
		onebpbeforestart = findseq(fa = fa,chr = chr,start = start)
		ref = '-'
		alt = 'wujie'
		if strand == '-':
			ref = TrRevSeq(ref)
			alt = TrRevSeq(alt)
	elif re.search(r'dup[ATGC]+',cds_mut):
		ref_all_must_given = re.search(r'dup([ATGC]+)',cds_mut).group(1)
		end = start + (len(ref_all_must_given))
		start = int(start) + 1
		ref,onebpbeforestart = findseq(fa = fa,chr = chr,start = start,end = end)
		alt = ref*2 
	elif re.search(r'dup\d+',cds_mut):
		dup_len = re.search(r'dup(\d+)',cds_mut).group(1)
		end = start + int(dup_len)
		start = start+1;
		ref,onebpbeforestart = findseq(fa = fa,chr = chr,start = start,end = end)
		alt = ref*2
	else:
		onebpbeforestart = 'xxxx'
		ref = 'xxxx'
		alt = 'xxxx'
	ref_tmp = '' if ref == '-' else ref
	alt_tmp = '' if alt == '-' else alt
	retAnnovar = '\t'.join([chr,str(int(start)-1),'.',onebpbeforestart+ref_tmp,onebpbeforestart+alt_tmp,'.','.','.'])+'\n'
	return ref,alt,retAnnovar

#def annovar(unfoundMutFunc,analysisPath):
def annovar(cds_mut,chr,start,end,gene,analysisPath):
	vcfList = []
	pre1 = analysisPath+'/list1';pre2 = analysisPath+'/list2'
	F4 = open('%s.annovar.hg19_multianno.txt'%(pre2),'a')
	F4.write('\t'.join(["Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene","Func.wgEncodeGencodeBasicV19","Gene.wgEncodeGencodeBasicV19","GeneDetail.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","Otherinfo"])+'\n')

	ref,alt,retAnnovar = getrefalt(cds_mut = cds_mut,chr = chr,start = start,end = end,gene = gene) 
	vcfList.append(retAnnovar)
	print(ref,alt,retAnnovar)
	with open(pre1+'.vcf','w') as F1:
		F1.write(retAnnovar)
	with open(pre2+'.vcf','a') as F2:
		F2.write(retAnnovar)
	p1 = subprocess.Popen('perl /share/public/database/Annovar/2016Feb01/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile %s.avinput %s.vcf'%(pre1,pre1),shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	p1.wait()
	p2 = subprocess.Popen('perl /share/public/database/Annovar/2016Feb01/table_annovar.pl --buildver hg19 --thread 4  --remove --otherinfo --protocol refGene,wgEncodeGencodeBasicV19 -operation g,g -nastring . %s.avinput /share/public/database/Annovar/2016Feb01/humandb --outfile %s.annovar'%(pre1,pre1),shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	p2.wait()
	p3 = subprocess.Popen('perl /share/public/database/Annovar/2016Feb01/annotate_variation.pl --geneanno --buildver hg19 --hgvs --thread 4 %s.avinput /share/public/database/Annovar/2016Feb01/humandb --outfile %s.hgvs'%(pre1,pre1),shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	p3.wait()
#		ann = subprocess.Popen('sh annotation.sh',shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	with open('%s.annovar.hg19_multianno.txt'%(pre1),'r') as F3:
		F3.readline()
		F4.write(F3.readline())
	with open('%s.annovar.hg19_multianno.txt'%(pre1),'r') as outF:
		outF.readline()
		mutFuncAnnovar = outF.readline().strip().split('\t',9)[8]
	os.system('mkdir -p %s/annovarTmp'%(analysisPath))
		#os.system('mv %s* %s* %s/annovarTmp'%(pre1,pre2,analysisPath))
	os.system('mv %s* %s/annovarTmp'%(pre1,analysisPath))
	print(mutFuncAnnovar)
	return mutFuncAnnovar

#-------------将BZ666QC结果添加到csmart QC结果中-----------------------------
#从result.txt中提取需要添加到excel表格中的信息.
info = collections.OrderedDict()
Header = collections.deque()
Value = collections.deque()
with open(BZ_QC,'r') as f:
	for line in f:
		head,value = line.strip().split(":")
		if head in ['Q20','Q30','Mapped_bases_ratio','Mapped_bases_ratio_to_target','Reads_uniquely_Mapped_bases_ratio_to_target','Bases_ratio_of_target_covered_at_least_5x','Bases_ratio_of_target_covered_at_least_10x']:
			value = float(value)*100
			value = '{0:.2f}'.format(value)
		Header.append(head)
		Value.append(value)
		info[head] = value
info_num = len(info)
#print( Value)
#将数据写入excel表格'QC'sheet中.
QC_rb = xlrd.open_workbook(CS_QC)
QC_rs = QC_rb.sheet_by_name('QC')
nrows = QC_rs.nrows
ncols = QC_rs.ncols
QC_wb = copy(QC_rb)
QC_ws = QC_wb.get_sheet("QC")
for col in range(ncols,int(ncols)+int(info_num)):
	QC_ws.write(0,col,Header.popleft())
	QC_ws.write(1,col,Value.popleft())
#QC_wb.save(out)				

#----------------在all.must.given中添加NM号和mutation_type列-----------------
#从mustgiven_LC_add.txt文件中获取染色体/起始位点/终止位点/cds_mut/NM编号/mutation_type信息.
NMF=open(must_given_LC_add,'r')
title=NMF.readline().strip().split('\t')
given_dict={}
for line in NMF:
	[chr,start,end,cds_mut,NM,mutation_type]=line.strip().split('\t')[3:]
	cds_mut = re.sub(r'(c\..*(?:ins|del))(?:\d+|\w+)',r'\1',cds_mut)
	key='\t'.join([chr,start,end,cds_mut])
	given_dict[key]=[NM,mutation_type]
	
#对于mustgiven_LC_add.txt文件中没有mutation_type的位点,从Z18L05895.brief.xls中寻找信息.
briefF=open(CS_brief,'r')
title=briefF.readline()
for line in briefF:
	chr,start,end = line.strip().split('\t')[0:3]
	mutation_type,cds_mut_info = line.strip().split('\t')[7:9]
	if cds_mut_info != '-':
		NM=re.search(r'(NM_\d+)',cds_mut_info).group(1)
		print(cds_mut_info)
		if re.search(r'(c.\S+):p.',cds_mut_info) != None:
			cds_mut = re.search(r'(c.\S+):p.',cds_mut_info).group(1)
			cds_mut = re.sub(r'(c\..*(?:ins|del|dup))(?:\d+|\w+)',r'\1',cds_mut)
		else:
			cds_mut = re.search(r'(c.\S+)',cds_mut_info).group(1)
			cds_mut = re.sub(r'(c\..*(?:ins|del|dup))(?:\d+|\w+)',r'\1',cds_mut)
		key="\t".join([chr,start,end,cds_mut])
		#print(key)
		given_dict[key]=[NM,mutation_type]

#将NM和mutation_type信息写入excel.
#Given_rb = xlrd.open_workbook('181207_TPNB500312_0040_AH75JJBGX7_20.xls')
Given_rb=QC_rb
Given_rs = Given_rb.sheet_by_name('all.must.given')
nrows = Given_rs.nrows
ncols = Given_rs.ncols
#Given_wb = copy(Given_rb)
Given_wb = QC_wb
Given_ws = Given_wb.get_sheet('all.must.given')
Given_ws.write(0,ncols,'NM')
Given_ws.write(0,ncols+1,'mutation_type')

for m,n in given_dict.items():
	print(m,n)
print('finish')
	
for i in range(1,nrows):
	chr,start,end,cds_mut = Given_rs.row_values(i)[4:8]
	cds_mut_NM = deepcopy_from_copy(cds_mut)	
	print('cds_mut:')
	gene = cds_mut.split(':')[0]
	print(gene)
	cds_mut = re.sub(r'(c\..*(?:ins|del))(?:\d+|[A-Z]+)',r'\1',cds_mut)
	start=str(int(start))
	end=str(int(end))
	key = '\t'.join([chr,start,end,cds_mut])
	NM,mutation_type = given_dict.get(key,['',''])
	if NM == '' and mutation_type == '':
		print('kong')
		print(cds_mut)
		cds_mut_2 = deepcopy_from_copy(cds_mut)
		cds_mut = re.search(r':(c\..*):p\..*',cds_mut).group(1) #TP53:NM_000546:exon7:c.708_709delCA:p.M237Vfs*2提取c.708_709delCA
		print(cds_mut)
		print('here?')
		print(gene)
		print(cds_mut,chr,start,end,gene,os.getcwd())
		mutation_type = annovar(cds_mut,chr,start,end,gene,os.getcwd())
		NM = cds_mut_NM.split(':')[1] 
		key = '\t'.join([chr,start,end,cds_mut])
		for each in given_dict.keys():
			if re.match(key,each) != None:
				NM,mutation_type = given_dict.get(each,['',''])
				print(NM,mutation_type)

	Given_ws.write(i,ncols,NM)
	Given_ws.write(i,ncols+1,mutation_type)
#print(given_dict)
Given_wb.save(OUT)

sample = re.search(r'(.*).brief.xls',CS_brief).group(1)
os.system("python %s/ResultTideUp.py --sample %s --xls-name %s "%(pipeline,sample,OUT))

