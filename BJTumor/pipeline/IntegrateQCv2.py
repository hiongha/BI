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
OUT = re.match(r'(.*)_\d+.xls',CS_QC).group(1)+'_cSMART_%s_bz.xls'%(sample_type_brief)

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
		cds_mut=re.search(r'(c.\S+):p.',cds_mut_info).group(1)
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
	cds_mut = re.sub(r'(c\..*(?:ins|del))(?:\d+|\w+)',r'\1',cds_mut)
	start=str(int(start))
	end=str(int(end))
	key = '\t'.join([chr,start,end,cds_mut])
	NM,mutation_type = given_dict.get(key,['',''])
	if NM == '' and mutation_type == '':
		#print(key)
		#print(cds_mut)
		#cds_mut = re.search(r'(c\..*(?:del|ins))(\d+|\w+):p.*',cds_mut).group(1) #TP53:NM_000546:exon7:c.708_709delCA:p.M237Vfs*2提取c.708_709delCA
		print('kong')
		print(cds_mut)
		cds_mut_2 = deepcopy_from_copy(cds_mut)
		cds_mut = re.search(r':(c\..*):p\..*',cds_mut).group(1) #TP53:NM_000546:exon7:c.708_709delCA:p.M237Vfs*2提取c.708_709delCA
		print(cds_mut)
		#if re.search(r'(del|ins)',cds_mut) != None:
		#	cds_mut = re.search(r':(c\..*(?:del|ins))(\d+|\w+):p\..*',cds_mut_2).group(1) #TP53:NM_000546:exon7:c.708_709delCA:p.M237Vfs*2提取c.708_709delCA
		if re.search(r'(del|ins)',cds_mut) != None:
			cds_mut = re.search(r':(c\..*(?:del|ins))(\d+|\w+):p\..*',cds_mut_2).group(1) #TP53:NM_000546:exon7:c.708_709delCA:p.M237Vfs*2提取c.708_709delCA
			print(cds_mut)
		key = '\t'.join([chr,start,end,cds_mut])
		for each in given_dict.keys():
			if re.match(key,each) != None:
				NM,mutation_type = given_dict.get(each,['',''])
				print(NM,mutation_type)

	Given_ws.write(i,ncols,NM)
	Given_ws.write(i,ncols+1,mutation_type)
#print(given_dict)
Given_wb.save(OUT)


