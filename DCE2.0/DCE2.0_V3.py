#!/usr/bin/python
#coding=utf-8
#xinghe <xinghe3223@berryoncology.com>

from __future__ import division
from datetime import datetime
import os
import re
import glob
import sys
import xlrd
import xlwt
from xlutils.copy import copy
from xlutils.filter import process,XLRDReader,XLWTWriter
from collections import OrderedDict
from copy import deepcopy
import argparse
from Bio import SeqIO
import subprocess


def read_mustgiven_add(tissuetype = ''):
	#file = '/share/Onc_Pro/cSMART/xinghe/dceChangeNM/DCE2.0_new_mustgiven_%s_add.v2.txt'%(tissuetype)
	file = '/share/public/pipeline/cSMART170503/supplement/DCE2.0_new_mustgiven_%s_add.v2.txt'%(tissuetype)
	#file = '/share/public/pipeline/cSMART170503/supplement/DCE2.0_new_mustgiven_%s_add.v2.txt'%(tissuetype)
	dic = {}
	dic2 = {}
	dic3 = {}
	dic4 = {}	
	f = open(file,'r')
	for line in f:
		gene,aa,exon,chr,start,end,cds_mut,NM,mutation_type,cds_mut_DCE2_format,snv_or_indel,aa_DCE2_format = line.strip().split('\t')	
		key = '_'.join([chr,start,end,gene,exon,cds_mut])
		dic[key] = mutation_type
		dic2[key] = cds_mut_DCE2_format
		dic3[key] = snv_or_indel
		dic4[key] = aa_DCE2_format				
	return dic,dic2,dic3,dic4			

def changeCDSAAFormat1(cds_mut,aa):
	try:
		aa = re.sub(r'(Exon.*)','.',aa)	
	except:
		pass	
	if re.search('>',cds_mut):
		try:
			alp1,alp2 = re.search(r'c.\d+_?\d+(?:\+?|\-?)(?:\d+)?([ATGC]*)>([ATGC]*)',cds_mut).groups()
			ref = alp1
			alt = alp2
		except:
			print('模式:')
			print(cds_mut,aa)
		if len(alp1) == 1 and len(alp2) ==1:
			return 'SNV',cds_mut,cds_mut,aa,aa
		elif len(alp1) == 0 and len(alp2) >=1 :
			new_cds_mut = re.sub(r'(c.\d+_?\d+\+?(?:\d+)?)[ATGC]*>([ATGC]*)',r'\1ins\2',cds_mut)
			try:
				new_aa = re.sub(r'>','ins',aa) #L747_A750>P
			except:
				new_aa = aa
			return 'Indel',cds_mut,new_cds_mut,aa,new_aa	
		else:
			new_cds_mut = re.sub(r'(c.\d+_?\d+\+?(?:\d+)?)[ATGC]*>([ATGC]*)',r'\1delins\2',cds_mut)
			try:
				new_aa = re.sub(r'>','delins',aa) #L747_A750>P
			except:
				new_aa = aa
			return 'Indel',cds_mut,new_cds_mut,aa,new_aa	
	elif re.search('del|ins|dup',cds_mut):
		return 'Indel',cds_mut,cds_mut,aa,aa

def changeCDSAAFormat(cds_mut):
	if re.search('>',cds_mut):
		try:
			alp1,alp2 = re.search(r'c.\d+_?\d+(?:\+?|\-?)(?:\d+)?([ATGC]*)>([ATGC]*)',cds_mut).groups()
			ref = alp1
			alt = alp2
		except:
			print('模式:')
			print(cds_mut)
		if len(alp1) == 1 and len(alp2) ==1:
			return 'SNV',cds_mut,cds_mut
		elif len(alp1) == 0 and len(alp2) >=1 :
			new_cds_mut = re.sub(r'(c.\d+_?\d+\+?(?:\d+)?)[ATGC]*>([ATGC]*)',r'\1delins\2',cds_mut)
			return 'Indel',cds_mut,new_cds_mut
		else:
			new_cds_mut = re.sub(r'(c.\d+_?\d+\+?(?:\d+)?)[ATGC]*>([ATGC]*)',r'\1delins\2',cds_mut)
			return 'Indel',cds_mut,new_cds_mut	
	elif re.search('del|ins|dup',cds_mut):
		return 'Indel',cds_mut,cds_mut

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

def TrRevSeq(seq):
	dic = {'A':'T','T':'A','G':'C','C':'G'}	
	seq = seq[::-1]
	tr_seq = []
	for i in seq:
		tr_seq.append(dic[i.upper()])
	return ''.join(tr_seq)		

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

def annovar(unfoundMutFunc,analysisPath):
	vcfList = []
	new_unfoundMutFunc = deepcopy(unfoundMutFunc)
	pre1 = analysisPath+'/list1';pre2 = analysisPath+'/list2'
	F4 = open('%s.annovar.hg19_multianno.txt'%(pre2),'a')
	F4.write('\t'.join(["Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene","Func.wgEncodeGencodeBasicV19","Gene.wgEncodeGencodeBasicV19","GeneDetail.wgEncodeGencodeBasicV19","ExonicFunc.wgEncodeGencodeBasicV19","AAChange.wgEncodeGencodeBasicV19","Otherinfo"])+'\n')
	for each in unfoundMutFunc:
		chr,start,end,gene,exon,cds_mut = unfoundMutFunc[each]['info']  
		ref,alt,retAnnovar = getrefalt(cds_mut = cds_mut,chr = chr,start = start,end = end,gene = gene) 
		vcfList.append(retAnnovar)
		new_unfoundMutFunc[each]['vcf'] = retAnnovar
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
		new_unfoundMutFunc[each]['mutFunc'] = mutFuncAnnovar
		os.system('mkdir -p %s/annovarTmp'%(analysisPath))
		#os.system('mv %s* %s* %s/annovarTmp'%(pre1,pre2,analysisPath))
		os.system('mv %s* %s* %s/annovarTmp'%(pre1,pre2,analysisPath))
	return new_unfoundMutFunc


def transvar(chr,start,end,cds_mut,gene,analysisPath):
	os.system('mkdir -p %s/transvarTmp'%(analysisPath))
	pre1 = analysisPath+'/temp.one.vcf';pre2 = analysisPath+'/transvarTmp/temp.all.vcf'
	ref,alt,retAnnovar = getrefalt(cds_mut = cds_mut,chr = chr,start = start,end = end,gene = gene)
	with open(pre1,'w') as F1:
		F1.write(retAnnovar)	
	with open(pre2,'a') as F2:
		F2.write(retAnnovar)
	p1 = subprocess.Popen('/share/public/software/Python-2.7.13/bin/transvar ganno --refseq --vcf %s > %s.result '%(pre1,pre1),shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	p1.wait()
	flag = False
	F4 = open(pre2+'.result','a')
	with open(pre1+'.result','r') as F3:
		for line in F3:
			if re.search(r'NM_\d+',line.strip().split('\t')[8]):
				NM_transvar = re.search(r'(NM_\d+)',line.strip().split('\t')[8]).group(1)
				print('NM_transvar:'+NM_transvar)
				if NM_transvar == gene2NM[gene]:
					flag = True 
					aa_transvar = line.strip().split('\t')[11].split('/')[2]
					cds_mut_transvar = line.strip().split('\t')[11].split('/')[1]
					F4.write(line)
					break	
	F4.close()							
	os.system('mv %s* %s/transvarTmp'%(pre1,analysisPath))
	if flag == True:
		return aa_transvar,cds_mut_transvar
	else:
		return 'no TranvarNM matches','no TranvarNM matches'	

def copy2(rb):
	wf = XLWTWriter()
	process(XLRDReader(rb,'foo.xls'),wf)
	return wf.output[0][0],wf.style_list

def parseBrief(file = ''):

	#key = '\t'.join([sample,chr,start,end,gene,exon,aa,cds_mut])
	dic = {}
	with open(file,'r') as F:
		F.readline()
		for line in F:
			Chr,Start,End,Ref,Alt,Func_refGene,Gene_refGene,ExonicFunc_refGene,AAChange_refGene,cosmic70,csmart_lung,Otherinfo = line.strip().split('\t',11)
			try:
				Gene = re.search(r'(.*)\(.*\)',Gene_refGene).group(1)
			except:
				pass
			if csmart_lung != '-':
				Gene,exon,aa,pos,cds_mut = re.search(r'Gene=(.*);exon=(.*);aa=(.*);Pos=(.*);CDS=(.+?);',csmart_lung).groups()
				try:
					aa = re.search(r'(Exon\d+_skipping)_?\d+',aa).group(1)
				except:
					pass
				key = '_'.join([Chr,Start,End,Gene,exon,cds_mut])
				if AAChange_refGene != '-':
					try:
						Gene1,nm1,exon1,cds_mut1,new_aa = AAChange_refGene.split(':') #ALK:NM_004304:exon25:c.3817A
					except:
						new_aa = '.'
				else:
					new_aa = '.'	
				if key not in dic.keys(): 
					dic[key] = {}
				dic[key]['aa'] = new_aa
				dic[key]['MutFunc'] = ExonicFunc_refGene
				
			elif AAChange_refGene != '-':
			#不同的Func_refGene,Gene_refGene的提取方式不同:
				try:
					Gene,nm,exon,cds_mut,aa = AAChange_refGene.split(':') #ALK:NM_004304:exon25:c.3817A
					try:
						exon = re.search(r'exon(\d+)',exon).group(1)
					except:
						pass
				except:
					Gene,nm,exon,cds_mut = AAChange_refGene.split(':') #ALK:NM_004304:exon25:c.3817A
					aa = '.'	
				#key = '_'.join([Chr,Start,End,Gene,exon,aa,cds_mut])
				key = '_'.join([Chr,Start,End,Gene,exon,cds_mut])
				if key not in dic.keys(): 
					dic[key] = {}
				dic[key]['aa'] = aa
				dic[key]['MutFunc'] = ExonicFunc_refGene
	for each in dic.keys():
		key1 = '--'
		if re.search(r'delins',each):
			key1 = re.sub(r'(.*delins).*',r'\1',each)
		elif re.search(r'del',each):
			key1 = re.sub(r'(.*del).*',r'\1',each)
		elif re.search(r'ins',each):
			key1 = re.sub(r'(.*ins).*',r'\1',each)
		elif re.search(r'dup',each):
			key1 = re.sub(r'(.*dup).*',r'\1',each)
		if key1 != '--':  
			dic[key1] = dic[each] 
	return dic

def gene2nm(file = ''):
	dic = {}
	dic2 = {}
	with open(file,'r') as F:
		for line in F:
			gene,NM,strand = line.strip().split('\t')
			dic[gene] = NM
			dic2[gene] = strand
	return dic,dic2	

def getBriefFile(sample='',analysisPath=''):
	briefList = []
	print(analysisPath)
	#briefList = glob.glob(analysisPath+'/*/*brief*.xls')
	briefList = glob.glob(analysisPath+'/*brief*.xls')
	if briefList == []:
		briefList = glob.glob(analysisPath+'/*brief*.xls')
	print(briefList)
	dic = {}
	for each in briefList:
		samp = os.path.basename(each).strip('brief.xls')
		dic[samp] = each
	print(dic)
	return dic

def code(List):
	newList = []
	for i in List:
		try:
			newList.append(i.encode())
		except:
			newList.append(i)	
	return newList

def toInt(List):
	newList = []
	for i in List:
		try:
			newList.append(str(int(float(i))))
		except:
			newList.append(i)
	return newList

#20190806
def getNewExon(newNM,startPos):
	exonP = subprocess.Popen('python /share/public/pipeline/cSMART170503/supplement/getcdsnum_V5.py --transid %s --posi %s'%(newNM,startPos),shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	for line in exonP.stdout:
		newExon = re.search(r'exon(\d+):',line.strip()).group(1)	
	return newExon					
	
def getSheetContent(mustgiven_mutFunc , workbook = '',sheetname='all.must.given'):
	
	rb = xlrd.open_workbook(workbook,formatting_info = True,on_demand = True)
	rs = rb.sheet_by_name(sheetname)
	nrow = rs.nrows
	ncol = rs.ncols
	title = rs.row_values(1)
	title = code(title)
	##sample	gene	aa	exon	chr	start	end	cds_mut	cosmic70	total	mutation	ratio	mutPosInTempates	nReads	tempID	cosmic	CIGAR	biaoji	comment	
	#sample	gene	exon	CDSmutation	proteinMutation	ReadDepth_Var/ReadDepth_Total	VarAlleleFrac	type	class	chr	MajorTranscriptID	MutFunc
	#sample	gene	exon	cds_mut	aa	mutation/total	ratio	type	class	chr	MajorTranscriptID	MutFunc
	#=====read information from brief file.
	samples = set()
	allSampleSet = set()
	for i in range(1,nrow):
		rowList = rs.row_values(i)
		allSampleSet.add(rowList[0])
		#if int(float(rowList[10])) == 0:
		#if rowList[18] == '-':
		if rowList[18] != 'positive':
			continue
		else:
			rowList = code(rowList)
			samples.add(rowList[0])
	briefInfos = {}
	for i in samples:
		samplebriefPath = briefPath[i]
		briefInfos[i] = parseBrief(file=samplebriefPath)
	
	newInfo = OrderedDict()
	unfoundMutFunc = OrderedDict()
	for i in range(1,nrow):
		rowList = rs.row_values(i)
		rowList = code(rowList)
		#if int(float(rowList[10])) == 0:
		if rowList[18] == '-':
			continue
		sample,gene,aa,exon,chr,start,end,cds_mut,cosmic70,total,mutation,ratio,mutPosInTempates,nReads,tempID,cosmic,CIGAR,biaoji,comment = rowList[::]
		exon,start,end,mutation = toInt( [exon,start,end,mutation])

		cds_mut_ori = deepcopy(cds_mut)
		aa_ori = deepcopy(aa)

		#对others位点或者整合的位点重新定义gene/NM/exon/cds_mut/aa等信息
		if aa_ori == 'others' or re.search(r':',cds_mut_ori):
			#TP53:NM_000546:exon7:c.747G>T:p.R249S
			#对于不是others位点,但是cds_mut列是这种形式的: EGFR:NM_005228:exon20:c.2312_2313insGAC:p.N771>KT
			try:
				gene,NM,exon,cds_mut,aa = re.search(r'(.+?):(NM_.+?):exon(\d+):(c\..+?):(p\..*)',cds_mut).groups()
			except:
				gene,NM,exon,cds_mut = re.search(r'(.+?):(NM_.+?):exon(\d+):(c\..+)',cds_mut).groups()
				aa = '-'
		elif re.search(r'Exon(\d+)_',aa_ori):
			pass
		else:
			aa = 'p.'+aa

		key = '_'.join([chr,start,end,gene,exon,cds_mut]) #这里没有加入氨基酸坐标
		newInfoKey = sample+'_'+key		
		newInfo[newInfoKey] = OrderedDict()
		#==========需要更新的变量======================================:
			#==需要从brief文件中提取
		#aa: exon_skipping和复杂位点的（带*数字的）需要更改aa信息：
		if re.search(r'Exon\d+_skipping',aa_ori,re.I):
			try:
				aa = briefInfos[sample][key]['aa']
			except:
				aa = '.' 
		#20190903 change fs and "*" situation
		elif re.search(r'(.*?)\*fs\*\d+$',cds_mut_ori): 
			aa = re.search(r'(.*?)\*fs\*\d+$',aa).group(1)+'Xfs'
		elif re.search(r'(.*)\*\d+$',cds_mut_ori):
			aa = re.search(r'(.*)\*\d+$',aa).group(1)
		elif re.search(r'(p..*)\*$',cds_mut_ori):
			aa = re.search(r'(p..*)\*$',aa).group(1)+'X'

		#for aa:氨基酸形式中含有'>'的需要去掉,如果存在于mustgiven表中就直接替换掉，如果不存在就用transvar注释一下。对于Exon_skipping位点没有氨基酸突变的现在改成skpping了。
		if key in mustgiven_aa_DCE2_format.keys():
			aa = mustgiven_aa_DCE2_format[key]
		elif re.search(r'>',aa):
			aa_transvar,cds_mut_transvar = transvar(chr=chr,start = start,end = end,cds_mut = cds_mut,gene = gene,analysisPath=current_path)	
			#aa_transvar = transvar(chr=chr,start = start,end = end,cds_mut = cds_mut,gene = gene,analysisPath=analysisPath)	
			if aa_transvar == 'no TranvarNM matches':
				pass
			else:
				aa = aa_transvar
		elif (re.search(r'skipping',aa_ori) or re.search(r'skipping',cds_mut_ori)) and aa == '.':
			aa = 'splicing'		
				
		#mutFunc: others位点和合并的位点需要知道mutFunc信息:
		MutFunc = '--'
		if key in mustgiven_mutFunc.keys():
			MutFunc = mustgiven_mutFunc[key]
		elif key in briefInfos[sample].keys():
			MutFunc = briefInfos[sample][key]['MutFunc']
		elif re.search(r'fs\*',cds_mut_ori): 
			#aa = re.search(r'(p\..*fs)\*.*',aa).group(1)
			if re.search(r'del',cds_mut) and re.search(r'ins',cds_mut):
				MutFunc = 'frameshift substitution'
			elif re.search(r'del',cds_mut) and (not re.search(r'ins',cds_mut)):
				MutFunc = 'frameshift deletion'
			if re.search(r'dup',cds_mut) or (re.search(r'ins',cds_mut) and not re.search(r'del',cds_mut)):
				MutFunc = 'frameshift insertion'
#		elif key1 in briefInfos[sample].keys():
		elif re.search(r'(.*)(delins|del|ins|dup)(\d+|[ATGC]+)',key):
			key1 = re.sub(r'(.*)(delins|del|ins|dup)(\d+|[ATGC]+)',r'\1\2',key)
			#print(key1)
			if key1 in briefInfos[sample].keys():
				MutFunc = briefInfos[sample][key1]['MutFunc']
			else:
				if newInfoKey not in unfoundMutFunc.keys():
					unfoundMutFunc[newInfoKey] = OrderedDict() 
				unfoundMutFunc[newInfoKey]['info'] = [chr,start,end,gene,exon,cds_mut]
				print('没找到这个others位点的aa和MutFunc信息:')
				print(rowList[0:8])
				print(key,sample,analysisPath,inputfilexls)
				print(briefInfos[sample])
				MutFunc = '--'
	#			sys.exit('\n')
		else:
			if newInfoKey not in unfoundMutFunc.keys():
				unfoundMutFunc[newInfoKey] = OrderedDict() 
			unfoundMutFunc[newInfoKey]['info'] = [chr,start,end,gene,exon,cds_mut]
			print('没找到这个others位点的aa和MutFunc信息:')
			print(rowList[0:8])
			print(key,sample,analysisPath,inputfilexls)
			print(briefInfos[sample])
			MutFunc = '--'
	#		sys.exit('\n')
	
			#===不需要从brief文件中提取的变量
		#for mutation/total
		mutation_total = str(int(float(mutation)))+'/'+str(int(float(total)))
		#for type:
		if key in mustgiven_cds_mut_DCE2_format.keys():
			cds_mut_after= mustgiven_cds_mut_DCE2_format[key]
			cds_mut_before = cds_mut_ori
			muttype = mustgiven_mutType[key]
		else:		
			muttype,cds_mut_before,cds_mut_after= changeCDSAAFormat(cds_mut)
		#20190806: to adapt changed NM ID:
		if re.search(r'other',aa_ori) and gene in ['KRAS','MET','HRAS','AKT1']:
			aa_transvar,cds_mut_transvar = transvar(chr=chr,start = start,end = end,cds_mut = cds_mut,gene = gene,analysisPath=current_path)
			new_exon = getNewExon(gene2NM[gene],int(start))
			exon = new_exon
			print('new_exon:'+str(new_exon))	
			if aa_transvar == 'no TranvarNM matches':
				cds_mut_after = '%s was not matched in transvar result'%(gene2NM[gene])
			else:
				cds_mut_after = cds_mut_transvar
				aa = aa_transvar
				#20190903 change fs and "*" situation
				if re.search(r'(.*?)\*fs\*\d+$',aa_transvar): 
					aa = re.search(r'(.*?)\*fs\*\d+$',aa_transvar).group(1)+'Xfs'
				elif re.search(r'(.*)\*\d+$',aa_transvar):
					aa = re.search(r'(.*)\*\d+$',aa_transvar).group(1)
				elif re.search(r'(p..*)\*$',aa_transvar):
					aa = re.search(r'(p..*)\*$',aa_transvar).group(1)+'X'
				if (re.search(r'skipping',aa_ori) or re.search(r'skipping',cds_mut_ori)) and aa == '.':
					aa = 'splicing'	
		#for class:
		classa = 'somatic'
		#for MajorTranscriptID
		try:
			NM = gene2NM[gene]
		except:
			NM = 'not in all_gene.txt'
	
		#fortest: 
		#newInfo[newInfoKey].update({'sample':sample,'gene':gene,'exon':exon,'cds_mut':cds_mut_ori,'cds_mut_after':cds_mut_after,'aa':aa_ori,'aa_after':aa,'mutation_total':mutation_total,'ratio':ratio,'type':muttype,'class':classa,'chr':chr,'MajorTrans':NM,'MutFunc':MutFunc})
		#newInfo[newInfoKey].update({'sample':sample,'gene':gene,'exon':exon,'cds_mut_after':cds_mut_after,'aa_after':aa,'mutation_total':mutation_total,'ratio':ratio,'type':muttype,'class':classa,'chr':chr,'MajorTrans':NM,'MutFunc':MutFunc})
		#20190813:
		newInfo[newInfoKey].update({'sample':sample,'gene':gene,'exon':exon,'cds_mut_after':cds_mut_after,'aa_after':aa,'mutation_total':mutation_total,'ratio':ratio,'type':muttype,'class':classa,'chr':chr,'MajorTrans':NM,'MutFunc':MutFunc,'cds_mut_before':cds_mut_ori,'aa_before':aa_ori})
#sample	gene	exon	CDSmutation	proteinMutation	ReadDepth_Var/ReadDepth_Total	VarAlleleFrac	type	class	chr	MajorTranscriptID	MutFunc
	return newInfo,samples,allSampleSet,unfoundMutFunc


def update_newInfo(newInfo,unfoundMutFunc):
	for each in unfoundMutFunc.keys():
		print(each,unfoundMutFunc[each]['mutFunc'])
		newInfo[each]['MutFunc'] = unfoundMutFunc[each]['mutFunc']	 

def check(newInfo):
	print(newInfo.keys())
	for each in newInfo.keys():
		if newInfo[each]['MutFunc'] == '--':
			sys.exit('MutFunc information of %s was not found!'%s(each))


def writeToExcel(newInfo,samples,allSampleSet,workbook,sheetname='SnvIndel',output = 'try.xls'):

	rb = xlrd.open_workbook(workbook,formatting_info = True,on_demand=True)
	wb = copy(rb)
	ws = wb.add_sheet(sheetname)
#	rb.release_resources()
	#sample	gene	exon	CDSmutation	proteinMutation	ReadDepth_Var/ReadDepth_Total	VarAlleleFrac	type	class	chr	MajorTranscriptID	MutFunc
	i = 0
	p = 0
	#for test: 
	#for t in ["sample","gene","exon","CDSmutation",'cds_mut_after',"proteinMutation",'aa_after',"ReadDepth_Var/ReadDepth_Total","VarAlleleFrac","type","class","chr","MajorTranscriptID","MutFunc"]:
	for t in ["sample","gene","exon","CDSmutation","proteinMutation","ReadDepth_Var/ReadDepth_Total","VarAlleleFrac","type","class","chr","MajorTranscriptID","MutFunc","cds_mut_before","aa_before"]:
		ws.write(i,p,t)
		p += 1
	i += 1	
	for each in newInfo.keys():
#		newList = []
		j = 0
		#for test:
		#for col in ['sample','gene','exon','cds_mut','cds_mut_after','aa','aa_after','mutation_total','ratio','type','class','chr','MajorTrans','MutFunc']:
		#for col in ['sample','gene','exon','cds_mut_after','aa_after','mutation_total','ratio','type','class','chr','MajorTrans','MutFunc']:
		#20190813:
		for col in ['sample','gene','exon','cds_mut_after','aa_after','mutation_total','ratio','type','class','chr','MajorTrans','MutFunc','cds_mut_before','aa_before']:
			ws.write(i,j,newInfo[each][col])
			j += 1
		i += 1
	negSet = allSampleSet-samples
	print(allSampleSet)
	print(samples)
	for neg in negSet:
		ws.write(i,0,neg)
		ws.write(i,1,'Negative')
		for j in range(2,14):
			ws.write(i,j,'-')
			j += 1	
		i += 1
	wb.save(output)

def CreateNewVanSheet(workbook,sheetname='SnvIndel',output = 'try.xls'):
	rb = xlrd.open_workbook(workbook,formatting_info = True,on_demand=True)
	wb = copy(rb)
	ws = wb.add_sheet(sheetname)
	wb.save(output)

#---------add Sheet "time"------------------------------------------
def getSampleCfgFile(report_of_sample):
	sample_path = os.path.dirname(report_of_sample)
	sample_cfg = glob.glob(sample_path + '/sample.cfg.*')[0]
	#finish_file = glob.glob(sample_path + '/finish.txt')[0]
	finish_file = glob.glob(sample_path + '/final_result.xls')[0]
	samples = []
	with open(sample_cfg,'r') as F:
		for eachline in F:
			sample = eachline.strip().split('\t')[3]
			samples.append(sample)
	return samples,sample_cfg,finish_file

def getFileTime(startfile,endfile,final_report):
	un_starttime = os.path.getmtime(startfile)
	un_endtime = os.path.getmtime(endfile)
	starttime = datetime.fromtimestamp(int(un_starttime)).strftime('%Y/%m/%d,%H:%M:%S')
	endtime = datetime.fromtimestamp(int(un_endtime)).strftime('%Y/%m/%d,%H:%M:%S')
	delta_time = round((un_endtime - un_starttime)/60/60,2)
	un_report_complete_time = os.path.getmtime(final_report)
	report_complete_time = datetime.fromtimestamp(int(un_report_complete_time)).strftime('%Y/%m/%d,%H:%M:%S')
	return (starttime,endtime,delta_time,report_complete_time)

def addTimeSheet(inputfile,outputfile):
	#rb_file = xlrd.open_workbook(inputfile,formatting_info = True,on_demand = True)
	rb_file = xlrd.open_workbook(outputfile,formatting_info = True,on_demand = True)
	#wb_file = xlrd.open_workbook(outputfile,formatting_info = True,on_demand = True)
	wb_file = copy(rb_file)
	ws_time = wb_file.add_sheet('time')
	title_time = ['Sample','start_time','end_time','time_used(hour)','report_complete_time']
	for i,j in enumerate(title_time):
		ws_time.write(0,i,j)
	inputfile_chai = inputfile.split(',')	
	for ith,each in enumerate(inputfile_chai):
		print(each)
		samples,sample_cfg,finish_file = getSampleCfgFile(each)
		print('samples:')
		print(samples)
		for sample in samples:
			print('sample:')
			print(sample)
			starttime,endtime,delta_time,report_complete_time = getFileTime(sample_cfg,finish_file,outputfile)
			info = [sample,starttime,endtime,delta_time,report_complete_time]
			for m,n in enumerate(info):
				ws_time.write(ith+1,m,n)
			ith=ith+1
	wb_file.save(outputfile)

#------------add Sheet "cis_trans"----------
def add_cis_trans_sheet(outputfile):
	rb_file = xlrd.open_workbook(outputfile,formatting_info=True,on_demand=True)
	wb_file = copy(rb_file)
	ws_cis_trans = wb_file.add_sheet('egfr_cis_trans')
	ith = 0 
	cis_trans_file = os.path.dirname(outputfile)+'/before/cis_tran.xls'
	with open(cis_trans_file) as F:
		for line in F:
			jth = 0
			for each in line.strip().split('\t'):
				ws_cis_trans.write(ith,jth,each)
				jth +=1
			ith +=1	
	wb_file.save(outputfile)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(prog = '',description = 'add DCE2.0 sheet to cSMART')
	parser.add_argument("-p","--analysis-path", help = 'analysis path',required = True)
	parser.add_argument("-i","--inputfile", help = 'input file', required = True)
	parser.add_argument("-o","--outputfile", help = 'output file name',default='try.xls')
	parser.add_argument("-s","--sample-source", help = 'sample type: lung,crc_2.1,not contain SCLC,LYM17,TM21,',required=True)
	args = parser.parse_args()
	fa = '/share/public/database/genome/Homo_sapiens/UCSC-hg19/download/hg19.chr.fa'
	analysisPath = args.analysis_path
#	for test:
	current_path = os.getcwd()
	inputfilexls = args.inputfile
	outputfilexls = args.outputfile
	sample_source = args.sample_source

	if sample_source != 'lung' and sample_source != 'crc_2.1':
		print('this is not lung or crc_2.1')
		CreateNewVanSheet(workbook = inputfilexls,sheetname='SnvIndel',output = outputfilexls)
		add_cis_trans_sheet(outputfile = outputfilexls)
		addTimeSheet(inputfile = inputfilexls ,outputfile = outputfilexls)
	else:		
		gene2NM,gene2strand = gene2nm(file = '/share/public/pipeline/cSMART170503/supplement/DCE2.0_all_gene_strand.v2.txt')
		#gene2NM,gene2strand = gene2nm(file = '/share/public/pipeline/cSMART170503/supplement/DCE2.0_all_gene_strand.v2.txt')
		mustgiven_mutFunc,mustgiven_cds_mut_DCE2_format,mustgiven_mutType,mustgiven_aa_DCE2_format = read_mustgiven_add(tissuetype = sample_source)
		#print(gene2NM)
		#briefPath = getBriefFile(sample='',analysisPath='/share/Oncology/production/cSMART/cSMART/analysis/xinghe/DCE/test1')		
		briefPath = getBriefFile(sample='',analysisPath=analysisPath)		
		#newInfo,samples,allSampleSet = getSheetContent(workbook = '190514_TPNB500289_0124_AHK23GAFXY_LC.xls')
		newInfo,samples,allSampleSet,unfoundMutFunc = getSheetContent(mustgiven_mutFunc = mustgiven_mutFunc, workbook = inputfilexls)
		print(newInfo)
		if unfoundMutFunc != OrderedDict():
			print('unfoundMutFunc:')
			print(unfoundMutFunc)
			current_path = os.getcwd()
			#current_path = analysisPath
			
			new_unfoundMutFunc = annovar(unfoundMutFunc,current_path)
			#new_unfoundMutFunc = annovar(unfoundMutFunc,analysisPath)
			print('new_unfoundMutFunc')
			print(new_unfoundMutFunc)
			update_newInfo(newInfo = newInfo,unfoundMutFunc = new_unfoundMutFunc)
	#		print('\n')
		check(newInfo = newInfo)
		writeToExcel(newInfo = newInfo,samples = samples,allSampleSet=allSampleSet,workbook = inputfilexls,output=outputfilexls)
		print(inputfilexls)
		print(outputfilexls)
		add_cis_trans_sheet(outputfile = outputfilexls)
		addTimeSheet(inputfile = inputfilexls ,outputfile = outputfilexls)
		#writeToExcel(newInfo = newInfo,samples = samples,allSampleSet=allSampleSet,workbook = inputfilexls,output=outputfilexls)


