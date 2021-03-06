#!/usr/bin/python
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>
'''
csmart数据转换成北肿要求的csv文件中的格式.
'''
from __future__ import division
import copy
import re
from Bio import SeqIO
import collections
import os
import argparse
import xlrd
import csv

parser = argparse.ArgumentParser(prog = '', description = 'csmart数据转换成北肿要求的csv文件格式')
parser.add_argument('--brief-file', help = 'csmart中的sample.brief.xls文件', required = True)
parser.add_argument('--must-given', help = 'csmart中的sample.must.given.xls文件', required = True)
parser.add_argument('--final-report', help = 'csmart中的最终报告,如190111_NS500511_0178_AH3T5GAFXY_1.xls', required = True)
parser.add_argument('--out-put', help = '输出文件sample.mutations.csv', required = True)
parser.add_argument('--type', help = '样本类型,如lung/crc_2.1/...', required = True)
args=parser.parse_args()
print args

#--------------初始化变量及函数定义------------------------------------------------
pipeline = '/share/work1/xinghe/proc/BZQC/pipeline'
root_dir = os.path.abspath('./')
lung = pipeline + '/config/hg19_csmart_%s.txt'%(args.type)
NMfile = pipeline + '/config/all_gene.txt'
fa = '/share/public/database/genome/Homo_sapiens/UCSC-hg19/download/hg19.chr.fa'
strand_file = pipeline + '/config/all_gene_strand.txt'
brief = root_dir + "/" + args.brief_file
sample_must_given_xls = root_dir + '/'+ args.must_given
mutcsv_output = root_dir + '/' + args.out_put
final_report = root_dir + '/' + args.final_report

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

def findpos(fa ,chr, start):
	'''
	寻找输入位点start的前一个位点的碱基
	'''
	with open(fa, 'r') as fa:
		for record in SeqIO.parse(fa,'fasta'):
			if record.id == chr:
				seq = record.seq
				return seq[int(start)-1-1]

#----------氨基酸对应表-------------------------------------------------------------
sgl_abbr = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
tri_abbr = ["Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"]
aa_abbr_dict = dict(zip(sgl_abbr,tri_abbr))

#------------从报告中提取positive位点-----------------------------------------------
def positive_find(final_report):
	'''
	从xls报告中提取positive位点.
	'''
	rw = xlrd.open_workbook(final_report)
	rs = rw.sheet_by_name('all.must.given')
	nrows = rs.nrows
	ncols = rs.ncols
	positive_site_list = []
	for i in range(1,nrows):
		line = rs.row_values(i)[0:8]
		positive_flag = rs.row_values(i)[18]
		if positive_flag == 'positive':
			new_line = []
			for each in line: 
				if isinstance(each,float) == True:
					new_line.append(str(int(each)))
				else:
					new_line.append(str(each))
			new_line = '\t'.join(new_line)
			positive_site_list.append(new_line)
		else:
			pass
	return positive_site_list

positive_site_list = positive_find(final_report)
#----------------------格式处理-----------------------------------------------------------
#基因正负链信息.
strandF = open(strand_file,'r')
strandDict={}
for line in strandF:
	gene_name = line.strip().split()[0]
	strand = line.strip().split()[2]
	strandDict[gene_name] = strand

#根据必报点文件获取Ref和Alt序列信息
#lung='/share/public/pipeline/cSMART170503/third-party/annovar/humandb/hg19_csmart_lung.txt'
vcfF=open('list1.vcf','w')
mutcsvF=open(mutcsv_output,'w')
mutcsvF.write("Gene\tGene.ID\tAA.Change\tChr.start\tChr.end\tRef\tAlt\tDP(ref:alt)\tStrand\tAF\tFunc\tdbSNP\tCosmicID\tCosmic.Occurence\t1000g\t1000gEAS\tExAC\tSIFT.score\tSIFT.pred\tPolyPhen.score\tPolyPhen.pred\tCLIN_SIG\tHGVS\n")
lungF=open(lung,'r')
info_dict={}
vcf_csv_dict=collections.OrderedDict()
for line in lungF:
	info=line.strip().split('\t')
	chr='chr'+info[0]
	start,end,ref,alt,mut=info[1:]
	cds_mut=re.search(r'CDS=(.*?)(?:;)',mut).group(1)
	key=[chr,start,end,cds_mut]
	key='_'.join(key)
	value=[ref,alt]
	info_dict[key]=value

#gene与NM号对应.
#NMfile='/share/public/pipeline/170711_Amplicon_programV2/all_gene.txt'
NM_dict={}
NM_F=open(NMfile,'r')
for line in NM_F:
	gene,NM=line.strip().split('\t')
	NM_dict[gene]=NM	

#brief文件信息提取
#brief="/share/work1/xinghe/proc/BZQC/Z18L06516/Z18L06516.brief.xls"
brieF=open(brief,'r')
brieF.readline()
info_dict_brief={}
for line in brieF:
	list=line.strip().split('\t')
	chr,start,end,ref,alt=list[0:5]
	cds_mut=''
	cds_mut_info=list[8]  #ALK:NM_004304:exon25:c.3834C>T:p.Y1278Y  
	if cds_mut_info=='-':
		continue
	elif re.search(r':(c\..*):',cds_mut_info):
		cds_mut=re.search(r':(c\..*):',cds_mut_info).group(1)
	elif re.search(r':(c\..*)',cds_mut_info):
		cds_mut=re.search(r':(c\..*)',cds_mut_info).group(1)
	key=[chr,start,end,cds_mut]
	key='_'.join(key)
	value=[ref,alt]
	info_dict_brief[key]=value

#格式处理
reverse_list={'A':'T','T':'A','G':'C','C':'G'}
allGivenF=open(sample_must_given_xls,'r')
title=allGivenF.readline().strip().split('\t')
for line in allGivenF:
	positive_flag = '\t'.join(line.strip().split('\t')[0:8])
	if positive_flag in positive_site_list:
		positive_flag = True
		print('flagdaodishisha')
	else:
		positive_flag = False
		print('flagdaodishishaxxxxx')
	print(line.strip())
	line = line.replace('\tunique_cons_test','')
	[sample,gene,aa,exon,chr,start,end,cds_mut,cosmic70,total,mutation,ratio,mutPosInTempates,nReads,tempID,cosmic,CIGAR]=line.strip().split('\t')
	##AA.stat	#YES|T
	aa_ori = copy.deepcopy(aa)
	Gene=gene	#SMAD4
	NM=NM_dict.get(Gene,"")
	Gene_ID=gene+':'+NM+':exon'+exon	#SMAD4:NM_005359.5:exon9,添加NM编号
	if re.search(r'others|Exon',aa,re.I)!=None or re.search(r':',cds_mut):
		AA_Change=cds_mut	#p.R361S (c.C1081A)
	else:
		AA_Change='p.'+aa+'(%s)'%(cds_mut) #cds_mut还需要转换一下格式.
	
	Chr_start=chr+':'+start	#chr18:48591918
	Chr_end=chr+':'+end	#chr18:48591918

	cds_mut_ori=copy.deepcopy(cds_mut)	
	cds_mut=re.split('/|:',cds_mut)
	cds_mut=[j for j in cds_mut if re.match('c.',j)]
	Ref=''
	Alt=''
	for each in cds_mut:
		key='_'.join([chr,start,end,each])
		if re.search(r'ins\d+',key)!=None:
			key=re.sub(r'ins\d+','ins',key).strip()
			r=r'%s'%key
			for x in info_dict.keys():
				if re.match(key,x)!=None:
					Ref,Alt=info_dict[x]
		else:
			Ref=info_dict.get(key,['',''])[0]	
			Alt=info_dict.get(key,['',''])[1]
		if int(start)>int(end): 
			Ref=Ref[::-1]	#C #反向互补
			Alt=Alt[::-1]	#A #反向互补
			Ref=''.join([reverse_list[each] for each in Ref])
			Alt=''.join([reverse_list[each] for each in Alt])
		if Ref!='' or Alt!='':
			break
		else:
			continue
	if Ref=='' and Alt=='':
		exon_new=re.search(r'(exon\d+)',cds_mut_ori).group(1)
		Gene_ID=gene+':'+NM+':'+exon_new	#SMAD4:NM_005359.5:exon9,添加NM编号
		cds_tmp,aa_tmp=re.search(r'(c\..*):(p\..*)',cds_mut_ori).groups()
		AA_Change=aa_tmp+"("+cds_tmp+")"
		
		cds_tmp=re.sub(r'(del|ins|dup)(\w+|\d+)',r'\1',cds_tmp)
		search='_'.join([chr,start,end,cds_tmp])
		for item in info_dict_brief.keys():
			if re.search(search,item)!=None:
				Ref,Alt=info_dict_brief.get(item,'empty')	
	DP=total+"|"+mutation	#5466(3484:1982)|2813(1815:998)
	Strand = strandDict[Gene]
	if int(total) == 0:
		AF = 0
	else:
		AF = int(mutation)/int(total)	#36.26%|35.48%	
	#for HGVS
	hgvs_nm = re.search(r'(NM_\d+)',Gene_ID).group(1)
	hgvs_gene = Gene	
	hgvs_cds_mut = re.search(r'(c\..+?)(?:$|:p.)',cds_mut_ori).group(1)
	#for hgvs_aa
	#hgvs_aa = copy.deepcopy(aa_ori)
	if re.search('Exon',aa_ori,re.I) != None:
		hgvs_aa = ''
	elif re.search('others',aa_ori,re.I) != None or re.search(r':',cds_mut_ori):
		hgvs_aa_others_site  = re.search(r':p\.(.+)',cds_mut_ori).group(1)
		aa_set = set(re.findall('[A-Z]',hgvs_aa_others_site))
		for each_aa in aa_set:
			hgvs_aa_others_site = re.sub(each_aa, aa_abbr_dict[each_aa], hgvs_aa_others_site)
		hgvs_aa = '(' + hgvs_aa_others_site + ')'
	else:
		aa_set = set(re.findall('[A-Z]',aa_ori))
		for each_aa in aa_set:
			aa_ori = re.sub(each_aa, aa_abbr_dict[each_aa], aa_ori)
		hgvs_aa = '(' + aa_ori + ')'

	HGVS =	hgvs_nm + '('+hgvs_gene+')'+':'+ hgvs_cds_mut +hgvs_aa
	print HGVS
	#print hgvs_cds_mut,'\t',hgvs_aa

	if positive_flag == True:
		print('daodishishei:')
		print(hgvs_cds_mut,chr,start,end,Gene)
#		ref_func,alt_func,retAnnovar = getrefalt(cds_mut,chr,start,end,Gene)
		ref_func,alt_func,retAnnovar = getrefalt(hgvs_cds_mut,chr,start,end,Gene)
#		pos=findpos(fa = fa, chr = chr, start = start)
#		ref_tmp = '' if Ref == '-' else Ref
#		alt_tmp = '' if Alt == '-' else Alt
#		vcfF.write('\t'.join([chr,str(int(start)-1),'.',pos+ref_tmp,pos+alt_tmp,'.','.','.'])+'\n')	
#		vcf_key ='\t'.join([chr,str(int(start)-1),'.',pos+ref_tmp,pos+alt_tmp,'.','.','.'])	
		vcfF.write(retAnnovar)
		vcf_key = retAnnovar.strip()
		csv_value = Gene+"\t"+Gene_ID+"\t"+AA_Change+"\t"+Chr_start+"\t"+Chr_end+"\t"+Ref+"\t"+Alt+"\t"+DP+"\t"+Strand+"\t"+str(AF)
		vcf_csv_dict[vcf_key] = (csv_value,HGVS)

vcfF.close()

#注释并解析注释结果
try:
	os.system('sh %s/annotation.sh' %(pipeline))
	annF = open('list1.annovar.hg19_multianno.txt','r')
	title=annF.readline().strip().split('\t')

	for line in annF:
		for each in vcf_csv_dict.keys():
			if re.search(each,line) != None:
				info=line.strip().split('\t')
				anninfo=dict(zip(title,info))			
				#[Func,dbSNP,CosmicID,Cosmic_Occurence,1000g,1000gEAS,ExAC,SIFT.score,SIFT.pred,PolyPhen_score,PolyPhen_pred,CLIN_SIG,HGVS] = [ anninfo.get(i,'empty') for i in ['ExonicFunc.refGene','avsnp147','cosmic70','cosmic70','1000g2015aug_all','1000g2015aug_eas','ExAC_ALL,SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','clinvar_20170130']]
				anninfo_used = [ anninfo.get(i,'empty') for i in ['ExonicFunc.refGene','avsnp147','cosmic70','cosmic70','1000g2015aug_all','1000g2015aug_eas','ExAC_ALL','SIFT_score','SIFT_pred','Polyphen2_HDIV_score','Polyphen2_HDIV_pred','clinvar_20170130']]
			#Func,dbSNP,CosmicID,Cosmic.Occurence,1000g,1000gEAS,ExAC,SIFT.score,SIFT.pred,PolyPhen.score,PolyPhen.pred,CLIN_SIG,HGVS
				cosmic_tmp = anninfo_used[2]
				if cosmic_tmp != '.':
					cosmic_id,cosmic_occur = re.search(r'ID=(.*);OCCURENCE=(.*)',cosmic_tmp).groups(0)
					cosmic_id = re.sub(r',',r'|',cosmic_id)
					cosmic_occur = re.sub(r',',r'|',cosmic_occur)
					anninfo_used[2] = cosmic_id
					anninfo_used[3] = cosmic_occur

				mutcsvF.write( vcf_csv_dict[each][0]+'\t'+'\t'.join(anninfo_used)+'\t'+vcf_csv_dict[each][1]+'\n')
except:
	pass
finally:
	mutcsvF.close()

mutcsv_output_tmp = mutcsv_output + '_tmp'
os.system('mv %s %s'%(mutcsv_output,mutcsv_output_tmp))

#将mut.csv文件转换成真正的csv格式
def txt2csv(txt,csv_file, delimiter = '\t'):
	in_txt = csv.reader(open(txt,'rb'), delimiter = delimiter)
	out_csv = csv.writer(open(csv_file,'wb'))
	out_csv.writerows(in_txt)

txt2csv(mutcsv_output_tmp,mutcsv_output)


