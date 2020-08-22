#!/usr/bin/python
#coding=utf-8
#sample=$1
#Author: xinghe
#Contact: xingh3223@berryoncology.com
import os
import sys
import glob
import subprocess
import re

pipeline = '/share/work1/xinghe/proc/BZQC/pipeline'
root_dir = os.getcwd()
info = 'info'
info_txt = 'info.txt'
#==========整理待合并位点信息文件格式,使其标准化======================================
with open(root_dir+'/'+info,'r') as infoF:
	sample = infoF.readline().strip()
	with open(root_dir+'/'+info_txt,'w') as infoF1:
		infoF1.write(sample+'\n')
print('待分析的样品为: {sample}'.format(**locals()))

#=============寻找并链接数据到当前路径=============================================
#datapath1 = "/share/Oncology/production/cSMART/cSMART/analysis/sample_cfg_tmp/sample.cfg.*"
#datapath2 = "/share/Oncology/production/cSMART/cSMART/analysis"
#datapath3 = "/share/Oncology/production/cSMART/cSMART/analysis/sample_cfg_tmp/analysis_done_sample/sample.cfg.*"
#datapath1 = "/share/Onc_Pro/cSMART/sample_cfg_tmp/sample.cfg.*"
#datapath2 = "/share/Onc_Pro/cSMART"
#datapath3 = "/share/Onc_Pro/cSMART/sample_cfg_tmp/analysis_done_sample/sample.cfg.*"
datapath1 = "/share/Onc_SmallPanel/cSMART/sample_cfg_tmp/sample.cfg.*"
datapath2 = "/share/Onc_SmallPanel/cSMART"
datapath3 = "/share/Onc_SmallPanel/cSMART/sample_cfg_tmp/analysis_done_sample/sample.cfg.*"
primer_CRC = ['CRC_TM21_V1.3']
primer_LC = ['LC9_V4.10','LC9_V4.12','LC6_V5.4']
primer_SCLC = ['SCLC_V1.3']
primer_TM21 = ['TM21_V1.3']
primer_sfx = ['CRC','LC','SCLC','TM21']
org = ['crc_2.1','lung','SCLC','TM21']
primer_sfx2org = dict(zip(primer_sfx,org))

sample_cfgs = subprocess.Popen('grep %s %s'%(sample,datapath1),stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE,shell=True,universal_newlines=True)
ith = 0
for each in sample_cfgs.stdout:
	ith +=1 
	cfg = each.strip().split('\t')
	FC_short = cfg[0].split(':')[1]
	primer = cfg[2]
	for sfx in primer_sfx:
		prm_tmp = 'primer_'+sfx
		if primer in locals()[prm_tmp]:
			primer_prefix = sfx
			primer_prefix_org = primer_sfx2org[primer_prefix]
	print primer_prefix,primer_prefix_org	
	sample_type = cfg[5]
	all_dir = datapath2+'/'+'*'+FC_short+'*'+'/'+'*/sample.cfg.*'
	sample_cfg_in_analysis_dir = subprocess.Popen('grep -l %s %s'%(sample,all_dir),stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE,shell=True,universal_newlines=True)
	sample_cfg = []	
	for ipath in sample_cfg_in_analysis_dir.stdout:
		ipath=ipath.strip()
		sample_cfg.append(ipath)
		os.system('ln -fs %s'%(os.path.dirname(ipath)))
		csmart_analysis_path = os.path.dirname(ipath)
		FC = os.path.basename(os.path.dirname(ipath))
		R1 = os.path.basename(glob.glob(csmart_analysis_path+'/*.R1.clean.fastq.gz')[0])
		R2 = os.path.basename(glob.glob(csmart_analysis_path+'/*.R2.clean.fastq.gz')[0])
		report_xls_name = re.sub(r'(.*_.*_.*_.*)_.*',r'\1_cSMART_%s_bz.xls'%primer_prefix,os.path.basename(csmart_analysis_path)) #190121_NS500511_0183_AH3LHTAFXY_cSMART_LC_bz.xls 
		report_zip_name = re.sub(r'(.*_.*_.*_.*)_.*',r'\1_cSMART_%s.zip'%primer_prefix,os.path.basename(csmart_analysis_path)) #190121_NS500511_0183_AH3LHTAFXY_cSMART_LC.zip
		with open(root_dir+'/mail_tmp.sh','w') as F1:
			F1.write( "/share/public/software/Python-2.7.13/bin/python %s/mail_bz_check.py -s %s -f %s,%s -w %s \n" %(pipeline,sample,report_xls_name,report_zip_name,os.getcwd()))
		with open(root_dir+'/check_jobs.sh','w') as F2:
			F2.write("source /home/xinghe/.bashrc;/share/public/software/Python-2.7.13/bin/python %s/check.py -u %s,csmartpro -s runbz1.sh\n"%(pipeline,os.getlogin()))
		with open(root_dir+'/runbz1.sh','w') as F3:
			F3.write( "samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam |"%(FC,sample)) 
			F3.write("awk '{print $1"+'"\\t"$3"\\t"$4"\\t"$6"\\t"$10"\\t"$(NF-2)}'+"'")
			F3.write("| sort -k1,1 |perl %s/area.pl > reformed_bam_pos.txt\n" %(pipeline)) #该步骤分析特别耗内存和存储，建议文件目录放开所有的权限，用csmartpro在计算节点进行运算；#area.pl提取出M的部分
			F3.write("%s/run2QCv2.sh %s %s %s %s %s %s\n"%(pipeline,sample,FC,R1,R2,primer_prefix,primer_prefix_org))
		print('\n'.join([csmart_analysis_path,FC,report_xls_name,report_zip_name,R1,R2,primer_prefix,primer_prefix_org]))
if ith == 0:
	print('sample_cfg_analysis文件夹中不存在该样本信息!')
elif ith >=2:
	print('Warning: 该样本可能存在加测!')


