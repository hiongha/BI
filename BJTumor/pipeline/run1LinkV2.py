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
from collections import OrderedDict

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
datapath1 = "/share/Onc_Pro/cSMART/sample_cfg_tmp/sample.cfg.*"
datapath2 = "/share/Onc_Pro/cSMART"
datapath3 = "/share/Onc_Pro/cSMART/sample_cfg_tmp/analysis_done_sample/sample.cfg.*"
primer_CRC = ['CRC_TM21_V1.3']
primer_LC = ['LC9_V4.10','LC9_V4.12','LC6_V5.4']
primer_SCLC = ['SCLC_V1.3']
primer_TM21 = ['TM21_V1.3']
primer_sfx = ['CRC','LC','SCLC','TM21']
org = ['crc_2.1','lung','SCLC','TM21']
primer_sfx2org = dict(zip(primer_sfx,org))


sample_cfgs = subprocess.Popen('grep %s %s'%(sample,datapath1),stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE,shell=True,universal_newlines=True)
processedFile = OrderedDict()

for each in sample_cfgs.stdout:
	percfgPath = each.strip().split(':')[0]
	if percfgPath in processedFile:
		continue
	else:
		processedFile[percfgPath] = each

ith = 0
for cfgPath in processedFile.keys():
	ith +=1
	each = processedFile[cfgPath]
	cfg = each.strip().split('\t')
	FC_short = cfg[0].split(':')[1]
	primer = cfg[2]
	print(primer)
	for sfx in primer_sfx:
		prm_tmp = 'primer_'+sfx
		if primer in locals()[prm_tmp]:
			primer_prefix = sfx
			primer_prefix_org = primer_sfx2org[primer_prefix]
	print primer_prefix,primer_prefix_org	
	sample_type = cfg[5]

	if re.search('ADD',cfgPath) == None:
		all_dir = datapath2+'/'+'*'+FC_short+'_'+primer_prefix+'/'+'*/sample.cfg.*'
	else:
		ftmp = open(cfgPath,'r')
		FC_short = ftmp.readline().strip().split('\t')[0]
		all_dir = datapath2+'/'+'*'+FC_short+'_'+primer_prefix+'_ADD'+'/'+'sample.cfg.*'
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
		#shell脚本中其他变量设置. 已有的变量:root_dir/sample/pipeline/FC/R1/R2/primer_prefix/primer_prefix_org
		fc_short = re.search(r'(.*)_(.*)_(.*)_(.*)_(.*)',FC).group(4)
		cs_qc = re.sub(r'(.*)_(.*)_(.*)_(.*)_(.*)',r'\1_\2_\3_\4_\5',FC)+'.xls'
		brief_file = sample+".brief.xls"
		must_given_xls = sample+'.must.given.xls'
		fq_R1_gz = FC+'/' + R1   #181123_TPNB500289_0074_AHGTK2BGX9_14/Z18L06604_L1_H511.R1.clean.fastq.gz
		fq_R1 = re.search(r'(.*.clean.fastq).gz',R1).group(1)   #Z18L06604_L1_H511.R1.clean.fastq
		fq_R2_gz = FC + '/' +R2   #181123_TPNB500289_0074_AHGTK2BGX9_14/Z18L06604_L1_H511.R2.clean.fastq.gz
		fq_R2 = re.search(r'(.*.clean.fastq).gz',R2).group(1)   #Z18L06604_L1_H511.R1.clean.fastq
		report_dir = re.sub(r'(.*_.*_.*_.*)_(.*)',r'\1_cSMART_%s'%(primer_prefix),FC)
		IntervalFile = 'egfr.'+primer_prefix_org+'_location.txt_'+primer

		with open(root_dir+'/mail_tmp.sh','w') as F1:
			F1.write( "/share/public/software/Python-2.7.13/bin/python %s/mail_bz_check.py -s %s -f %s,%s -w %s \n" %(pipeline,sample,report_xls_name,report_zip_name,os.getcwd()))
		with open(root_dir+'/check_jobs.sh','w') as F2:
			F2.write("source /home/xinghe/.bashrc;/share/public/software/Python-2.7.13/bin/python %s/check.py -u %s,csmartpro -s runbz1.sh\n"%(pipeline,os.getlogin()))
		with open(root_dir + '/run_draw_format%s.sh'%(ith),'w') as F4:
			F4.write('''
sh {pipeline}/base_quality.sh
python {pipeline}/qual.py #输出sample_quality.txt/error.txt 
/share/public/software/R-3.3.3/bin/Rscript {pipeline}/gc.R base.txt {sample}.gc.png
/share/public/software/R-3.3.3/bin/Rscript {pipeline}/qual.R sample_quality.txt {sample}.quality.png
#/share/public/software/R-3.3.3/bin/Rscript {pipeline}/error.R error.txt {sample}_error.png
#/share/public/software/R-3.3.3/bin/Rscript {pipeline}/insert.R insert_bp_freq.txt {sample}_insert.png
/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 {pipeline}/IntegrateQC.py --bz-qc result.txt --cs-qc {cs_qc} --brief-file {brief_file} --type {primer_prefix}
python {pipeline}/mutations_csv.py --brief-file {brief_file} --must-given {must_given_xls} --final-report {cs_qc} --type {primer_prefix_org} --out-put {sample}.mutations.csv >mut.log\n'''.format(**locals()))

		with open(root_dir+'/runbz%s.sh'%(ith),'w') as F3:
			F3.write('''samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '{print $1"\\t"$3"\\t"$4"\\t"$6"\\t"$10"\\t"$(NF-2)}'| perl %s/area.pl > reformed_bam_pos.txt\n\n'''%(FC,sample,pipeline))	

			F3.write('''#链接数据
echo -e "Sample:{sample}\\nType:plasma" >> result.txt
ln -fs {FC}/*{fc_short}*.xls .
ln -fs {FC}/{sample}.brief.xls .
ln -fs {FC}/{sample}.must.given.xls .
gunzip -c {root_dir}/{fq_R1_gz} > {root_dir}/{fq_R1} &
gunzip -c {root_dir}/{fq_R2_gz} > {root_dir}/{fq_R2}
#计算碱基分布
mkdir -p R1_base; cd R1_base; perl {pipeline}/ATGCN.pl {root_dir}/{fq_R1} & 
cd {root_dir};mkdir -p R2_base; cd R2_base; perl {pipeline}/ATGCN.pl {root_dir}/{fq_R2} &
cd {root_dir}
#计算总reads数及Q2_Q30碱基数:
perl {pipeline}/Q20_Q30_33.pl {root_dir}/{fq_R1_gz} > {root_dir}/Q20_Q30.R1
perl {pipeline}/Q20_Q30_33.pl {root_dir}/{fq_R2_gz} > {root_dir}/Q20_Q30.R2
wait
total_seq_num=$(wc -l < {root_dir}/Q20_Q30.R1);total_seq_num=$(bc <<< "scale=4;2*($total_seq_num-2)")
echo -e "Rawreads:$total_seq_num" >> result.txt
seq_strategy=`tail -1 Q20_Q30.R1 |cut -f1`;
Q20_R1=`tail -1 Q20_Q30.R1 |cut -f2`;Q20_R2=`tail -1 Q20_Q30.R2 |cut -f2`
Q30_R1=`tail -1 Q20_Q30.R1 |cut -f3`;Q30_R2=`tail -1 Q20_Q30.R2 |cut -f3`
let Q20=Q20_R1+Q20_R2;let Q30=Q30_R1+Q30_R2; echo -e "Q20碱基数\\t$Q20" ;echo -e "Q30碱基数\\t$Q30" 
Q20_pct=$(bc <<< "scale=5;$Q20/($total_seq_num*$seq_strategy)");echo -e "Q20:$Q20_pct" >> result.txt
Q30_pct=$(bc <<< "scale=5;$Q30/($total_seq_num*$seq_strategy)");echo -e "Q30:$Q30_pct" >> result.txt\n\n'''.format(**locals()))

			#1.插入片段大小按照.reformed.fastq 来统计,插入片段大小是最大数量的片段长度；
			F3.write('''#计算InsertSize
perl -e 'my $n=1;my $count=0;while(<>){chomp;my @a=split;if($n==1){$count=$a[4];$count=~s/XF:i://g;$n++;}elsif($n==2){my $tmp1=length($a[0])."\\n";print ($tmp1 x $count)."\\n";$n++;}elsif($n==3){$n++;}elsif($n==4){$n=1;}}' %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq | sort -k1,1n | uniq -c > reformed_length_distribution.txt\n'''%(FC,sample))
			F3.write('''insert_size=$(sort -k1,1nr reformed_length_distribution.txt|sed -r 's/\s+/\\t/g'|head -1|cut -f3)
echo -e "Insertsize:$insert_size"  >> result.txt #取得插入片段数量最大的长度作为结果\n\n''')

			F3.write('''export totalInsertFragNum=`less reformed_length_distribution.txt|awk '{sum=sum+$1}END{print sum}'`\n''')
			#echo -e "总共有多少插入片段\t$totalInsertFragNum" >> result.txt
			F3.write('''less reformed_length_distribution.txt|perl -e 'while(<>){chomp;my @a=split;$bp=$a[1];$freq=$a[0]/$ENV{'totalInsertFragNum'};print $bp."\\t"."$freq"."\\n";}' > %s/insert_bp_freq.txt\n\n'''%(root_dir))

			F3.write('''matched_to_human_read_num=$(samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '{print $1"\\t"$(NF-2)}' | sort - | uniq - | sed 's/XF:i://g' | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[1];}print $sum."\\n";');echo -e "Mapped_reads:$matched_to_human_read_num" >> result.txt\n\n'''%(FC,sample))

			F3.write('''total_matched_bp_to_human=$(samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '{print $1"\\t"length($10)"\\t"$(NF-2)}' | sed 's/XF:i://g' | sort - | uniq - | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[1]*$a[2];}print $sum."\\n";'); #echo -e "比对到人基因组上的总碱基数(包含没有match到基因组的位点):$total_matched_bp_to_human"  >> result.txt\n\n'''%(FC,sample))

			F3.write('''echo "total_matched_bp_to_human:$total_matched_bp_to_human"
samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam| awk '{print $6"\\t"$(NF-2)}'| sed 's/XF:i://g' | sort - > reform_bam_ciga_XF.txt\n'''%(FC,sample))

			F3.write('''matched_bp_to_human=$(perl -e 'while(<>){chomp;my @a=split;$a[0]=~s/(\d+)D//g;$a[0]=~s/(\d+)S//g;$a[0]=~s/(\d+)I//g;print $a[0]."\\t".$a[1]."\\n";}' reform_bam_ciga_XF.txt | perl -e 'my $sumcheck=0;while(<>){chomp;my @a=split;my @b=split(/M/,$a[0]);my $sum=0;for(my $m=0;$m<@b;$m++){$sum=$sum+$b[$m];}$sumcheck=$sumcheck+$sum*$a[1];}print $sumcheck."\\n";')
echo "Mapped_bases:$matched_bp_to_human"  >> result.txt #(只包含匹配结果为M的位点)
matched_bp_to_human_pct=$(bc <<< "scale=5;$matched_bp_to_human/$total_matched_bp_to_human")
#echo "Mapped_bases_ratio:$matched_bp_to_human_pct"  >> result.txt
python %s/Mapped_bases_ratio.py %s %s  #上一句echo替换成这一句.\n\n'''%(pipeline,FC,sample))

			F3.write('''awk '{print $2"\\t"$3"\\t"$4}' reformed_bam_pos.txt | perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=$a[0]."\\t".$a[1];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[2];}else{$hash{$key}=$a[2];}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' > reformed_bam_pos.txt_count\n\n''')
			F3.write('''#5.比对到预定区域的碱基数 （bp）——引物附近区域
perl {pipeline}/get_location.pl {pipeline}/config/{IntervalFile} {root_dir}/reformed_bam_pos.txt > reformed_bam_pos.txt_getlocation.txt\n\n'''.format(**locals()))
			F3.write('''matched_to_reserved_primer=$(awk '{sum=sum+$4}END{print sum}' reformed_bam_pos.txt_getlocation.txt )
echo "Mapped_bases_to_target:$matched_to_reserved_primer"  >> result.txt
matched_to_reserved_primer_pct=$(bc <<< "scale=5;$matched_to_reserved_primer/$matched_bp_to_human")
echo "Mapped_bases_ratio_to_target:$matched_to_reserved_primer_pct"  >> result.txt\n\n''')

			F3.write('''samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '{print $1"\\t"$3"\\t"$4"\\t"$6"\\t"$10"\\t"$(NF-2)}' | sort -k1,1 | perl {pipeline}/area_all.pl > %s/reformed_bam_pos_all.txt\n\n'''%(FC,sample,root_dir))

			F3.write('''perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=join("\\t",@a[0..2]);if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[3];}else{$hash{$key}=$a[3];}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' reformed_bam_pos_all.txt > reformed_bam_pos_all_count.txt
perl -e 'my %hash;while(<>){chomp;my @a=split;for(my $n=$a[1];$n<=$a[2];$n++){my $key=$a[0]."\\t".$n;if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[3];}else{$hash{$key}=$a[3];}}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' reformed_bam_pos_all_count.txt > reformed_bam_pos_all_count.txt_count1\n\n''')

			F3.write('''samtools view %s/%s.noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.sorted.bam | awk '{print $1"\\t"$3"\\t"$4"\\t"$6"\\t"$10"\\t"$(NF-2)}' | sort -k1,1 | perl %s/area_unique.pl > unique_bam_pos.txt\n\n'''%(FC,sample,pipeline))

			F3.write('''perl {pipeline}/get_location.pl {pipeline}/config/{IntervalFile} unique_bam_pos.txt > unique_bam_pos.txt_getlocation.txt\n'''.format(**locals()))


			F3.write('''matched_to_reserved_location_bp=`wc -l < unique_bam_pos.txt_getlocation.txt`
echo "Reads_uniquely_Mapped_bases_to_target:$matched_to_reserved_location_bp"  >> result.txt
matched_to_reserved_location_bp=`wc -l < unique_bam_pos.txt_getlocation.txt`
matched_to_reserved_location_bp_total=`wc -l < unique_bam_pos.txt`
matched_to_reserved_location_bp_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp/$matched_to_reserved_location_bp_total")
echo "Reads_uniquely_Mapped_bases_ratio_to_target:$matched_to_reserved_location_bp_pct"  >> result.txt\n\n''')

			F3.write('''awk '{print $2"\\t"$3}' unique_bam_pos.txt_getlocation.txt | sort -k1,1 -k2,2n | uniq -c > unique_bam_pos.txt_getlocation.txt_count
matched_to_reserved_location_bp_times5=$(awk '$1>=5{print}' unique_bam_pos.txt_getlocation.txt_count|wc -l)
echo "Bases_of_target_covered_at_least_5x:$matched_to_reserved_location_bp_times5"  >> result.txt
matched_to_reserved_location_bp_times10=$(awk '$1>=10{print}' unique_bam_pos.txt_getlocation.txt_count|wc -l)
matched_to_reserved_location_bp_all=`wc -l < unique_bam_pos.txt_getlocation.txt_count`
matched_to_reserved_location_bp_times5_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp_times5/$matched_to_reserved_location_bp_all")
matched_to_reserved_location_bp_times10_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp_times10/$matched_to_reserved_location_bp_all")
echo "Bases_ratio_of_target_covered_at_least_5x:$matched_to_reserved_location_bp_times5_pct"  >> result.txt
echo "Bases_of_target_covered_at_least_10x:$matched_to_reserved_location_bp_times10"  >> result.txt
echo "Bases_ratio_of_target_covered_at_least_10x:$matched_to_reserved_location_bp_times10_pct"  >> result.txt\n\n''')

#7.统计reform.unique_cons.fastq的每一个片段的长度以及其长度重复数；
			F3.write('''perl -e 'my $n=1;while(<>){chomp;my @a=split;if($n==1 or $n==3){$n++;}elsif($n==2){print length($a[0]);print "\\n\\n";$n++;}elsif($n==4){$n=1;}}' %s/%s.noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq | sort -k1,1n | uniq -c > unique_cons_length_distribution.txt\n\n'''%(FC,sample))
			F3.write('''awk '{print ($1*$2)}' unique_cons_length_distribution.txt | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];}print $sum."\\n\\n";' \n\n''')

			F3.write('''perl -e 'my %hash;my ($file1,$file2)=@ARGV;open IN1,"<$file1" or die "Can not open IN1";while(<IN1>){chomp;my @a=split;$hash{$a[4]}=$a[3];}close IN1;open IN2,"<$file2" or die "Can not open IN2";while(<IN2>){chomp;my @a=split;if(exists $hash{$a[0]}){print $hash{$a[0]}."\\t".$_."\\n";}else{print "-\\t".$_."\\n";}}close IN2;' ''' +  "%s/%s"%(FC,sample) + '''.noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.detail.txt unique_bam_pos.txt_getlocation.txt > unique_bam_pos.txt_getlocation.txt_addnum\n\n''')

			F3.write('''awk '{print $2"\\t"$3"\\t"$4}' reformed_bam_pos.txt_getlocation.txt | perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=$a[0]."\\t".$a[1];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[2];}else{$hash{$key}=$a[2];}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' > reformed_bam_pos.txt_getlocation.txt_count\n\n''')

			F3.write('''perl -e 'my $sum=0;while(<>){chomp;my @a=split(/:|-/,$_);print join("\\t",@a);print "\\t";print ($a[2]-$a[1]+1);print "\\n";$sum=$sum+($a[2]-$a[1]+1);}print $sum."\\n";' %s/config/%s\n\n'''%(pipeline,IntervalFile))

			F3.write('''awk '{print $2"\\t"$3}' reformed_bam_pos.txt | sort - | uniq - | wc -l
grep 'XP:Z:' %s/%s.noMS.pairPrimer.NSC.join.reformed.fastq | awk '{print $3"\\t"$5}' | sed 's/XF:i://g;s/XP:Z://g' | perl -e 'my %hash;while(<>){chomp;my @a=split;if(exists $hash{$a[0]}){$hash{$a[0]}=$hash{$a[0]}+$a[1];}else{$hash{$a[0]}=$a[1];}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' > reform_fastq_primer_count\n\n''')

			F3.write('''perl -e 'my %hash;my ($file1,$file2)=@ARGV;open IN1,"<$file1" or die "Can not open IN1";while(<IN1>){chomp;my @a=split;$hash{$a[4]}=$a[3];}close IN1;open IN2,"<$file2" or die "Can not open IN2";while(<IN2>){chomp;my @a=split;if(exists $hash{$a[0]}){print $hash{$a[0]}."\\t".$_."\\n";}else{print "-\\t".$_."\\n";}}close IN2;' ''' + "%s/%s"%(FC,sample)+'''.noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.detail.txt unique_bam_pos.txt > unique_bam_pos.txt_addnum\n\n''')

			F3.write('''perl -e 'my %hash;my $sum=0;while(<>){chomp;my @a=split;my $key=$a[2]."\\t".$a[3];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[0];}else{$hash{$key}=$a[0];}}foreach(keys %hash){print $_."\\t".$hash{$_}."\\n";}' unique_bam_pos.txt_getlocation.txt_addnum > unique_bam_pos.txt_getlocation.txt_addnum_count\n''')

			F3.write('''awk '{print $3"\\t"$4}' unique_bam_pos.txt_getlocation.txt_addnum | sort - | uniq -c > unique_bam_pos.txt_getlocation.txt_addnum_count1\n\n
#8.每个碱基的平均不重复覆盖数 
awk '{sum=sum+$1}END{print sum}' unique_bam_pos.txt_getlocation.txt_addnum_count1
wc -l unique_bam_pos.txt_getlocation.txt_addnum_count1
total_bp=$(perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];}print $sum."\\n";' unique_bam_pos.txt_getlocation.txt_count)
total_reads=$(wc -l < unique_bam_pos.txt_getlocation.txt_count)
non_average_cov_of_each_base=$(bc <<< "scale=5;$total_bp/$total_reads")
#echo "Average_depth_on_bases_uniquely:$non_average_cov_of_each_base"  >> result.txt #放到了最后一行. 
''')
			F3.write('''#9.每个碱基的平均覆盖数
average_coverage_of_each_base=$(awk '{print $NF}' reformed_bam_pos.txt_getlocation.txt_count | perl -e 'my $sum=0;my $n=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];$n++;}print $sum."\\t".$n."\\t";print ($sum/$n);print "\\n";' |cut -f3)
echo "Average_depth_on_bases:$average_coverage_of_each_base"  >> result.txt
echo "Average_depth_on_bases_uniquely:$non_average_cov_of_each_base"  >> result.txt \n''')

			F3.write('''
sh run_draw_format.sh
mkdir -p {root_dir}/{sample}_report/{report_dir}
mv *_cSMART_{primer_prefix}_bz.xls {root_dir}/{sample}_report
mv {root_dir}/{sample}.mutations.csv {root_dir}/*.png {sample}_report/{report_dir}
cd {root_dir}/{sample}_report
zip {report_dir}.zip {report_dir}/*
cd {root_dir}
'''.format(**locals()))

		print('\n'.join([csmart_analysis_path,FC,report_xls_name,report_zip_name,R1,R2,primer_prefix,primer_prefix_org]))
if ith == 0:
	print('sample_cfg_analysis文件夹中不存在该样本信息!')
elif ith >=2:
	print('Warning: 该样本可能存在加测!')
