#!/bin/bash
#Created by suxiaoxing <suxiaoxing001@berryoncology.com>;chenxiaoyan <chenxiaoyan712@berryoncology.com>
#Updated by xinghe <xingh3223@berryoncology.com>

#变量设置
pipeline=/share/work1/xinghe/proc/BZQC/pipeline
root_dir=`pwd`
#链接数据:

sample=$1
FC=$2 #181123_TPNB500289_0074_AHGTK2BGX9_14
R1=$3 #Z18L06604_L1_H511.R1.clean.fastq.gz
R2=$4 #Z18L06604_L1_H511.R2.clean.fastq.gz
fc_short=$(sed -r 's/(.*)_(.*)_(.*)_(.*)_(.*)/\4/' <<< $FC)
ln -fs $FC/*$fc_short*.xls .
ln -fs $FC/$sample".brief.xls" .
ln -fs $FC/$sample".must.given.xls" .
cs_qc=$(ls -d $FC/*$fc_short*.xls|xargs -l basename)
brief_file=$(ls -d $sample".brief.xls")
must_given_xls=$(ls -d $sample".must.given.xls")

fq_R1_gz=$FC/$R1   #181123_TPNB500289_0074_AHGTK2BGX9_14/Z18L06604_L1_H511.R1.clean.fastq.gz
fq_R1=$(sed -r 's/(.*.clean.fastq).gz/\1/' <<< $R1)   #Z18L06604_L1_H511.R1.clean.fastq
fq_R2_gz=$FC/$R2   #181123_TPNB500289_0074_AHGTK2BGX9_14/Z18L06604_L1_H511.R2.clean.fastq.gz
fq_R2=$(sed -r 's/(.*.clean.fastq).gz/\1/' <<< $R2)
echo -e "Sample:$1\nType:plasma" >> result.txt

#生成run_draw_format.sh 但是在整个脚本的最后才运行
echo "#draw pictures
sample=$sample
pipeline=$pipeline
FC=$FC
cs_qc=$cs_qc
brief_file=$brief_file
must_given_xls=$must_given_xls
sh $pipeline/base_quality.sh
python $pipeline/qual.py #输出sample_quality.txt/error.txt 
/share/public/software/R-3.3.3/bin/Rscript $pipeline/gc.R base.txt $sample".gc.png"
/share/public/software/R-3.3.3/bin/Rscript $pipeline/qual.R sample_quality.txt $sample".quality.png"
#/share/public/software/R-3.3.3/bin/Rscript $pipeline/error.R error.txt $sample"_error.png"
#/share/public/software/R-3.3.3/bin/Rscript $pipeline/insert.R insert_bp_freq.txt $sample"_insert.png"
" > run_draw_format.sh 
echo "/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 $pipeline/IntegrateQC.py --bz-qc result.txt --cs-qc $cs_qc --brief-file $brief_file " >> run_draw_format.sh 
echo "python $pipeline/mutations_csv.py --brief-file $brief_file --must-given $must_given_xls --final-report $cs_qc --out-put $sample"".mutations.csv >mut.log" >> run_draw_format.sh

#解压数据:
gunzip -c $root_dir/$fq_R1_gz > $root_dir/$fq_R1 &
gunzip -c $root_dir/$fq_R2_gz > $root_dir/$fq_R2
#计算碱基分布

mkdir -p R1_base; cd R1_base; perl $pipeline/ATGCN.pl $root_dir/$fq_R1 & 
cd $root_dir;mkdir -p R2_base; cd R2_base; perl $pipeline/ATGCN.pl $root_dir/$fq_R2 &
cd $root_dir

#计算Q2_Q30碱基数:
perl $pipeline/Q20_Q30_33.pl $root_dir/$fq_R1_gz > $root_dir/Q20_Q30.R1
perl $pipeline/Q20_Q30_33.pl $root_dir/$fq_R2_gz > $root_dir/Q20_Q30.R2
wait
echo "start to calculate Q20 and Q30"
total_seq_num=$(wc -l < $root_dir/Q20_Q30.R1);total_seq_num=$(bc <<< "scale=4;2*($total_seq_num-2)")
#echo -e "总测序片段数\t$total_seq_num" >> result.txt
echo -e "Rawreads:$total_seq_num" >> result.txt
seq_strategy=`tail -1 Q20_Q30.R1 |cut -f1`;
Q20_R1=`tail -1 Q20_Q30.R1 |cut -f2`;Q20_R2=`tail -1 Q20_Q30.R2 |cut -f2`
Q30_R1=`tail -1 Q20_Q30.R1 |cut -f3`;Q30_R2=`tail -1 Q20_Q30.R2 |cut -f3`
let Q20=Q20_R1+Q20_R2;let Q30=Q30_R1+Q30_R2; echo -e "Q20碱基数\t$Q20" ;echo -e "Q30碱基数\t$Q30" 
Q20_pct=$(bc <<< "scale=5;$Q20/($total_seq_num*$seq_strategy)");echo -e "Q20:$Q20_pct" >> result.txt
Q30_pct=$(bc <<< "scale=5;$Q30/($total_seq_num*$seq_strategy)");echo -e "Q30:$Q30_pct" >> result.txt

#1.插入片段大小按照.reformed.fastq 来统计,插入片段大小是最大数量的片段长度；
perl -e 'my $n=1;my $count=0;while(<>){chomp;my @a=split;if($n==1){$count=$a[4];$count=~s/XF:i://g;$n++;}elsif($n==2){my $tmp1=length($a[0])."\n";print ($tmp1 x $count)."\n";$n++;}elsif($n==3){$n++;}elsif($n==4){$n=1;}}' $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq" | sort -k1,1n | uniq -c > reformed_length_distribution.txt
insert_size=$(sort -k1,1nr reformed_length_distribution.txt|sed -r 's/\s+/\t/g'|head -1|cut -f3)
#echo -e "插入片段大小估计(bp)\t$insert_size"  >> result.txt #取得插入片段数量最大的长度作为结果
echo -e "Insertsize:$insert_size"  >> result.txt #取得插入片段数量最大的长度作为结果
export totalInsertFragNum=`less reformed_length_distribution.txt|awk '{sum=sum+$1}END{print sum}'`
#echo -e "总共有多少插入片段\t$totalInsertFragNum" >> result.txt
less reformed_length_distribution.txt|perl -e 'while(<>){chomp;my @a=split;$bp=$a[1];$freq=$a[0]/$ENV{'totalInsertFragNum'};print $bp."\t"."$freq"."\n";}' > $root_dir/insert_bp_freq.txt


matched_to_human_read_num=$(samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam" | awk '{print $1"\t"$(NF-2)}' | sort - | uniq - | sed 's/XF:i://g' | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[1];}print $sum."\n";');echo -e "Mapped_reads:$matched_to_human_read_num" >> result.txt

#3536289

#比对到人基因组上的总碱基数(包含没有match到基因组的位点): #(一条比对回基因组的read中的全部碱基可能是完全match,也可能是没有必读回去.两种统计方式：perl -e或者awk)
total_matched_bp_to_human=$(samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam" | awk '{print $1"\t"length($10)"\t"$(NF-2)}' | sed 's/XF:i://g' | sort - | uniq - | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[1]*$a[2];}print $sum."\n";'); #echo -e "比对到人基因组上的总碱基数(包含没有match到基因组的位点):$total_matched_bp_to_human"  >> result.txt
echo "total_matched_bp_to_human:$total_matched_bp_to_human"
#samtools view 181123_TPNB500289_0074_AHGTK2BGX9_14/Z18L06604.noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '{print $1"\t"length($10)"\t"$(NF-2)}' | sed 's/XF:i://g' | sort - | uniq - |awk '{sum=sum+$2*$3}END{print sum}'   ####???181116_TPNB500289_0071_AHGHMJBGX9_6/Z18L02663不对吧.
#600927501 总碱基数，包含未比对的
samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam"| awk '{print $6"\t"$(NF-2)}'| sed 's/XF:i://g' | sort - > reform_bam_ciga_XF.txt
matched_bp_to_human=$(perl -e 'while(<>){chomp;my @a=split;$a[0]=~s/(\d+)D//g;$a[0]=~s/(\d+)S//g;$a[0]=~s/(\d+)I//g;print $a[0]."\t".$a[1]."\n";}' reform_bam_ciga_XF.txt | perl -e 'my $sumcheck=0;while(<>){chomp;my @a=split;my @b=split(/M/,$a[0]);my $sum=0;for(my $m=0;$m<@b;$m++){$sum=$sum+$b[$m];}$sumcheck=$sumcheck+$sum*$a[1];}print $sumcheck."\n";')
#echo "比对到基因组的碱基总数:$matched_bp_to_human"  >> result.txt #(只包含匹配结果为M的位点)
echo "Mapped_bases:$matched_bp_to_human"  >> result.txt #(只包含匹配结果为M的位点)
matched_bp_to_human_pct=$(bc <<< "scale=5;$matched_bp_to_human/$total_matched_bp_to_human")
#echo "Mapped_bases_ratio:$matched_bp_to_human_pct"  >> result.txt
python $pipeline/Mapped_bases_ratio.py $FC $sample  #上一句替换成这一句.

#暂时注释掉下面这行:
#samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam" | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$10"\t"$(NF-2)}' | sort -k1,1 | $pipeline/perl area.pl > reformed_bam_pos.txt & #该步骤分析特别耗内存和存储，建议文件目录放开所有的权限，用csmartpro在计算节点进行运算；#area.pl提取出M的部分

awk '{print $2"\t"$3"\t"$4}' reformed_bam_pos.txt | perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=$a[0]."\t".$a[1];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[2];}else{$hash{$key}=$a[2];}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' > reformed_bam_pos.txt_count

#5.比对到预定区域的碱基数 （bp）——引物附近区域
perl $pipeline/get_location.pl $pipeline/config/egfr.crc_location.txt $root_dir/reformed_bam_pos.txt > reformed_bam_pos.txt_getlocation.txt
matched_to_reserved_primer=$(awk '{sum=sum+$4}END{print sum}' reformed_bam_pos.txt_getlocation.txt )
echo "Mapped_bases_to_target:$matched_to_reserved_primer"  >> result.txt
matched_to_reserved_primer_pct=$(bc <<< "scale=5;$matched_to_reserved_primer/$matched_bp_to_human")
echo "Mapped_bases_ratio_to_target:$matched_to_reserved_primer_pct"  >> result.txt
#perl get_location.pl location_clean.txt reformed_bam_pos.txt > reformed_bam_pos.txt_getlocation.txt
#444640816

samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam" | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$10"\t"$(NF-2)}' | sort -k1,1 | perl $pipeline/area_all.pl > $root_dir/reformed_bam_pos_all.txt
#600642982

perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=join("\t",@a[0..2]);if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[3];}else{$hash{$key}=$a[3];}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' reformed_bam_pos_all.txt > reformed_bam_pos_all_count.txt

perl -e 'my %hash;while(<>){chomp;my @a=split;for(my $n=$a[1];$n<=$a[2];$n++){my $key=$a[0]."\t".$n;if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[3];}else{$hash{$key}=$a[3];}}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' reformed_bam_pos_all_count.txt > reformed_bam_pos_all_count.txt_count1


samtools view $FC/$sample".noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.sorted.bam" | awk '{print $1"\t"$3"\t"$4"\t"$6"\t"$10"\t"$(NF-2)}' | sort -k1,1 | perl $pipeline/area_unique.pl > unique_bam_pos.txt
#wc -l unique_bam_pos.txt
#34790459

#perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[3];}print $sum."\n";' reformed_bam_pos.txt

perl $pipeline/get_location.pl $pipeline/config/egfr.crc_location.txt unique_bam_pos.txt > unique_bam_pos.txt_getlocation.txt
matched_to_reserved_location_bp=`wc -l < unique_bam_pos.txt_getlocation.txt`
echo "Reads_uniquely_Mapped_bases_to_target:$matched_to_reserved_location_bp"  >> result.txt
#27882117比对到预定区域的碱基数（bp）（除去重复的测序片段——模板的碱基数）
#perl get_location.pl location_clean.txt unique_bam_pos.txt > unique_bam_pos.txt_getlocation.txt

#6.比对到预定区域的碱基比例（除去重复的测序片段）
matched_to_reserved_location_bp=`wc -l < unique_bam_pos.txt_getlocation.txt`
matched_to_reserved_location_bp_total=`wc -l < unique_bam_pos.txt`
matched_to_reserved_location_bp_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp/$matched_to_reserved_location_bp_total")
echo "Reads_uniquely_Mapped_bases_ratio_to_target:$matched_to_reserved_location_bp_pct"  >> result.txt


awk '{print $2"\t"$3}' unique_bam_pos.txt_getlocation.txt | sort -k1,1 -k2,2n | uniq -c > unique_bam_pos.txt_getlocation.txt_count
matched_to_reserved_location_bp_times5=$(awk '$1>=5{print}' unique_bam_pos.txt_getlocation.txt_count|wc -l)
echo "Bases_of_target_covered_at_least_5x:$matched_to_reserved_location_bp_times5"  >> result.txt
matched_to_reserved_location_bp_times10=$(awk '$1>=10{print}' unique_bam_pos.txt_getlocation.txt_count|wc -l)
matched_to_reserved_location_bp_all=`wc -l < unique_bam_pos.txt_getlocation.txt_count`
matched_to_reserved_location_bp_times5_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp_times5/$matched_to_reserved_location_bp_all")
matched_to_reserved_location_bp_times10_pct=$(bc <<< "scale=5;$matched_to_reserved_location_bp_times10/$matched_to_reserved_location_bp_all")
echo "Bases_ratio_of_target_covered_at_least_5x:$matched_to_reserved_location_bp_times5_pct"  >> result.txt
echo "Bases_of_target_covered_at_least_10x:$matched_to_reserved_location_bp_times10"  >> result.txt
echo "Bases_ratio_of_target_covered_at_least_10x:$matched_to_reserved_location_bp_times10_pct"  >> result.txt

#7.统计reform.unique_cons.fastq的每一个片段的长度以及其长度重复数；
perl -e 'my $n=1;while(<>){chomp;my @a=split;if($n==1 or $n==3){$n++;}elsif($n==2){print length($a[0]);print "\n";$n++;}elsif($n==4){$n=1;}}' $FC/$sample".noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq" | sort -k1,1n | uniq -c > unique_cons_length_distribution.txt

awk '{print ($1*$2)}' unique_cons_length_distribution.txt | perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];}print $sum."\n";'

perl -e 'my %hash;my ($file1,$file2)=@ARGV;open IN1,"<$file1" or die "Can not open IN1";while(<IN1>){chomp;my @a=split;$hash{$a[4]}=$a[3];}close IN1;open IN2,"<$file2" or die "Can not open IN2";while(<IN2>){chomp;my @a=split;if(exists $hash{$a[0]}){print $hash{$a[0]}."\t".$_."\n";}else{print "-\t".$_."\n";}}close IN2;' $FC/$sample".noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.detail.txt" unique_bam_pos.txt_getlocation.txt > unique_bam_pos.txt_getlocation.txt_addnum
#438843469

awk '{print $2"\t"$3"\t"$4}' reformed_bam_pos.txt_getlocation.txt | perl -e 'my %hash;while(<>){chomp;my @a=split;my $key=$a[0]."\t".$a[1];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[2];}else{$hash{$key}=$a[2];}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' > reformed_bam_pos.txt_getlocation.txt_count


perl -e 'my $sum=0;while(<>){chomp;my @a=split(/:|-/,$_);print join("\t",@a);print "\t";print ($a[2]-$a[1]+1);print "\n";$sum=$sum+($a[2]-$a[1]+1);}print $sum."\n";' $pipeline/config/egfr.crc_location.txt
#56297
#56297-20000

awk '{print $2"\t"$3}' reformed_bam_pos.txt | sort - | uniq - | wc -l
#175149

grep 'XP:Z:' $FC/$sample".noMS.pairPrimer.NSC.join.reformed.fastq" | awk '{print $3"\t"$5}' | sed 's/XF:i://g;s/XP:Z://g' | perl -e 'my %hash;while(<>){chomp;my @a=split;if(exists $hash{$a[0]}){$hash{$a[0]}=$hash{$a[0]}+$a[1];}else{$hash{$a[0]}=$a[1];}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' > reform_fastq_primer_count

perl -e 'my %hash;my ($file1,$file2)=@ARGV;open IN1,"<$file1" or die "Can not open IN1";while(<IN1>){chomp;my @a=split;$hash{$a[4]}=$a[3];}close IN1;open IN2,"<$file2" or die "Can not open IN2";while(<IN2>){chomp;my @a=split;if(exists $hash{$a[0]}){print $hash{$a[0]}."\t".$_."\n";}else{print "-\t".$_."\n";}}close IN2;' $FC/$sample".noMS.pairPrimer.NSC.join.reformed.unique_cons.fastq.detail.txt" unique_bam_pos.txt > unique_bam_pos.txt_addnum


perl -e 'my %hash;my $sum=0;while(<>){chomp;my @a=split;my $key=$a[2]."\t".$a[3];if(exists $hash{$key}){$hash{$key}=$hash{$key}+$a[0];}else{$hash{$key}=$a[0];}}foreach(keys %hash){print $_."\t".$hash{$_}."\n";}' unique_bam_pos.txt_getlocation.txt_addnum > unique_bam_pos.txt_getlocation.txt_addnum_count

awk '{print $3"\t"$4}' unique_bam_pos.txt_getlocation.txt_addnum | sort - | uniq -c > unique_bam_pos.txt_getlocation.txt_addnum_count1

#8.每个碱基的平均不重复覆盖数 
awk '{sum=sum+$1}END{print sum}' unique_bam_pos.txt_getlocation.txt_addnum_count1
wc -l unique_bam_pos.txt_getlocation.txt_addnum_count1
total_bp=$(perl -e 'my $sum=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];}print $sum."\n";' unique_bam_pos.txt_getlocation.txt_count)
total_reads=$(wc -l < unique_bam_pos.txt_getlocation.txt_count)
non_average_cov_of_each_base=$(bc <<< "scale=5;$total_bp/$total_reads")
#echo "Average_depth_on_bases_uniquely:$non_average_cov_of_each_base"  >> result.txt #放到了最后一行. 

#9.每个碱基的平均覆盖数
average_coverage_of_each_base=$(awk '{print $NF}' reformed_bam_pos.txt_getlocation.txt_count | perl -e 'my $sum=0;my $n=0;while(<>){chomp;my @a=split;$sum=$sum+$a[0];$n++;}print $sum."\t".$n."\t";print ($sum/$n);print "\n";' |cut -f3)
echo "Average_depth_on_bases:$average_coverage_of_each_base"  >> result.txt
echo "Average_depth_on_bases_uniquely:$non_average_cov_of_each_base"  >> result.txt 


sh run_draw_format.sh
report_dir=$(sed -r 's/(.*_.*_.*_.*)_(.*)/\1_cSMART_CRC/' <<< $FC)
mkdir -p $root_dir/$sample"_report"/$report_dir
mv *_cSMART_CRC_bz.xls $root_dir/$sample"_report"
mv $root_dir/$sample".mutations.csv" $root_dir/*.png $sample"_report"/$report_dir

cd $root_dir/$sample"_report"
zip $report_dir".zip" $report_dir/*
cd $root_dir
#tar -zcvf $sample'_report.tar.gz' $sample"_report"


