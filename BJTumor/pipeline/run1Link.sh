#sample=$1
#Author: xinghe
#Contact: xingh3223@berryoncology.com

pipeline=/share/work1/xinghe/proc/BZQC/pipeline
root_dir=`pwd`
info=info
#==========整理待合并位点信息文件格式,使其标准化======================================
sed -r 's/\s+/\t/g' $info > info.txt
sample=`head -1 info.txt|cut -f1`
echo -e "待分析的样品为: $sample\n"
#=============寻找并链接数据到当前路径=============================================
datapath1="/share/Oncology/production/cSMART/cSMART/analysis/sample_cfg_tmp/sample.cfg.*"
datapath2="/share/Oncology/production/cSMART/cSMART/analysis"
datapath3="/share/Oncology/production/cSMART/cSMART/analysis/sample_cfg_tmp/analysis_done_sample/sample.cfg.*"
sample_cfg=$(grep -l "$sample" $datapath1)
if [[ -z $sample_cfg ]];then
	sample_cfg=$(grep -l "$sample" $datapath3)
#elif sed 's/ /\n/g' <<< $sample_cfg|wc -l 
fi
echo -e "样品所在sample.cfg文件: $sample_cfg\n"

ipathNum=0
echo $sample_cfg
for ipath in $sample_cfg
do
	((ipathNum++))
	flowcell=$(cut -f1 $ipath|head -1)
	primer=$(cut -f3 $ipath|head -1)
	echo -e "引物是: $primer" #CRC_TM21_V1.3,LC9_V4.12
	primer_CRC=['CRC_TM21_V1.3']
	primer_LC=['LC9_V4.10','LC9_V4.12']
#	if [[ -n "${primer_LC[$primer]}" ]];then
#		suffix='LC'
#	elif [[ -n "${primer_CRC[$primer]}" ]];then
#		suffix='CRC'
#	fi
#	echo -e "$suffix\n"
	sample_type=$(sed -r 's/(.*)_.*/\1/' <<< $primer)
	echo -e "sample类型是:$sample_type"
	echo -e "样品所在的第$ipathNum个flowcell: $flowcell"
	sample_path=$(grep -l "$sample" $datapath2/*$flowcell*/*/sample.cfg.*|xargs -l dirname)
	echo -e "样品所在的第$ipathNum个路径: $sample_path"
	ln -fs $sample_path .
	R1=$(ls $sample_path/*.R1.clean.fastq.gz|xargs -l basename)
	R2=$(ls $sample_path/*.R2.clean.fastq.gz|xargs -l basename)
	#echo "R2 $R2"
	FC=$(basename $sample_path)
	report_xls_name=$(sed -r 's/(.*_.*_.*_.*)_.*/\1_cSMART_LC_bz\.xls/' <<< $FC) #190121_NS500511_0183_AH3LHTAFXY_cSMART_LC_bz.xls
	report_zip_name=$(sed -r 's/(.*_.*_.*_.*)_.*/\1_cSMART_LC\.zip/' <<< $FC) #190121_NS500511_0183_AH3LHTAFXY_cSMART_LC.zip
	echo $report_xls_name
	echo $report_zip_name

	echo -e "样品所在的第$ipathNum个runid: $FC\n"
	echo "/share/public/software/Python-2.7.13/bin/python $pipeline/mail_bz_check.py -s $sample -f $report_xls_name,$report_zip_name -w `pwd` " > 'mail_tmp.sh'
	echo "source /home/xinghe/.bashrc;/share/public/software/Python-2.7.13/bin/python $pipeline/check.py -u `whoami`,csmartpro -s runbz1.sh" > "check_jobs.sh"
	echo "samtools view $FC/$sample"".noMS.pairPrimer.NSC.join.reformed.fastq.sorted.bam | awk '" '{print $1"\t"$3"\t"$4"\t"$6"\t"$10"\t"$(NF-2)}' "'| sort -k1,1 |perl $pipeline/area.pl > reformed_bam_pos.txt" > "runbz"$ipathNum".sh"  #该步骤分析特别耗内存和存储，建议文件目录放开所有的权限，用csmartpro在计算节点进行运算；#area.pl提取出M的部分
	echo "$pipeline/run2QC.sh $sample $FC $R1 $R2 " >> "runbz"$ipathNum".sh"

	if [[ $ipathNum -ge 2 ]];then
		echo "Warning: 数据有加测或重测,请确认使用哪个flowcell上的数据"
	fi	
done


