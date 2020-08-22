#!/bin/bash
#author xingh xingh3223@berryoncology.com

#处理base结果文件:
readlength=$(awk '{print NF-1}' R2_base/base.txt|head -1)
awk -v readlength="$readlength" 'BEGIN{OFS="\t"}{if (NR==1) {for(i=2;i<=NF;i++){$i=$i+readlength}}print $0}' R2_base/base.txt |cut -f1 --complement > R2_base/base.txt.tmp
paste R1_base/base.txt R2_base/base.txt.tmp > base_paste.txt
sum=$(awk 'BEGIN{OFS="\t"}{if (NR>=2){for(i=2;i<=NF;i++){sum[i]+=$i}}}END{print sum[2]}' base_paste.txt)
awk -v sum="$sum" 'BEGIN{OFS="\t"}{if (NR==1){print $0}if(NR>=2){for(i=2;i<=NF;i++){$i=$i/sum*100}print $0}}' base_paste.txt > base.txt


#处理碱基质量结果的文件:
readlength=$(awk '{print NF-1}' R2_base/quality.txt |head -1)
awk -v readlength="$readlength" 'BEGIN{OFS="\t"}{if (NR==1) {for(i=2;i<=NF;i++){$i=$i+readlength}}print $0}' R2_base/quality.txt |cut -f1 --complement > R2_base/quality.txt.tmp

paste R1_base/quality.txt R2_base/quality.txt.tmp > quality.txt
