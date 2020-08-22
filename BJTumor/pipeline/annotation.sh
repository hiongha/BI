#/bin/bash
#suxiaoxing <suxiaoxing001@berryoncology.com>

perl /share/public/database/Annovar/2016Feb01/convert2annovar.pl --includeinfo --format vcf4 --allsample --outfile list1.avinput list1.vcf

perl /share/public/database/Annovar/2016Feb01/table_annovar.pl --buildver hg19 --thread 4  --remove --otherinfo --protocol refGene,wgEncodeGencodeBasicV19,mitimpact24,cytoBand,avsnp147,clinvar_20170130,cosmic70,icgc21,1000g2015aug_all,1000g2015aug_eas,exac03,esp6500siv2_all,popfreq_max_20150413,cadd,dann,gerp++gt2,dbnsfp33a,genomicSuperDups,wgRna,targetScanS,genomicSuperDups,dgvMerged,gwasCatalog,wgEncodeBroadHmmGm12878HMM,tfbsConsSites,phastConsElements46way,phastConsElements100way,Berry_400k,gnomad_exome,intervar_20170202 -operation g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r,r,r,r,r,r,f,f,f -nastring . list1.avinput /share/public/database/Annovar/2016Feb01/humandb --outfile list1.annovar

perl /share/public/database/Annovar/2016Feb01/annotate_variation.pl --geneanno --buildver hg19 --hgvs --thread 4 list1.avinput /share/public/database/Annovar/2016Feb01/humandb --outfile list1.hgvs
