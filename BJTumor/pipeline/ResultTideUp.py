#!/usr/bin/python
#coding=utf-8
import xlrd as xl
import csv
import os

def get_content_from_xls(file,sheetname):
	rb = xl.open_workbook(file)
	rs = rb.sheet_by_name(sheetname)
	nrow = rs.nrows
	ncol = rs.ncols
	content = []
	for i in range(0,nrow):
		line = rs.row_values(i)
		line1 = map(lambda x:str(x),line)
		new_line = '\t'.join(line1)
		content.append(new_line)
	if sheetname == 'summary':
		new_content = []
		for i in content:
			if i.startswith('#Sample'):
				new_content.append(i)
		content = new_content	
	return content

def write_to_csv(outputfile,content):
	mutcsvF=open(outputfile+'temp','w')
	for each in content:
		mutcsvF.write(each+'\n')
	mutcsvF.close()
	in_txt = csv.reader(open(outputfile+'temp','r'),delimiter='\t')
	out_csv = csv.writer(open(outputfile,'w'))	
	out_csv.writerows(in_txt)
	os.system('rm %s'%(outputfile+'temp'))

if __name__ == '__main__':
	import argparse
	parser=argparse.ArgumentParser(prog="",description="生成fusion.csv/SnvIndel.csv/QC.csv")
	parser.add_argument("-s","--sample",help='样本名称')
	parser.add_argument("-x","--xls-name",help='如191011_TPNB500210_0317_AHK5KWAFXY_cSMART_CRC_bz.xls')
	args=parser.parse_args()
	print(args)
#	file = '191011_TPNB500210_0317_AHK5KWAFXY_cSMART_CRC_bz.xls'
#	sample = 'Z19L03656'
	file = args.xls_name
	sample = args.sample
	QC = get_content_from_xls(file=file,sheetname='QC')
	write_to_csv(sample+'.QC.csv',QC)
	SnvIndel = get_content_from_xls(file=file,sheetname='SnvIndel')
	print(SnvIndel)
	write_to_csv(sample+'.SnvIndel.csv',SnvIndel)
	fusion = get_content_from_xls(file=file,sheetname='summary')
	write_to_csv(sample+'.fusion.csv',fusion)
		

