#/usr/bin/python
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>

from __future__ import division

f=open('quality.txt','r')
title=f.readline()
Q_sum=[0]*(len(title.strip().split())-1)
Read_sum=[0]*(len(title.strip().split())-1)
mean=[0]*(len(title.strip().split())-1)
error=[0]*(len(title.strip().split())-1)

for line in f:
	line_list=line.strip().split()
	qual=line_list[0]
	readnumPos=line_list[1:]
	readnumMultiplyQual=[int(qual)*int(each) for each in readnumPos]
	for i,j in enumerate(readnumMultiplyQual):
		Q_sum[i]+=j
	for i,j in enumerate(readnumPos):
		Read_sum[i]+=int(j)

for i,j in enumerate(Q_sum):
	mean[i]=Q_sum[i]/Read_sum[i]
	error[i]=10**(mean[i]/(-10))

#为绘制质量值分布输出输入文件:
f_out=open('sample_quality.txt','w') 
f_out.write(title)
f_out.write("mean"+"\t"+'\t'.join([str(each) for each in mean])+'\n')
f_out.close()

#为绘制错误率分布图输出输入文件:
f_out=open('error.txt','w') 
f_out.write(title)
f_out.write("error"+"\t"+'\t'.join([str(each) for each in error])+"\n")
f_out.close()


