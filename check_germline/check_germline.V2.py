#!/usr/bin/python
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>
import re
import os
import time
import glob
import subprocess
import sys
import smtplib,ssl
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.header import Header
from email.utils import formatdate
from email import encoders
from collections import OrderedDict
import paramiko


pipeline = '/share/work1/xinghe/proc/captools/check_germline'

def mail_1(dir,to_list,cc_list,subject,main_text,mail_user='oncology_bioinfo',mail_pass = 'herui2020.com',sender = 'oncology_bioinfo@berryoncology.com',formatting='plain'):
	mail_host = 'mail.berryoncology.com'
	receivers = to_list+cc_list
	message = MIMEMultipart()
	message['From'] = sender
	message['To'] = ';'.join(to_list)
	message['Cc'] = ';'.join(cc_list)
	message['Date'] = formatdate(localtime=True)
	message['Subject'] = Header(subject,'utf-8')
	message.attach(MIMEText(main_text,formatting,'utf-8'))

	try:
		smtpObj = smtplib.SMTP()
		smtpObj.connect(mail_host,25)
		isTls = False
		if isTls:
			smtpObj.starttls()
		smtpObj.login(mail_user,mail_pass)
		smtpObj.sendmail(sender,receivers,message.as_string())
		smtpObj.quit()
		print("邮件发送成功")
	except smtplib.SMTPException:
		print("Error:无法发送邮件")


#将已经发送germline审核邮件的flowcell写入done.txt文件中
def rewrite(FC,done_file):
	with open(done_file,'a') as F2:
		F2.write(FC+'\n')


#获取没有进行审核germline的分析路径：
def get_undone_fc(done_file,typ):
	done_fc = []
	with open(done_file) as F:
		for line in F:
			done_fc.append(line.strip())	
	dirs = glob.glob(root+'/*')
	dirs = map(os.path.basename,dirs)
	undone_fc = set(dirs)-set(done_fc)
	return undone_fc


def check_Somatic_interpretation(path):
	somatic_file = glob.glob(path+'/Somatic_interpretation/tumors_Clisig_reviewed.txt')	
	Tumor_check_list = []
	somatic_finish_flag = False	
	somatic_check_flag = False
	Tumor_check_line = ''
	Tumor_check_num = 0
	if somatic_file == []:
		somatic_finish_flag = False
	else:
		somatic_file = somatic_file[0]
		somatic_finish_flag = True
		with open(somatic_file,'r') as F:
			content = F.read()
			parts = content.split('--------------------')
		
		for each in parts:
#			if re.search('该基因涉及DCE2.0用药，请对该突变ACMG判读结果进行审核',each):
				Tumor_check_list.append(each)
		if len(Tumor_check_list) == 0:
			somatic_check_flag = False
			Tumor_check_line = ''
		elif len(Tumor_check_list) >=1:
			somatic_check_flag = True
			Tumor_check_line = '--------------------'.join(Tumor_check_list)
			Tumor_check_num = len(Tumor_check_list) - 1 
			print(Tumor_check_line)
		
	return somatic_finish_flag,somatic_check_flag,Tumor_check_line,Tumor_check_num


def QC_check(path):
	file = path + '/QC.xls.new'
	check_cols = [40,41,42,44,46]
	QC_check_dic = OrderedDict()
	if not os.path.exists(file):
		pass
	else:
		with open(file,'r') as F:
			F.readline()
			for line in F:
				arr = line.strip().split('\t')
				print(line)
				check,Pairs,Contaminated,Contamination_level,QC_other_check =  [arr[i] for i in check_cols]	
				print(check,Pairs,Contaminated,Contamination_level,QC_other_check)
				Nsamp,Psamp = Pairs.split(':')			
				QC_check_dic[Psamp] = [check,Pairs,Contaminated,Contamination_level,QC_other_check]		
	return QC_check_dic


def write_to_html(List):
	print(List)
	html = '<p>Hi all,</p><p></p><p>MSI results:</p><table style="font-size:13px" border="1" cellpadding="0" cellspacing="0" white-space:nowrap>'
	for each in List:
		html += '<tr>'
		for i in each:
			html +='<td> %s </td>'%(str(i))
		html += '</tr>'
	html+= '</table>'
	print('html=========')
	print(html)
	return html


def msi_check(path):
	file = path + '/msi/msi_combine.xls'
	msiList = [['sample','MSIscore','MSIscore(0.49)','check','Pairs','Contaminated','Contamination_level','QC_other_check']]
	QC_check_dic = QC_check(path)
	if not os.path.exists(file):
		pass
	else:
		with open(file,'r') as F:
			for line in F:
				if line.startswith('ProjectName') or re.search(r'Silico',line.strip()):
					continue
				else:
					arr = line.strip().split('\t')
					sample = arr[1]
					MSIscore = arr[2]
					if float(MSIscore) >= 0.49:
						flag = 'P'
					elif float(MSIscore) < 0.49:
						flag = 'N'
					QC = QC_check_dic.get(sample,['-','-','-','-','-'])
					msi_info = [sample,MSIscore,flag]
					msi_info.extend(QC)
					msiList.append(msi_info)
	html = write_to_html(msiList)	
	return html	
		

def multi_check(root,panel):
	for typ in ['germline','pair','finish']:
		type_undone = []
		type_done_file = pipeline+'/%s.done.txt'%(typ)
		undone_fc = get_undone_fc(done_file = type_done_file,typ = typ)
	
		if typ == 'germline':
			for dir in undone_fc:
				germline = root+'/%s/germline/normals_germline_reviewed.txt'%(dir)
				somatic_finish_flag,somatic_check_flag,Tumor_check_line,Tumor_check_num = check_Somatic_interpretation(root+'/'+dir)

				if (not os.path.exists(germline)) or (somatic_finish_flag==False):
					continue

				p = subprocess.Popen("grep 共计审核 %s"%(germline),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

				content_grep = ''.join([i.strip() for i in p.stdout])
				count = len(content_grep)
				if count == 0:
					pass
				else:
					shenhe_count = re.search(r'共计审核(\d+)个突变位点',content_grep).group(1)
					if int(shenhe_count) == 0 and not (somatic_finish_flag==True and somatic_check_flag==True):
						to_list = ['xingh3223@berryoncology.com']
						cc_list = []
					else:
						if re.search(r'86gene|654gene|457gene|31gene|WES_Plus',root):
							to_list = ['wun3623@berryoncology.com','zhangp4456@berryoncology.com']
							cc_list = ['baijian488@berryoncology.com','jiangdzh3403@berryoncology.com','yangrutao796@berryoncology.com','wangzx3872@berryoncology.com','wujq3870@berryoncology.com','xufl3252@berryoncology.com','wangwt4475@berryoncology.com','liull4492@berryoncology.com','fangmg3899@berryoncology.com','xingh3223@berryoncology.com']
					
					content = []
					flag = True
					with open(germline,'r') as F1:
						for each in F1:
							if flag == True:
								content.append(each)
					content = ''.join(content)
					print(content)
	
					pF = open('%s/%s/pair.txt'%(root,dir),'r')
					sample_count = str(len(pF.readlines()))
					

					subject = '[%s审]%s germline位点审核'%(panel.strip('gene'),dir)
					main_text = '\n\nHi ,\n\n\t%s批次germline分析结果,请审核。\n\t共%s对样本,%s个位点需要审核。\n\t数据路径：%s/%s\n\n\n%s\n\n ====================================\n\nSomatic审核:\n\n%s'%(dir,sample_count,str(int(shenhe_count)+Tumor_check_num),root,dir,content,Tumor_check_line)
					print('nimalegejide::::::::::::::::::::::::::::::::::::::::::::')
					print(str(shenhe_count))
					print(str(Tumor_check_num))
					if re.search(r'86gene|654gene|457gene|31gene|WES_Plus',root):
						mail_1(dir=dir,to_list = to_list,cc_list = cc_list,subject=subject,main_text=main_text,mail_user='oncology_mut_review',mail_pass = 'herui123.com',sender = 'oncology_mut_review@berryoncology.com')
					rewrite(FC=dir,done_file=type_done_file)
		if typ == 'cnv':
			for dir in undone_fc:
				CNV_file = root+'/%s/CNV.xls'%(dir)
				if not os.path.exists(CNV_file):
					continue
				subject = '[%sCNV]%s CNV审核'%(panel.strip('gene'),dir)
				main_text = '\n\nHi 付龙飞,\n\n\t%s批次CNV分析结果已完成,请审核。\n\t数据路径：%s/%s\n\n\n '%(dir,root,dir)
				cc_list = []
				mail_1(dir=dir,to_list = ['xingh3223@berryoncology.com'],cc_list = cc_list,subject=subject,main_text=main_text)
				rewrite(FC=dir,done_file=type_done_file)
		if typ == 'pair':
			for dir in undone_fc:
				pair_file = root+'/%s/getSamplecfg.log'%(dir)
				pair_txt = root+'/%s/pair.txt'%(dir)
				if os.path.exists(pair_txt):	
					pair_txt_content = ''.join([line for line in open(pair_txt,'r')])
				else:
					pair_txt_content = pair_txt
				if not os.path.exists(pair_file):
					continue
				p = subprocess.Popen("grep no %s"%(pair_file),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				sample_no_pair = [i.strip() for i in p.stdout]
				count_pair = len(sample_no_pair)
				if count_pair == 0:
					subject = '[%s配]%s 样本可以不match'%(panel.strip('gene'),dir)
					main_text = '\n\nHi 配对,\n\n\t%s批次不含有不配对样本。\n\npair.txt:\n%s\t数据路径：%s/%s\n\n\n '%(dir,pair_txt_content,root,dir)
					mail_1(dir=dir,to_list = ['xingh3223@berryoncology.com'],cc_list = [],subject=subject,main_text=main_text)
				else:
					subject = '[%s配]%s 样本可以不match'%(panel.strip('gene'),dir)
					main_text = '\n\nHi 菜鸟,\n\n\t%s批次可能含有不配对样本:\n%s。\n\npair.txt:\n%s数据路径：%s/%s\n\n\n '%(dir,'\n'.join(sample_no_pair),pair_txt_content,root,dir)
					mail_1(dir=dir,to_list = ['xingh3223@berryoncology.com'],cc_list = [],subject=subject,main_text=main_text)
				rewrite(FC=dir,done_file=type_done_file)

		if typ == 'finish':
			for dir in undone_fc:
				warn_file = root+'/%s/warnClosePos.txt'%(dir)
				final_result_file = root+'/%s/%s.xls'%(dir,dir)
				print(dir)
				dir2 = dir.replace('_sup','').replace('_LC','').replace('_CRC','').replace('_info','').replace('_auto','').replace('_info','').replace('_sup','')
				print('\n##########################################################')
				print("newDir:"+dir2)
				third_qsub = root+'/%s/%s.job.status.status.status'%(dir,dir2)
				print("statusFile:"+third_qsub)
				time_file = root +'/'+dir +'/' + dir2.replace('_sup','').replace('_auto','')+'.time.xls'


				print("timeFile:"+time_file)	
				if not (os.path.exists(third_qsub) and os.path.exists(time_file)):
					continue

				#flag1 status_done	
				p3 = subprocess.Popen("grep status %s"%(third_qsub),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				status_done = True
				for eachstatus in p3.stdout:
					if eachstatus.strip() != 'status done':
						status_done = False
						break
				if status_done == False:
					continue
				no_dce_file = glob.glob(root+'/'+dir+'/*.no_dce.xls')
				#200529_Novaseq_A00838_B_BH3MHGDSXY_csmart2.0_SPS_v20_31_LC.xls.no_dce.xls
				if re.search(r'LC_CRC',root) and no_dce_file == []:
					continue
				else:
					time.sleep(2)

				time.sleep(2)	
	
				#flag2 germline_flag2	
				germline_flag2 = False
				germline2 = root+'/%s/germline/normals_germline_reviewed.txt'%(dir)
				p2 = subprocess.Popen("grep 共计审核 %s"%(germline2),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				content_grep2 = ''.join([i.strip() for i in p2.stdout])
				if re.search(r'共计审核0个突变位点',content_grep2):
					germline_flag2 = True

				somatic_finish_flag2,somatic_check_flag2,Tumor_check_line,Tumor_check_num = check_Somatic_interpretation(root+'/'+dir)
				somatic_flag3 = False

				if somatic_finish_flag2 == True and somatic_check_flag2 == False:
					somatic_flag3 = True


				#flag5 total_close_sites	
				total_close_sites = []
				pos_list = []
				samp = '';chr='';start='';

				need_checked_list = []
				with open(warn_file,'r') as warnF:
					for line in warnF:
						if re.search(r'^PC|^NC|^Silico',line):
							continue
						elif line == '\n'  and len(need_checked_list) == 0 :
							continue
						elif line == '\n'  and need_checked_list[-1] == '\n' :
							continue
						else:
							need_checked_list.append(line.strip())
				if 1 == 1:					
					for line in need_checked_list:
						if line.strip() != '':
							array = line.strip().split('\t')
							total_close_sites.append(line.strip())
							print(array[0:3])
							samp,chr,start = array[0:3]
							this_gene = array[8].split(':')[0]
							zhushi_gene = NM_dic1.get(this_gene,'--') +'/'+ NM_dic2.get(this_gene,'--') +'/'+ NM_dic3.get(this_gene,'--') 
							pos_list.append(int(start))
							min_pos = str(min(pos_list)-1)
						elif line.strip() == '' :
							samcmd = '\n/share/public/software/samtools-1.3/samtools tview %s/%s/./aln/%s/*.cons.bam --reference /share/work3/capsmart/pipeline/capSMART/CAPcSMART/capSMART/cSMART170503/reference/hg19.fasta -p %s:%s|grep -i [1atgc*]|les\n\nMainTranscript:%s\n\n=================\n\n'%(root,dir,samp,chr,min_pos,zhushi_gene)
							total_close_sites.append(samcmd)
							pos_list = []
							samp = '';chr='';start='';

				time_sheet = dir2.replace('_sup','').replace('_auto','')+'.time.xls'
				if total_close_sites == []:
					mail_shell = '【请在合并为点后执行】:\ncd %s/%s;\n/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.test.py;\n/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.all.py;\n'%(root,dir)
				else:
					mail_shell = '【请在合并位点后执行】:\ncd %s/%s;\n/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.test.py -merge;\n/share/work2/lisuxing/suxx/software/Python-3.6.1/bin/python3 /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.all.py -merge;\n'%(root,dir)
	
				warn_text = '\n'.join(total_close_sites) + '\n' + mail_shell	
				subject = '[%s完]%s 【项目分析完成通知】'%(panel.strip('gene'),dir)
				main_text = '\n\nHi 分析完成,\n\n\t%s批次马上完成分析，待germline/CNV/需合并位点整理完毕后再向报告组反馈结果。\n\n\t数据路径：%s/%s\n\n【待合并的位点及命令】:\n\n%s\n '%(dir,root,dir,warn_text)
				mail_1(dir=dir,to_list = ['xingh3223@berryoncology.com'],cc_list = [],subject=subject,main_text=main_text)

				#flag6 msi_check
				if panel == '654gene' and re.search(r'_654T_',dir):
					msi_html = msi_check(root+'/'+dir)
					subject = '[%sMSI] %s '%(panel.strip('gene'),dir)
					msi_cc_list = ['jiangdzh3403@berryoncology.com','wangzx3872@berryoncology.com','wujq3870@berryoncology.com','wun3623@berryoncology.com','xingh3223@berryoncology.com' ]		
					mail_1(dir=dir,to_list = ['oncology_labreport@berryoncology.com'],cc_list = msi_cc_list,subject=subject,main_text=msi_html,formatting='html')

			
				rewrite(FC=dir,done_file=type_done_file)
				
				#flag7 hebing_flag2
				hebing_flag2 = False
				if total_close_sites == []:
					hebing_flag2 = True
			
				#flag8 pair_txt_list_flag
				pair_txt = root+'/%s/pair.txt'%(dir)
				pair_txt_list_flag = True
				if os.path.exists(pair_txt):	
					pair_txt_list = [line.strip() for line in open(pair_txt,'r') if not re.search(r'^Sili',line)]
				else:
					pair_txt_list = []
				if pair_txt_list == []:
					pair_txt_list_flag = False
				

				#judgement of different panel	
				print('status_done,hebing_flag2,germline_flag2,pair_txt_list_flag')
				print(status_done,hebing_flag2,germline_flag2,pair_txt_list_flag)
				
				if dir == '200529_Novaseq_A00838_B_BH3MHGDSXY_csmart2_0_SPL_v11_457_auto_sup':
					ssh = paramiko.SSHClient()
					ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
					ssh.connect(hostname='10.100.2.10',port=22,username='capsmart',password='berry2012')
					stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
					stdin, stdoutline, stderr = ssh.exec_command("export LANG=zh_CN.UTF-8;export LC_ALL=en_GB.UTF-8;cd %s/%s;/share/work1/xinghe/software/anaconda3/bin/python /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.all.py"%(root,dir))

				#不管是否有需要checkgermline的位点
				elif status_done == True and hebing_flag2 == True and re.search(r'_31gene',root) and pair_txt_list_flag == True and re.search(r'20080[67]',dir) and germline_flag2 == True and somatic_flag3 == True:
				#germline等都需要check
					ssh = paramiko.SSHClient()
					ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
					ssh.connect(hostname='10.100.2.10',port=22,username='capsmart',password='berry2012')
					stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
					stdin, stdoutline, stderr = ssh.exec_command('cd %s/%s;python /share/Onc_SmallPanel/capsmart_31gene/addReportTime.py %s ;grep reportSimple_production run.sh> aaa.sh;/bin/sh aaa.sh; /share/public/software/Onc_Soft/python/2.7.14/bin/python /share/Onc_SmallPanel/capsmart_31gene/mail_capsmart_new_v3.py -f %s %s.xls;echo finish > finish.txt'%(root,dir,time_sheet,dir,dir))
					mail31_err = stderr.read()
					if mail31_err != '':
						youjian = 'hi,\n\n该批次已完成,但是在发邮件时出了问题,请检查文件名是否正确:\ncd %s/%s;\npython /share/Onc_SmallPanel/capsmart_31gene/addReportTime.py %s\n;/share/public/software/Onc_Soft/python/2.7.14/bin/python /share/Onc_SmallPanel/capsmart_31gene/mail_capsmart_new.py -f %s %s.xls\nError内容:\n%s'%(root,dir,time_sheet,dir,dir,mail31_err)
						subject = '31gene报错'
						mail_1(dir=dir,to_list = ['xingh3223@berryoncology.com'],cc_list = [],subject=subject,main_text=youjian)

				elif status_done == True and hebing_flag2 == True and germline_flag2 == True and re.search(r'86gene',root) and pair_txt_list_flag == True and somatic_flag3 == True:
					ssh = paramiko.SSHClient()
					ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
					ssh.connect(hostname='10.100.2.10',port=22,username='capsmart',password='berry2012')
					stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
					stdin, stdoutline, stderr = ssh.exec_command("export LANG=zh_CN.UTF-8;export LC_ALL=en_GB.UTF-8;cd %s/%s;/share/work1/xinghe/software/anaconda3/bin/python /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.all.py"%(root,dir))
				elif status_done == True and hebing_flag2 == True and germline_flag2 == True and re.search(r'654gene',root) and pair_txt_list_flag == True and somatic_flag3 == True:
					ssh = paramiko.SSHClient()
					ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
					ssh.connect(hostname='10.100.2.10',port=22,username='capsmart',password='berry2012')
					stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
					stdin, stdoutline, stderr = ssh.exec_command("export LANG=zh_CN.UTF-8;export LC_ALL=en_GB.UTF-8;cd %s/%s;/share/work1/xinghe/software/anaconda3/bin/python /share/work1/xinghe/proc/captools/refresh_result/refresh.v1.all.py"%(root,dir))


def get_NM_dic():
	NM_list = ['/share/work3/capsmart/pipeline/capSMART/CAPcSMART/bed/gene_NM.list','/share/public/database/Gynecological_cancer_backup/GTDB/MainNM_yrt_20191111.txt','/share/work1/wulj/database/svdb/rna/anndb/gene2trs.tsv']
	NM_dic1 = {}
	NM_dic2 = {}
	NM_dic3 = {}
	i = 0
	for each in NM_list:
		i +=1
		for line in open(each,'r'):
			gene,NM = line.strip().split('\t')[0:2]
			dic_id = 'NM_dic'+str(i)
			locals()[dic_id][gene] = NM	
	return NM_dic1, NM_dic2, NM_dic3

NM_dic1, NM_dic2, NM_dic3 = get_NM_dic()	


roots = ['/share/Onc_SmallPanel/capsmart_86gene','/share/Onc_LargePanel/capsmart_457gene','/share/Onc_LargePanel/capsmart_654gene','/share/Onc_SmallPanel/capsmart_86gene/KYONF2019501','/share/Onc_SmallPanel/capsmart_31gene','/share/Onc_SmallPanel/capsmart_31gene/LC_CRC','/share/Onc_SmallPanel/WES_Plus']
for root in roots:
	if re.search(r'86gene',root):
		panel = '86gene'
	elif re.search(r'457gene',root):
		panel = '457gene'
	elif re.search(r'654gene',root):
		panel = '654gene'
	elif re.search(r'31gene',root) and not re.search('/LC_CRC',root):
		panel = '31gene'
	elif re.search(r'31gene',root) and re.search('/LC_CRC',root):
		panel = '13gene'
	elif re.search(r'WES_Plus',root):
		panel = 'WESPlus'
	else:
		panel = 'unknow_panel'
	multi_check(root,panel)
