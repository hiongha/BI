#!/usr/bin/python
#coding=utf-8
import sys
import re
import os
import time
import glob
import subprocess
import smtplib,ssl
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.header import Header
from email.utils import formatdate
from email import encoders
from collections import OrderedDict
import paramiko
import configparser

reload(sys)
sys.setdefaultencoding('utf-8')


def get_NM_dic():
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


def panel_check(root):
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
	return panel


def get_list(line):
	new_list = line.strip().split(',')
	return new_list


def mailer(dir,subject,main_text,formatting='plain',**kwargs):
	'''Just a mail box.'''
	for k,v in kwargs.items():
		exec("%s"%(k) + "=" + "'%s'"%(str(v)))
	to_list = to_list.strip().split(',')
	cc_list = cc_list.strip().split(',')
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
		print("Send mail successfully!")
	except smtplib.SMTPException:
		print("Error:can't send mail normally!")


def rewrite(FC,done_file):
	'''
	write dealed flowcell number into '*.done.txt'
	'''
	with open(done_file,'a') as F2:
		F2.write(FC+'\n')


def get_undone_fc(done_file,typ):
	'''
	get directories that germline result were not checked.
	'''
	done_fc = []
	with open(done_file) as F:
		for line in F:
			done_fc.append(line.strip())	
	dirs = glob.glob(root+'/*')
	dirs = map(os.path.basename,dirs)
	undone_fc = set(dirs)-set(done_fc)
	return undone_fc


def check_Somatic_interpretation(path):
	'''check if somatic sites exists.'''
	somatic_file = glob.glob(path+'/Somatic_interpretation/tumors_Clisig_reviewed.txt')	
	Tumor_check_list = []
	somatic_finish_flag = False	
	somatic_check_flag = False
	Tumor_check_line = ''
	Tumor_check_num = 0
	if somatic_file == []:
		somatic_finish_flag = False

	elif somatic_file != []:
		somatic_file = somatic_file[0]
		p = subprocess.Popen("grep 共计审核 %s"%(somatic_file),shell=True,
							stdin=subprocess.PIPE,stdout=subprocess.PIPE,
							stderr=subprocess.PIPE,universal_newlines=True)
		shifouwancheng = [m.strip() for m in p.stdout]
		if len(shifouwancheng) == 0:
			somatic_finish_flag = False
		else:
			somatic_finish_flag = True
			with open(somatic_file,'r') as F:
				content = F.read()
				parts = content.split('--------------------')
		
			for each in parts:
				Tumor_check_list.append(each)
			if len(Tumor_check_list) == 0:
				somatic_check_flag = False
				Tumor_check_line = ''
			elif len(Tumor_check_list) >=1:
				somatic_check_flag = True
				Tumor_check_line = '--------------------'.join(Tumor_check_list)
				Tumor_check_num = len(Tumor_check_list) - 1 
	return somatic_finish_flag,somatic_check_flag,Tumor_check_line,Tumor_check_num


def QC_check(path):
	'''get QC result.'''
	file = path + '/QC.xls.new'
	check_cols = check_cols_ini
	QC_check_dic = OrderedDict()
	if not os.path.exists(file):
		pass
	else:
		with open(file,'r') as F:
			F.readline()
			for line in F:
				arr = line.strip().split('\t')
				check,Pairs,Contaminated,Contamination_level,QC_other_check =  [arr[i] for i in check_cols]	
				Nsamp,Psamp = Pairs.split(':')			
				QC_check_dic[Psamp] = [check,Pairs,Contaminated,Contamination_level,QC_other_check]		
	return QC_check_dic


def write_to_html(List):
	'''wirte MSI results into html format'''
	css = '''
			<style type="text/css">
				#msi{
					border-collapse:collapse;
				}
				#msi th,#msi tr,#msi td {
					font-size:12px;
					border:1px solid #a1a1a1;
					white-space:nowrap;
				}
			</style>
			'''
	html = css
	html += ''' <p>Hi all,</p>
				</br>
				<p>MSI results:</p>
				<table id="msi">'''
	for each in List:
		html += '<tr>'
		for i in each:
			html +='<td> %s </td>'%(str(i))
		html += '</tr>'
	html += '</table>'
	return html


def msi_check(path):
	'''
	check MSI result.		
	'''
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


def get_final_qsub(all_status_files):
	'''
	sometimes the job would failed, and was qsub to clusters, so we need to find the final status file.	
	'''
	max = 0
	max_dir = all_status_files[0]
	for i in all_status_files:
		if len(i.strip().split('.status')) > max:
			max_dir = i.strip()
			max = len(i.strip().split('.status'))
	return max_dir 



def multi_check(root,panel):
	'''
	check pair/germline/analysis status.

	'''

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

				p = subprocess.Popen("grep 共计审核 %s"%(germline),shell=True,stdin=subprocess.PIPE,
									stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

				content_grep = ''.join([i.strip() for i in p.stdout])
				count = len(content_grep)
				if count == 0:
					pass
				else:
					shenhe_count = re.search(r'共计审核(\d+)个突变位点',content_grep).group(1)
					if int(shenhe_count) == 0 and not (somatic_finish_flag==True and somatic_check_flag==True):
						germ_receiver = default_receiver
					else:
						if re.search(r'86gene|654gene|457gene|31gene|WES_Plus',root):
							germ_receiver = germline_receiver
					content = []
					flag = True
					with open(germline,'r') as F1:
						for each in F1:
							if flag == True:
								content.append(each)
					content = ''.join(content)
	
					pF = open('%s/%s/pair.txt'%(root,dir),'r')
					sample_count = str(len(pF.readlines()))

					subject = '%s-%s审-%s'%(area,panel.strip('gene'),dir)
					main_text = '\n\nHi ,\n\n\tPlease check germline results of %s. \n\n\t 共%s对样本,%s个位点. \
								\n\n\tData path: %s/%s\n\n\n%s\n\n ====================================\n\n【Somatic Check】:\
								\n\n%s'%(dir,sample_count,str(int(shenhe_count)+Tumor_check_num),root,dir,content,Tumor_check_line)
					if re.search(r'86gene|654gene|457gene|31gene|WES_Plus',root):
						mailer(dir=dir,subject=subject,main_text=main_text,**dict(germline_sender,**germ_receiver))
					rewrite(FC=dir,done_file=type_done_file)

		if typ == 'cnv':
			for dir in undone_fc:
				CNV_file = root+'/%s/CNV.xls'%(dir)
				if not os.path.exists(CNV_file):
					continue
				subject = '%s-%sCNV-%s'%(area,panel.strip('gene'),dir)
				main_text = '\n\nHi, \n\n\t%s CNV analysis is over, please check it. \n\tanalysis path%s/%s\n\n\n '%(dir,root,dir)
				cc_list = []
				mailer(dir=dir,subject=subject,main_text=main_text,**default_sender.update(default_revceiver))
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
				p = subprocess.Popen("grep no %s"%(pair_file),shell=True,stdin=subprocess.PIPE,
									stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				sample_no_pair = [i.strip() for i in p.stdout]
				count_pair = len(sample_no_pair)
				if count_pair == 0:
					subject = '%s-%s配-%s NOT EXISTS'%(area,panel.strip('gene'),dir)
					main_text = '\n\nHi,\n\n\t%s batch NOT EXISTS non-matched samples.\n\npair.txt:\
								\n%s\tData path is %s/%s\n\n\n '%(dir,pair_txt_content,root,dir)
					mailer(dir=dir,subject=subject,main_text=main_text,**dict(default_sender,**default_receiver))
				else:
					subject = '%s-%s配-%s EXISTS'%(area,panel.strip('gene'),dir)
					main_text = '\n\nHi sample match result,\n\n\t%s batch EXISTS non-matched samples \n%s.\n\npair.txt:\
								\n%sData path is %s/%s\n\n\n '%(dir,'\n'.join(sample_no_pair),pair_txt_content,root,dir)
					mailer(dir=dir,subject=subject,main_text=main_text,**dict(default_sender,**default_receiver))
				rewrite(FC=dir,done_file=type_done_file)

		if typ == 'finish':
			for dir in undone_fc:
				warn_file = root+'/%s/warnClosePos.txt'%(dir)
				final_result_file = root+'/%s/%s.xls'%(dir,dir)
				dir2 = dir.replace('_sup','').replace('_LC','').replace('_CRC','').replace('_info','').replace('_auto','').replace('_info','').replace('_sup','')
				
				all_status_files = glob.glob(root+'/%s/%s.job.*status'%(dir,dir2))
				if len(all_status_files) < 3:
					continue
				elif len(all_status_files) == 3:
					final_qsub = root+'/%s/%s.job.status.status.status'%(dir,dir2)
				elif len(all_status_files) > 3:
					final_qsub = get_final_qsub(all_status_files)

				time_file = root +'/'+dir +'/' + dir2.replace('_sup','').replace('_auto','')+'.time.xls'
				if not (os.path.exists(final_qsub) and os.path.exists(time_file)):
					continue

				print('\n------------------------------------------------------------------------------')
				print("New Analysis Directory: \n" + root + '/' + dir2)
				#flag1 status_done: if job.status was finished.	
				p3 = subprocess.Popen("grep status %s"%(final_qsub),shell=True,stdin=subprocess.PIPE,
									stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				status_done = True
				for eachstatus in p3.stdout:
					if eachstatus.strip() != 'status done':
						status_done = False
						break
				if status_done == False:
					continue
				no_dce_file = glob.glob(root+'/'+dir+'/*.no_dce.xls')
				if re.search(r'LC_CRC',root) and no_dce_file == []:
					continue
				else:
					time.sleep(2)

				time.sleep(2)	
	
				#flag2 germline_flag2: whether exists germline sites that need to be checked.	
				germline_flag2 = False
				germline2 = root+'/%s/germline/normals_germline_reviewed.txt'%(dir)
				p2 = subprocess.Popen("grep 共计审核 %s"%(germline2),shell=True,stdin=subprocess.PIPE,
										stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
				content_grep2 = ''.join([i.strip() for i in p2.stdout])
				if re.search(r'共计审核0个突变位点',content_grep2):
					germline_flag2 = True

				#flag3 somatic_check_flag2: whether somatic sites exists. 
				somatic_finish_flag2,somatic_check_flag2,Tumor_check_line,Tumor_check_num = check_Somatic_interpretation(root+'/'+dir)

				#flag4 somatic_flag3: whether somatic sites exists and somatic analysis is finished or not.
				somatic_flag3 = False
				if somatic_finish_flag2 == True and somatic_check_flag2 == False:
					somatic_flag3 = True

				#flag5 total_close_sites: check warnClose.txt, if sites exist, then create samtools command and mail it.	
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
				   
				if True:					
					for line in need_checked_list:
						if line.strip() != '':
							array = line.strip().split('\t')
							total_close_sites.append(line.strip())
							samp,chr,start = array[0:3]
							this_gene = array[8].split(':')[0]
							zhushi_gene = NM_dic1.get(this_gene,'--') +'/'+ NM_dic2.get(this_gene,'--') +'/'+ NM_dic3.get(this_gene,'--') 
							pos_list.append(int(start))
							min_pos = str(min(pos_list)-1)
						elif line.strip() == '' :
							samcmd = '\n%s tview %s/%s/./aln/%s/*.cons.bam --reference %s -p %s:%s|grep -i [1atgc*]|les\
									\n\nMainTranscript:%s\n\n--------------------------\n\n'%(samtools,root,dir,samp,reference,chr,min_pos,zhushi_gene)
							total_close_sites.append(samcmd)
							pos_list = []
							samp = '';chr='';start='';

				time_sheet = dir2.replace('_sup','').replace('_auto','')+'.time.xls'
				if total_close_sites == []:
					mail_shell = '【refresh after merge】:\ncd %s/%s;\n%s %s;\n%s %s;\n'%(root,dir,python3,refreshtestpy,python3,refreshallpy)
				else:
					mail_shell = '【refresh after merge】:\ncd %s/%s;\n%s %s -merge;\n%s %s -merge;\n'%(root,dir,python3,refreshtestpy,python3,refreshallpy)
	
				warn_text = '\n'.join(total_close_sites) + '\n' + mail_shell	
				subject = '%s-%s完-%s'%(area,panel.strip('gene'),dir)
				main_text = '\n\nHi,\n\n\tAnalysis of %s is over \
							\n\n\tData path is %s/%s\n\n【merging commands】:\n\n%s\n '%(dir,root,dir,warn_text)
				mailer(dir=dir,subject=subject,main_text=main_text,**dict(default_sender,**default_receiver))

				#flag6 msi_check: mail the MSI result to report group.
				if panel == '654gene' and re.search(r'_654T_',dir):
					msi_html = msi_check(root+'/'+dir)
					subject = '%s-%sMSI-%s '%(area,panel.strip('gene'),dir)
					mailer(dir=dir,subject=subject,main_text=msi_html,formatting='html',**dict(default_sender,**msi_receiver))

			
				rewrite(FC=dir,done_file=type_done_file)
				
				#flag7 hebing_flag2: check if there are sites that need to be merged.
				hebing_flag2 = False
				if total_close_sites == []:
					hebing_flag2 = True
			
				#flag8 pair_txt_list_flag: whether only Silico sample exists.
				pair_txt = root+'/%s/pair.txt'%(dir)
				pair_txt_list_flag = True
				if os.path.exists(pair_txt):	
					pair_txt_list = [line.strip() for line in open(pair_txt,'r') if not re.search(r'^Sili',line)]
				else:
					pair_txt_list = []
				if pair_txt_list == []:
					pair_txt_list_flag = False

				if status_done == True and hebing_flag2 == True and germline_flag2 == True \
							and pair_txt_list_flag == True and somatic_flag3 == True:
					ssh = paramiko.SSHClient()
					ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
					ssh.connect(hostname='10.100.2.10',port=22,username='capsmart',password='berry2012')
					if re.search(r'_31genexxx',root):
						stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
						stdin, stdoutline, stderr = ssh.exec_command('cd %s/%s;python %s %s ;grep reportSimple_production run.sh> aaa.sh;\
										/bin/sh aaa.sh; %s %s -f %s %s.xls;echo finish > finish.txt'%(root,dir,addreporttime,time_sheet,python2,mail_capsmart,dir,dir))
					elif re.search(r'86gene|654gene',root):
						stdin,stdout,stderr = ssh.exec_command('ls %s/%s'%(root,dir))
						stdin, stdoutline, stderr = ssh.exec_command("export LANG=zh_CN.UTF-8;export LC_ALL=en_GB.UTF-8;cd %s/%s;%s %s"%(root,dir,anaconda3,refreshallpy))


'''make configures '''
if __name__ == '__main__':
	file = os.path.split(os.path.realpath(__file__))[0]
	config_file = file + '/config.ini'
	config = configparser.ConfigParser()
	config.read(config_file)
	solo_path = dict(config['path'].items())
	prog = dict(config['prog'].items())
	globals().update(solo_path)
	globals().update(prog)
	NM_list = get_list(nm_list)
	roots = get_list(roots)
	#senders
	default_sender = config['default_sender']
	germline_sender = config['germline_sender']
	report_sender = config['report_sender']
	#receivers
	default_receiver = config['default_receiver']
	germline_receiver = config['germline_receiver']
	msi_receiver = config['msi_receiver']
	#others
	check_cols_ini = [int(i) for i in get_list(config['qc_check']['check_cols_ini'])]

	NM_dic1, NM_dic2, NM_dic3 = get_NM_dic()

	for root in roots:
		panel = panel_check(root)
		multi_check(root,panel)


