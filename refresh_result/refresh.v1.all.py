#!/usr/bin/python
#coding=utf-8
import os
import glob
import sys
import re
import sys
import smtplib,ssl
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.header import Header
from email.utils import formatdate
from email import encoders
import xlrd
from collections import OrderedDict
import argparse
import subprocess

parser = argparse.ArgumentParser(prog="",description="过滤需要删除的位点")
#parser.add_argument("-c","--cnv",help='要删除的cnv位点所在的行')
parser.add_argument("-m","--merge",help='如果进行了合并操作，那么就添加这个参数')
parser.add_argument("-d","--dir",help='要重新生成的目录,如200111_Novaseq_A00838_A_AH327GDSXY_csmart2_0_SPM_v20_86_auto')
parser.add_argument("-n","--no",help='只发送邮件不重新刷新脚本')
args = parser.parse_args()

#reportSimple_production = '/share/work2/liuyan/git_code/capSMART2.0/csmart/report/scripts/reportSimple_production.v3.3.py'
addReportTime = os.path.dirname(os.getcwd())+'/addReportTime.py'
reportSimple_production = 'reportSimple_production'
getSNVIndel = 'getSNVIndel.v2.4.py'
getSNVIndel_right_align = 'getSNVIndel.right.align'
mergeSomaticGermline = 'mergeSomaticGermline'
getTMB = 'getTMB'
TMB_density = 'TMB-density.R .* -project TMB.xls'
TMB_rank_copy = 'TMB_rank_copy'
getSomaticInterpretationInfo = 'getSomaticInterpretationInfo'

def get_refresh_cmds(panel):

	if panel == '86':
		p1 = subprocess.Popen('grep %s ./data/report.sh'%(getSNVIndel),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
		p2 = subprocess.Popen('grep %s ./data/report.sh'%(getSNVIndel_right_align),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
		p3 = subprocess.Popen('grep %s ./data/report.sh'%(getSomaticInterpretationInfo),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
		
		p4 = subprocess.Popen('grep %s ./data/report.sh'%(mergeSomaticGermline),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
		p5 = subprocess.Popen('grep %s run.sh'%(reportSimple_production),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

		s1 = [m1 for m1 in p1.stdout]				
		s2 = [m2 for m2 in p2.stdout]				
		s3 = [m3 for m3 in p3.stdout]				
		s4 = [m4 for m4 in p4.stdout]				
		s5 = [m5 for m5 in p5.stdout]				
		
		if len(s1) > 1 or len(s2) > 1 or len(s3) > 1 or len(s4) > 1 or len(s5)>1:
			sys.exit('可能grep出多行命令..............')
		elif len(s1) == 0 or len(s2) == 0 or len(s3) == 0 or len(s4) == 0 or len(s5) == 0:
			sys.exit('可能有命令没有搜索到................')
		elif len(s1) == 1 and len(s2) == 1 and len(s3) == 1 and len(s4) == 1 and len(s5) ==1:
			cmd1 = '\n'.join(s1 + s2 + s3 + s4)
			cmd2 = '\n'.join(s5)
			return [cmd1, cmd2 ]

	elif panel == '457' or panel == '654':
		p1 = subprocess.Popen('grep %s ./data/report.sh'%(getSNVIndel),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(第一句刷新somatic_snv_indel.xls)
		p2 = subprocess.Popen('grep %s ./data/report.sh'%(getSNVIndel_right_align),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(第二句刷新somatic_snv_indel.right_aln.xls)
		p3 = subprocess.Popen('grep %s ./data/report.sh'%(getTMB),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(刷新TMB值)
		p4 = subprocess.Popen("grep -E '%s' ./data/report.sh"%(TMB_density),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(刷新TMB图)
		p5 = subprocess.Popen('grep %s ./data/report.sh'%(getSomaticInterpretationInfo),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
		p6 = subprocess.Popen('grep %s ./data/report.sh'%(mergeSomaticGermline),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(刷新SnvIndel.xls)
		p7 = subprocess.Popen('grep %s ./data/transfer_data.sh'%(TMB_rank_copy),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True) ##(重新上传TMB图)(可放最后上传))
		p8 = subprocess.Popen('grep %s run.sh'%(reportSimple_production),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)

		s1 = [m1 for m1 in p1.stdout]				
		s2 = [m2 for m2 in p2.stdout]				
		s3 = [m3 for m3 in p3.stdout]				
		s4 = [m4 for m4 in p4.stdout]				
		s5 = [m5 for m5 in p5.stdout]				
		s6 = [m6 for m6 in p6.stdout]				
		s7 = [m7 for m7 in p7.stdout]				
		s8 = [m8 for m8 in p8.stdout]				
		
		if len(s1) > 1 or len(s2) > 1 or len(s3) > 1 or len(s4) > 1 or len(s5) > 1 or len(s6) > 1 or len(s7) > 1 or len(s8)>1:
			#print(s1,s2,s3,s4,s5,s6,s6,s8)
			for t in [s1,s2,s3,s4,s5,s6,s6,s8]:
				if len(t)>1:
					print(t)
			sys.exit('可能grep出多行命令..............')
		elif len(s1) == 0 or len(s2) == 0 or len(s3) == 0 or len(s4) == 0  or len(s5) == 0 or len(s6) == 0 or len(s7) == 0 or len(s8)==0:
			sys.exit('可能有命令没有搜索到................')
		elif len(s1) == 1 and len(s2) == 1 and len(s3) == 1 and len(s4) == 1 and len(s5) == 1 and len(s6) == 1 and len(s7) == 1 and len(s8) == 1:
			#cmd1 = '\n'.join(s1 + s2 + s3 + s4 + s5 + s6 + s7)
			cmd1 = '\n'.join(s1 + s2 + s3 + s4  + s6 + s7)
			cmd2 = '\n'.join(s8)
			return [cmd1, cmd2 ]
	elif panel == '83':
		return ['No cmd','No cmd']
	else:
		sys.exit('Error:没有这个panel........')



def get_QC(xls,sheet,msg=''):

	rw = xlrd.open_workbook(auto_xls)
	rs = rw.sheet_by_name(sheet)
	nrow = rs.nrows
	ncol = rs.ncols
	content = []
	content_string = '<p>{msg}</p><table style="font-size:13px" border="1" cellspacing="0" cellpadding="0" white-space:nowrap>'.format(**locals())
	for i in range(nrow):
		tmp = ''
		#for j in [38,39,40,42,44]:
		for j in [40,41,42,44,46]:
			each_cell = str(rs.cell_value(i,j))
			mei_hang = str('\t'.join(map(str,rs.row_values(i))))
			if (re.search(r'^NOPASS',each_cell) and not re.search(r'NC_T3_V2_G86|NC_T_B19_4|NC_T_G654',mei_hang)) or (each_cell == 'high_contaminated' and not re.search(r'NC_T3_V2_G86|NC_T_B19_4|NC_T_G654',mei_hang)) or (each_cell == 'mid_contaminated' and not re.search(r'NC_T3_V2_G86|NC_T_B19_4|NC_T_G654',mei_hang)):
				val = '<td bgcolor="#FF0000" >'+each_cell+'</td>'
			else:
				val = '<td >'+each_cell+'</td>'
			tmp +=val
		tmp = '<tr>'+tmp+'</tr>'
		content_string += tmp
	content_string+='</table>'
	return content_string



def get_Pscore(xls,msg=''):

	content = []
	f = open(xls,'r')
	content_string = '<p>{msg}</p><table border="1" cellspacing="0" cellpadding="0">'.format(**locals())
	for i in f:
		arr = i.strip().split('\t')
		tmp = ''
		if re.search(r'^Silico',arr[0]):
			continue
		for j in arr:
			each_cell = str(j)
			if (re.search(r'^NOPASS',each_cell) and not re.search(r'NC_T3_V2_G86|NC_T_G654',each_cell)): 
			#or each_cell == 'high_contaminated' or each_cell == 'mid_contaminated' or each_cell == 'mismatch':
				val = '<td bgcolor="#FF0000">'+each_cell+'</td>'
			#elif re.search(r'^NOPASS',each_cell) or each_cell == 'low_contaminated':
			#	val = '<td bgcolor="#FFA07A">'+each_cell+'</td>'
			else:
				val = '<td>'+each_cell+'</td>'
			tmp +=val
		tmp = '<tr>'+tmp+'</tr>'
		content_string += tmp
	content_string+='</table>'
	return content_string


def to_mail(panel,QC_string,test=True):

	mail_host = 'mail.berryoncology.com'
	mail_user = 'xingh3223'
	mail_pass = 'iiis_1234'
	sender = 'xingh3223@berryoncology.com'

	#mail_list:
	if test == True:
		to_list = ['xingh3223@berryoncology.com']
		cc_list = []
	else:
		print('panel..........................')
		print(panel)
		if panel == '86' and not re.search(r'KYONF2019501',os.getcwd()):
			to_list = ['oncology_labreport@berryoncology.com']
			cc_list = ['oncology_ycfxs@berryoncology.com','baijian488@berryoncology.com','jiangdzh3403@berryoncology.com','liuyan725@berryoncology.com','lilinwei001@berryoncology.com','donghansheng041@berryoncology.com','dingj3639@berryoncology.com','oncology_qau@berryoncology.com','zhangwj3075@berryoncology.com','tangcong935@berryoncology.com','xufl3252@berryoncology.com','wangluxi746@berryoncology.com','wanghu894@berryoncology.com','xiliying196@berryoncology.com','zhaoyuanyuan@berryoncology.com','rulanlan@berryoncology.com','fulf3657@berryoncology.com','wun3623@berryoncology.com','wangzx3872@berryoncology.com','wujq3870@berryoncology.com','guoyushuai712@berryoncology.com','liull4492@berryoncology.com','yangrutao796@berryoncology.com','zhangp4456@berryoncology.com','xingh3223@berryoncology.com']
		elif panel == '83' and re.search(r'KYONF2019501',os.getcwd()):
			to_list = ['oncology_labreport@berryoncology.com']
			cc_list = ['oncology_ycfxs@berryoncology.com','baijian488@berryoncology.com','jiangdzh3403@berryoncology.com','wushx3173@berryoncology.com','wun3623@berryoncology.com','yangrutao796@berryoncology.com','zhangp4456@berryoncology.com','xingh3223@berryoncology.com']
		elif panel == '457':	
			to_list = ['oncology_labreport@berryoncology.com']
			cc_list = ['oncology_ycfxs@berryoncology.com','baijian488@berryoncology.com','jiangdzh3403@berryoncology.com','liuyan725@berryoncology.com','lilinwei001@berryoncology.com','donghansheng041@berryoncology.com','dingj3639@berryoncology.com','oncology_qau@berryoncology.com','zhangwj3075@berryoncology.com','tangcong935@berryoncology.com','xufl3252@berryoncology.com','wangluxi746@berryoncology.com','wanghu894@berryoncology.com','xiliying196@berryoncology.com','zhaoyuanyuan@berryoncology.com','rulanlan@berryoncology.com','fulf3657@berryoncology.com','wun3623@berryoncology.com','wangzx3872@berryoncology.com','wujq3870@berryoncology.com','liull4492@berryoncology.com','yangrutao796@berryoncology.com','zhangp4456@berryoncology.com','xingh3223@berryoncology.com']
		elif panel == '654':
			to_list = ['oncology_labreport@berryoncology.com']
			cc_list = ['oncology_ycfxs@berryoncology.com','baijian488@berryoncology.com','jiangdzh3403@berryoncology.com','liuyan725@berryoncology.com','lilinwei001@berryoncology.com','donghansheng041@berryoncology.com','dingj3639@berryoncology.com','oncology_qau@berryoncology.com','zhangwj3075@berryoncology.com','tangcong935@berryoncology.com','xufl3252@berryoncology.com','wangluxi746@berryoncology.com','wanghu894@berryoncology.com','xiliying196@berryoncology.com','zhaoyuanyuan@berryoncology.com','rulanlan@berryoncology.com','wujq3870@berryoncology.com','wangzx3872@berryoncology.com','guoyushuai712@berryoncology.com','liuk3940@berryoncology.com','wun3623@berryoncology.com','liull4492@berryoncology.com','yangrutao796@berryoncology.com','zhangp4456@berryoncology.com','xingh3223@berryoncology.com']
	
	receivers = to_list+cc_list
	
	#mailed files list:
	if panel == '86' :	
		mailedfiles = glob.glob(bz666dir+'/*.zip')+glob.glob(auto_xls)
		print(mailedfiles)
	elif panel == '83' :	
		mailedfiles = glob.glob(bz666dir+'/*.zip')+glob.glob(auto_xls) + glob.glob('./for_negative.zip')
		print(mailedfiles)
	elif panel == '457':	
		os.system('/usr/bin/zip TMB_rank.zip ./TMB_rank/*.png')
		mailedfiles = glob.glob(bz666dir+'/*.zip')+glob.glob(auto_xls)+glob.glob('TMB_rank.zip')
		print(mailedfiles)
	elif panel == '654':
		pair_num = len([line.strip() for line in open('pair.txt') if line.startswith('Silico')])
		if pair_num >= 0:
			os.system('/usr/bin/zip TMB_rank.zip ./TMB_rank/*.png')
			os.system('/usr/bin/zip cnv.zip ./cnv/*/*.png')
			mailedfiles = glob.glob(bz666dir+'/*.zip')+glob.glob(auto_xls)+glob.glob('TMB_rank.zip')+glob.glob('./cnv.zip')
		else:
			mailedfiles = glob.glob(bz666dir+'/*.zip')+glob.glob(auto_xls)+glob.glob('./TMB_rank/*.png')+glob.glob('./cnv/*/*.png')

		print(mailedfiles)
	else:
		sys.exit('Error:please choose panel: 86/83/457/654.....')

	qc_file = 'QC.xls.new'
	buhege = subprocess.Popen('''grep 'NOPASS:Fraction of target covered at least 0.2 x Average depth' %s|cut -f1'''%(qc_file),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
	buhege_list = ['for_negative/'+line.strip()+'_for_negative/'+line.strip()+'_hotspots.xls' for line in buhege.stdout]
		#for_negative/Z19N00088-B1TA_for_negative/Z19N00088-B1TA_hotspots.xls
	if buhege_list != []:
		mailedfiles.extend(buhege_list)

	######构造邮件
	message = MIMEMultipart()
	message['From'] = 'xingh3223@berryoncology.com'
	message['To'] = ';'.join(to_list)
	message['Cc'] = ';'.join(cc_list)
	message['Date'] = formatdate(localtime=True)
	subject = '[%sgene]%s report'%(panel,subdir)
	message['Subject'] = Header(subject,'utf-8')
	
	message.attach(MIMEText(QC_string, 'html', 'utf-8'))
	for each in mailedfiles:
		print(each)
		filename = os.path.basename(each)
		part1 = MIMEBase('application','octet-stream')
		part1.set_payload(open(each,'rb').read())
		encoders.encode_base64(part1)
		part1.add_header('Content-Disposition','attachment; filename="%s"'%filename)
		message.attach(part1)
	
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

	if panel == '457' or panel == '654':
		try:
			os.system('rm TMB_rank.zip')
			os.system('rm cnv.zip')
		except:
			pass


def get_qc_other(qc_file='QC_other.xls.new'):

	p_qc = subprocess.Popen('grep NOPASS %s'%(qc_file),shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
	nopass_white = []
	for line in p_qc.stdout:
		if re.search(r'NOPASS',line):
			nopass_white.append(line.strip().split('\t')[0])	
	qc_other_nopass = '\n'.join(nopass_white)
	if qc_other_nopass == '':
		qc_other_nopass = 'QC_other:PASS'
	else:
		qc_other_nopass = 'QC_other:NOPASS samples:\n\n' + qc_other_nopass+'\n'
	return qc_other_nopass

#parse arguments:
if args.dir == None:
	path = os.path.dirname(os.getcwd())
	subdir = os.path.basename(os.getcwd())
else:
	path = os.path.dirname(args.dir) 
	subdir = os.path.basename(args.dir)

if re.search(r'86gene',path) and not re.search(r'KYONF2019501',os.getcwd()):
	panel = '86'
elif re.search(r'86gene',path) and re.search(r'KYONF2019501',os.getcwd()):
	panel = '83'
elif re.search(r'457gene',path):
	panel = '457'
elif re.search(r'654gene',path):
	panel = '654'
else:
	sys.exit('请填写panel:86/457/654..........')


print('Panel:'+panel+'\n---------------------------------')
#subdir2 = re.search(r'(.*)(_auto|_info|_sup)',subdir).group(1)
subdir1 = subdir.replace('_info','').replace('_sup','')
subdir2 = subdir.replace('_auto','').replace('_info','').replace('_sup','')
analysis_dir = path+'/'+subdir
#auto_xls = './'+subdir2+'_auto.xls'  #最后的结果文件
auto_xls = './'+subdir1+'.xls'  #最后的结果文件
Pscore_xls = './Pscore.fmt.xls'  #最后的结果文件
time_sheet = analysis_dir+'/'+subdir2+'.time.xls'
maildir = analysis_dir+'/Mailed'
maildir_bz666 = maildir+'/bz666'
bz666dir = analysis_dir+'/bz666'

print('\nBackup important files\n--------------------------------------')
bak_num = len(glob.glob('%s/backup*'%(analysis_dir))) + 1
backupdir = analysis_dir + '/backup_' + str(bak_num)
print('mkdir -p %s;cp -r %s %s;cp %s/* %s'%(backupdir,bz666dir,backupdir,analysis_dir,backupdir))
os.system('mkdir -p %s;cp -r %s %s;cp %s/* %s'%(backupdir,bz666dir,backupdir,analysis_dir,backupdir))
os.system('/share/public/software/Python-2.7.13/bin/python %s %s'%(addReportTime,time_sheet))  #在time.xls中加上report_complete_time这一列

#get cmds that need refresh
cmd1,cmd2 = get_refresh_cmds(panel=panel)
if args.merge == None:
	pass 
else:
	print('\nClose sites was merged, so it is refreshing these files\n-------------------------------------')
	print(cmd1)
	os.system(cmd1) #if close sites was merged
if args.no == None:
	print('\nRefreshing final report\n----------------------------------------')
	print(cmd2)
	os.system(cmd2)
	print('\nSending mail\n-------------------------------------')
else:
	print('\nNo refreshed commands, just Sending mail\n-------------------------------------')
###to mail
QC_string = get_QC(xls=auto_xls,sheet='QC',msg='\n\n@报告组,\n\n附件是%sgene分析结果。请核对结果，如有问题请及时反馈。QC:\n'%(panel))
QC_string = QC_string+get_qc_other()
to_mail(panel=panel,QC_string=QC_string,test=False)


