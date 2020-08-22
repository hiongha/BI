#!/usr/bin/python
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>

import argparse
import sys
import smtplib,ssl
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.header import Header
from email.utils import formatdate
from email import encoders
from collections import OrderedDict 

parser=argparse.ArgumentParser(prog="",description="将北肿QC结果发送到xingh3223@berryoncology.com")
parser.add_argument("-m","--message-from-check",help="check_jobs.sh脚本中传出的信息",default='everything goes normally')
parser.add_argument("-s","--sample",help="sample name",required = True)
parser.add_argument("-w","--work-dir",help="样品分析路径",required = True)
parser.add_argument("-f","--uploaded-file",help="邮件附件",required = True)
args=parser.parse_args()

message_from_check = args.message_from_check
sample = args.sample
root_dir = args.work_dir
uploaded_file_tmp = args.uploaded_file.strip().split(',')
report_dir = root_dir + '/%s_report'%(sample)
upload_files = OrderedDict()
for m in uploaded_file_tmp:
	m_long = report_dir +'/'+ m
	upload_files[m] = m_long
	print(m)
	print(m_long)

mail_host = 'mail.berryoncology.com'
mail_user = 'xingh3223'
mail_pass = 'iiis_1234'

sender = 'xingh3223@berryoncology.com'
receivers = ['xingh3223@berryoncology.com','974110177@qq.com']
#cc = ['zoujianing911@berryoncology.com','weiyuyao001@berryoncology.com','suzhq3421@berryoncology.com','niewd3409@berryoncology.com','yujun001@berryoncology.com','wangzhiqiong915@berryoncology.com','chenxiaoyan712@berryoncology.com','suxiaoxing001@berryoncology.com','sunfuming615@berryoncology.com','baijian488@berryoncology.com','xufl3252@berryoncology.com','oncology_labreport@berryoncology.com']
cc = ['oncology_labreport@berryoncology.com','suxiaoxing001@berryoncology.com']
carboncopy = ';'.join(cc)
#创建一个带附件的实例
message = MIMEMultipart()
message['From'] = Header('xingh3223@berryoncology.com','utf-8')
message['To'] = Header('xingh3223@berryoncology.com','utf-8')
message['Date'] = formatdate(localtime=True)
subject = '【BZ666样本】-%s'% sample
message['Subject'] = Header(subject,'utf-8')

#邮件正文内容
message.attach(MIMEText('Hi all,\n\n附件是北肿样本%s的北肿QC结果, 请查收。\n%s\n%s\n\n祝好,\n邢鹤'%(sample,carboncopy,message_from_check),'plain','utf-8'))

for name,filepath in upload_files.items():
	part1 = MIMEBase('application','octet-stream')
	part1.set_payload(open(filepath,'rb').read())
	encoders.encode_base64(part1)
	part1.add_header('Content-Disposition','attachment; filename="%s"'%(name))
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


