#!/usr/bin/python
#coding=utf-8
#xinghe <xingh3223@berryoncology.com>

from __future__ import division
import os
import subprocess
import glob
import time
import re
import argparse

my_env = os.environ.copy()
print(my_env)
#my_env["PATH"] = "/usr/sbin:/sbin:" + my_env["PATH"]
#subprocess.Popen(my_command, env=my_env)

def checked_list(path = None, shell_name = 'runbz1.sh'):
	'''
	待检查的job_id列表.
	'''
	checked_list_id = []
	if path == None:
		root_dir = os.getcwd()
	else:
		root_dir = os.path.abspath(path)
	time.sleep(30)	
	stderr_list = glob.glob(root_dir +'/'+ '%s.e*'%shell_name)
	flag = True
	total_wait_time = 30
	while flag == True:
		stderr_list = glob.glob(root_dir +'/'+ '%s.e*'%shell_name)
		if stderr_list == []:
			time.sleep(30)
			total_wait_time += 30
			#if total_wait_time >= 900:
			if total_wait_time >= 100:
				#raise Exception('超过15分钟在该路径下未检测到投递的job_ids')
				print('超过15分钟在该路径下未检测到投递的job_ids')
		else:
			flag = False

	for	i in stderr_list :
		i = os.path.basename(i).split('.')[-1].lstrip('e')
		checked_list_id.append(i)

	return checked_list_id


def check_job(checked_user_list,checked_list=[],root_dir='.'):
	'''
	检查jobs是否仍然在进行.
	'''
	root_dir = os.path.abspath('./')
	print('root_dir: '+root_dir)
	check = True
	wait_time = 0
	message = ''
	print('checked_user_list:qian')
	print(checked_user_list)

	while check == True:

		job_ids = []	
		for each_user in checked_user_list:
			obj = subprocess.Popen('/opt/gridengine/bin/linux-x64/qstat -u %s'%each_user,shell=True,universal_newlines=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=my_env)
			print('obj:')
			for line in obj.stdout:
				print('line:')
				print(line)
				line = line.strip()
				if line.startswith('---'):
					pass
				else:
					job_id = re.split(r'\s+',line)[0]
					print('111:'+job_id)
					job_id_status = re.split(r'\s+',line)[4]
					job_ids.append(job_id)
			for b in obj.stderr:
				print('err:')
				print(b)

		checked_list_flag = []

		print('all:')
		print(job_ids)
		print('checked_list:')
		print(checked_list)

		for each in checked_list:
			if each not in job_ids:
				each_flag = True
			else:
				each_flag = False
			checked_list_flag.append(each_flag)

		finish_flag = all(checked_list_flag)
		print(checked_list_flag)
		print('finish_flag: '+str(finish_flag))

		if finish_flag == True:
			check=False
			message = '任务已经结束,路径:%s'%root_dir
			return message
		else:
			time.sleep(1800) #1800 
			wait_time += 1800
			if wait_time >= 28800:
				message = '分析时间已经超过8小时,请查看,路径:%s'%root_dir
				return message
			continue	


if __name__ == '__main__':

	parser=argparse.ArgumentParser(prog="",description="检测qsub上去的任务是否完成")
	parser.add_argument("-u","--user-list",help='被检查的user列表,如:xinghe,csmartpro',required = True)
	parser.add_argument("-s","--shell-name",help="被投递上去的脚本",required = True)
	args=parser.parse_args()

	pipeline_dir = '/share/work1/xinghe/proc/BZQC/pipeline' 

	if args.user_list == None:
		checked_user_list = []
		local_user = subprocess.Popen('whoami',shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		local_user = local_user.stdout.read().strip()
		checked_user_list.extend([local_user,'csmartpro'])
	else:
		checked_user_list = args.user_list.strip().split(',') #xinghe,csmartpro
	
	print("checked_user_list: "+','.join(checked_user_list))

	if args.shell_name == None:
		checked_list_id = checked_list(shell_name = 'runbz1.sh')
	else:
		checked_list_id = checked_list(shell_name = args.shell_name)
		
	print(checked_list_id)
	message = check_job(checked_user_list = checked_user_list,checked_list=checked_list_id)
	print(message)
	root_dir = os.path.abspath('./')
	mail_tmp_file = open('%s/mail_tmp.sh'%root_dir,'r')
	content_tmp = mail_tmp_file.read().strip()
	mail_tmp_file.close()
	mail_file = open('%s/mail.sh'%root_dir,'w')
	mail_file.write(content_tmp+' -m '+message+'\n')
	mail_file.close()
	os.system('sh %s/mail.sh'%(root_dir))



