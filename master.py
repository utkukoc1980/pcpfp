#!/usr/bin/env python3
import pickle
import socket
import sys
import random
import subprocess
import os
import shlex
import os.path
import glob
from time import sleep

def fresh_start(mydir,start):
	for fname in os.listdir(mydir):
		if fname.startswith(start):
			os.remove(os.path.join(mydir,fname))

PORT1=65531
PORT2=65532
HOST1='10.2.1.121'
HOST2='10.2.1.122'

fresh_start("/home/kocu/3501g/DV3","zzz")

subprocess.Popen(["ssh" ,"kocu@10.2.1.121", "cd 3501g/DV3;rm zzz*"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

subprocess.Popen(["ssh" ,"kocu@10.2.1.122", "cd 3501g/DV3;rm zzz*"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

command = "./pcpfp -I 0 -W 3";
subprocess.Popen(shlex.split(command), shell=False)

subprocess.Popen(["ssh" ,"kocu@10.2.1.121", "cd 3501g/DV3; touch zzz_SLAVE1_signal.in"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

subprocess.Popen(["ssh" ,"kocu@10.2.1.122", "cd 3501g/DV3;  touch zzz_SLAVE2_signal.in"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

subprocess.Popen(["ssh" ,"kocu@10.2.1.121", "cd 3501g/DV3; touch zzz_SLAVE1_object.in"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

subprocess.Popen(["ssh" ,"kocu@10.2.1.122", "cd 3501g/DV3;  touch zzz_SLAVE2_object.in"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

subprocess.Popen(["ssh", "kocu@10.2.1.121", "cd 3501g/DV3; ./pcpfp -I 1 -W 3"],
		shell=False,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE)

subprocess.Popen(["ssh", "kocu@10.2.1.122", "cd 3501g/DV3; ./pcpfp -I 2 -W 3"],
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

signals_to_slave1=[]
signals_to_slave2=[]
objects_to_slave1=[]
objects_to_slave2=[]

sleep(5)

j=0
while j==0:

	while os.path.exists("/home/kocu/3501g/DV3/zzz_SLAVE0_signal.out") is False:
		continue
	else:
		with open("/home/kocu/3501g/DV3/zzz_SLAVE0_signal.out","r") as file:
			signal_out=file.read()

		
	signal_out=signal_out.split("\n")
	signal_out.pop(-1)
	
	for i in range(0,len(signal_out)):
		signal_out[i]=signal_out[i].split(" ")

	
	while os.path.exists("/home/kocu/3501g/DV3/zzz_SLAVE0_object.out") is False:
		continue


	with open("/home/kocu/3501g/DV3/zzz_SLAVE0_object.out","r") as file:
		object_out=file.read()
	
	object_out=object_out.split("\n")
	object_out.pop(-1)	

	for i in range(0,len(object_out)):
        	object_out[i]=object_out[i].split(" ")

	
	for i in range(0,len(signal_out)):
		if signal_out[i][0]=='-111':
			if signals_to_slave1.count(signal_out[i])==0:
				signals_to_slave1.append(signal_out[i])
			if signals_to_slave2.count(signal_out[i])==0:
				signals_to_slave2.append(signal_out[i])
		elif signal_out[i][0]=='1':
			if signals_to_slave1.count(signal_out[i])==0:
				signals_to_slave1.append(signal_out[i])
		elif signal_out[i][0]=='2':
			if signals_to_slave2.count(signal_out[i])==0:
				signals_to_slave2.append(signal_out[i])

	for i in range(0,len(object_out)):
		if object_out[i][0]=='-111':
			if objects_to_slave1.count(object_out[i])==0:
				objects_to_slave1.append(object_out[i])
			if objects_to_slave2.count(object_out[i])==0:			
				objects_to_slave2.append(object_out[i])
		elif object_out[i][0]=='1':
			if objects_to_slave1.count(object_out[i])==0:
				objects_to_slave1.append(object_out[i])
		elif signal_out[i][0]=='2':
			if objects_to_slave2.count(object_out[i])==0:
				objects_to_slave2.append(object_out[i])			
	
	sock_1=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
	sock_1.connect((HOST1,PORT1))
	sock_1.send(pickle.dumps(signals_to_slave1+objects_to_slave1))
	sock_1.shutdown(socket.SHUT_WR)

	sock_2=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
	sock_2.connect((HOST2,PORT2))
	sock_2.send(pickle.dumps(signals_to_slave2+objects_to_slave2))
	sock_2.shutdown(socket.SHUT_WR)

	data1=pickle.loads(sock_1.recv(1024))
	data2=pickle.loads(sock_2.recv(1024))

	for i in range(0,len(data1)):
		data1[i]=data1[i][0].replace("0","1")

	for i in range(0,len(data2)):
		data2[i]=data2[i][0].replace("0","2")

	data=data1+data2

	for item in range(0,len(data)):
		if len(data[item])==3:
			with open("home/kocu/3501g/DV3/zzz_SLAVE0_signal.in","a") as file:
				file.write(str(data[item][0])+ " " + str(data[item][1]) + " " + str(data[item][2]) + "\n")       
		else:
			with open("home/kocu/3501g/DV3/zzz_SLAVE0_object.in","a") as file:
				for element in range(0,len(data[item])):
					file.write(str(data[item][element]) + " ")
				file.write("\n")

