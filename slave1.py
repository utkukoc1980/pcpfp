#!/usr/bin/env python3

import pickle
import socket
import sys
import random
import shlex
import subprocess
from time import sleep
import os



input=sys.stdin.buffer.read()
slave1_input=pickle.loads(input)
slave1_input[0]="0"

while os.path.exists("/home/kocu/3501g/DV3/zzz_SLAVE1_signal.in") is False:
        continue

with open("home/kocu/3501g/DV3/zzz_SLAVE1_signal.in","a") as file:
        file.write(str(slave1_input[0])+ " " + str(slave1_input[1]) + " " + str(slave1_input[2]) + "\n")


while os.path.exists("/home/kocu/3501g/DV3/zzz_SLAVE1_signal.out") is False:
        continue

with open("/home/kocu/3501g/DV3/zzz_SLAVE1_signal.out","r") as file:
        slave1_signal_output=file.read()

slave1_signal_output=slave1_signal_output.split("\n")
slave1_signal_output.pop(-1)
sys.stdout.buffer.write(pickle.dumps(slave1_signal_output))

