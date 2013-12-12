#!/usr/bin/python

import os

environment=os.environ

#environement.items()

for (k,v) in environment.items():
	print k,"=",v
