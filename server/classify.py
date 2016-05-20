#! /usr/bin/env python

import argparse
import h5py
from main import events_numpy_to_dict
import os
import requests
import json
import timeit
import glob 

parser = argparse.ArgumentParser(description='Classify a fast5 file')
parser.add_argument('--files', type=str, nargs='+',help='fast5 files')
parser.add_argument('--dir', type=str,help='fast5 files')
parser.add_argument('--host', type=str,help='fast5 files', default = "http://localhost:8001") #"http://40.85.98.117/"
args = parser.parse_args()

if args.files:
	filelist = args.files
else:
	filelist = glob.glob(args.dir + "/*")
for f in filelist:
	t0 = timeit.default_timer()
	try:
		fast5 = h5py.File(f,'r')
	except IOError:
		pass
	else:
		try:
			read=fast5["Raw"]["Reads"].keys()[0]
			try:
				events= fast5["Analyses"]["EventDetection_000"]["Reads"][str(read)]["Events"]
			except KeyError:
				events = fast5["Analyses"]["Basecall_RNN_1D_000"]["BaseCalled_template"]["Events"]
			d =  events_numpy_to_dict(events[100:400])
			d["id"] = os.path.basename(f)
			t1 = timeit.default_timer()	
			# print t1-t0
			response = requests.post(args.host, json=json.dumps(d))
			t2 = timeit.default_timer()
			# print t2-t0
			try:
				print response.json()
			except:
				print {"id" : d["id"], "error" : response.content}
		except KeyError:
			pass
			print {"id" : d["id"], "error" : "No Analyses Folder"}
			# print fast5.keys()
			if "Analyses" in fast5:
				print fast5["Analyses"]["Basecall_RNN_1D_000"]["BaseCalled_template"]["Events"].keys()
				print str(fast5["Analyses"]["Basecall_RNN_1D_000"]["BaseCalled_template"]["Fastq"])
