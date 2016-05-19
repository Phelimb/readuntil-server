#! /usr/bin/env python

import argparse
import h5py
from main import events_numpy_to_dict
import os
import requests
import json

parser = argparse.ArgumentParser(description='Classify a fast5 file')
parser.add_argument('files', type=str, nargs='+',help='fast5 files')
args = parser.parse_args()

for f in args.files:
	fast5 = h5py.File(f,'r')
	read=fast5["Raw"]["Reads"].keys()[0]
	d =  events_numpy_to_dict(fast5["Analyses"]["EventDetection_000"]["Reads"][str(read)]["Events"])
	d["id"] = os.path.basename(f)
	response = requests.post("http://127.0.0.1:5000/", json=json.dumps(d))
	print response.json()
