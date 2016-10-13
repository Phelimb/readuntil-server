import h5py
import os
import glob
import sys
import array
import struct

if len(sys.argv)!=3:
	print "Wrong number of arguments"
	print sys.argv[0],"[PATH_TO_FAST5_FOLDER] [IndexName]"
	sys.exit(-1)

if not os.path.isdir(sys.argv[1]):
	print sys.argv[0],"[PATH_TO_FAST5_FOLDER] [IndexName]"
	sys.exit(-1)

if os.path.isfile(sys.argv[2]):
	if raw_input("Index already exists, overwrite? (N/y)").strip().lower()!='y': sys.exit(0)

readsummary=[]
mean_ary=array.array('d')
stdv_ary=array.array('d')
start_ary=array.array('l')
length_ary=array.array('l')
totevs=0


files = glob.glob(sys.argv[1]+"/*.fast5" )
print "Processing",len(files),"files"
for i,f in enumerate(files):
	if not i%1000:
		sys.stderr.write("\r{0}".format(i))
		sys.stderr.flush()
	try:
		fast5 = h5py.File(f,'r')
		read=fast5["Raw"]["Reads"].keys()[0]
		start = fast5["Raw"]["Reads"][str(read)].attrs.get("start_time")
		duration = fast5["Raw"]["Reads"][str(read)].attrs.get("duration")
		evs=fast5['Analyses']['EventDetection_000']['Reads'][read]['Events']
#		evs=fast5['Analyses']['Basecall_RNN_1D_000']['BaseCalled_template']['Events']
	
		
		mean_ary.extend(evs['mean'])
		stdv_ary.extend(evs['stdv'])
		if evs['start'][0].__class__.__name__=="float64": start_ary.extend((evs['start']*4000).astype(int))
		else: start_ary.extend(evs['start'])
		if evs['length'][0].__class__.__name__=="float64": length_ary.extend((evs['length']*4000).astype(int))
		else: length_ary.extend(evs['length'])

		readsummary.append([start,duration,totevs,len(evs),os.path.basename(f)])
		totevs+=len(evs)
	
		fast5.close()			
	except KeyError as e: 
		sys.stderr.write("\r{0} Error in file: {1}".format(i,f.strip()))
		sys.stderr.flush()


readsummary.sort(key=lambda x: x[0])
readsummary=zip(*readsummary)

output=open(sys.argv[2],'w')
output.write(struct.pack('l',len(readsummary[0])))
output.write(array.array('l',readsummary[0]).tostring())
output.write(array.array('l',readsummary[1]).tostring())
output.write(array.array('l',readsummary[2]).tostring())
output.write(array.array('l',readsummary[3]).tostring())
output.write(struct.pack('l',len(mean_ary)))
output.write(mean_ary.tostring())
output.write(stdv_ary.tostring())
output.write(start_ary.tostring())
output.write(length_ary.tostring())
output.write("\n".join(readsummary[4]))
output.close()


