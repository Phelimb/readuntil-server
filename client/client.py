import struct
import array
import time
import threading
import json
import urllib2
import argparse

class NanoClient:
	def __init__(self,indexfile='index',trimlength=100,sampEvSize=300,server="http://localhost:8001/"):
		print "Loading index..."
		a=open(indexfile)
		sumlen=int(struct.unpack('l',a.read(8))[0])
		### summary ###
		self.start=array.array('l',a.read(8*sumlen))
		self.duration=array.array('l',a.read(8*sumlen))
		self.pos=array.array('l',a.read(8*sumlen))
		self.nevs=array.array('l',a.read(8*sumlen))
		###############
		evlen=int(struct.unpack('l',a.read(8))[0])
		### events ###
		self.mean_ary=array.array('d',a.read(8*evlen))
		self.stdv_ary=array.array('d',a.read(8*evlen))
		self.start_ary=array.array('l',a.read(8*evlen))
		self.length_ary=array.array('l',a.read(8*evlen))
		##############
		assert len(self.mean_ary)== len(self.stdv_ary)== len(self.start_ary)== len(self.length_ary) and not a.read()
		print "Index loaded"

		self.threads=[]
		self.trimlength=trimlength
		self.sampEvSize=sampEvSize
		self.server=server

	def simulateread(self,n):
		# Time when the read starts and should have finished
		startRead=(self.start[n])/float(4000)
		endRead=(self.start[n]+self.duration[n])/float(4000)

		# Current time
		curTime=time.time()-self.simStartTime

		# Pos of the first event in the event arrays
		pos=self.pos[n]
		# If we don't have enought events to trim, we ignore the read
		if self.nevs[n]<self.trimlength: return

		# if we don't have enough events to create a sample read, we calculate the end, to crceate a shorter read
		endSampEvSize=min(self.nevs[n],self.trimlength+self.sampEvSize)

		# We calculate the time at which the set of events will have to be sent to the server
		sendTime=(sum(self.length_ary[pos:pos+endSampEvSize])+self.start[n])/float(4000)

		# Json containing the data of the events
		data=json.dumps({'start':list(self.start_ary[pos+self.trimlength:pos+endSampEvSize]),
						 'length':list(self.length_ary[pos+self.trimlength:pos+endSampEvSize]),
						 'mean':list(self.mean_ary[pos+self.trimlength:pos+endSampEvSize]),
						 'stdv':list(self.stdv_ary[pos+self.trimlength:pos+endSampEvSize]),
						 'id':'-'
						 })

		while curTime<endRead:
			time.sleep(.1)
			curTime=time.time()-self.simStartTime
			if curTime>sendTime:
				try:
					req = urllib2.Request(self.server, data, {'Content-Type': 'application/json', 'Content-Length': len(data)})
					f = urllib2.urlopen(req)
					response = json.loads(f.read())
					curTime=time.time()-self.simStartTime
					print float(curTime-startRead)*100/(endRead-startRead),self.nevs[n]
					if response["is_tb"]: print "IS TB!!"
				except: print "Error!"
				break


	def cleanThreads(self):
		done=[i for i in self.threads if not i.isAlive()]
		for i in done:
			i.join()
		self.threads=[i for i in self.threads if i not in done]

	def simulate(self,nreads=0):
		self.simStartTime=time.time()-self.start[0]/float(4000)
		lastThreadClean=0
		

		readn=0
		if nreads: totreads=nreads
		else: totreads=len(self.start)

		while readn<totreads:
			curTime=time.time()-self.simStartTime
			while self.start[readn]/float(4000)<curTime:
				th=threading.Thread(target=self.simulateread,kwargs={'n':readn})
				th.start()
				self.threads.append(th)
				readn+=1

			if curTime-lastThreadClean>5:
				print len(self.threads)
				self.cleanThreads()
				print len(self.threads),"\n"
				
				lastThreadClean=curTime

			time.sleep(.1)

		for i in self.threads:
			i.join()
		print  len(self.threads)


if __name__=='__main__':
	parser = argparse.ArgumentParser(description='client')
	parser.add_argument("-i",dest="index",help="Fast5 Index file",required=True)
	parser.add_argument("-t",dest="trimlen",help="Event triming length [100]",default=100, type=int)
	parser.add_argument("-r",dest="readlen",help="Read length (measured in events) [300]",default=300, type=int)
	parser.add_argument("-n",dest="nreads",help="#reads to simulate [0=all]",default=0, type=int)
	parser.add_argument("-s",dest="server",help="Http Server url [http://localhost:8001/]",default="http://localhost:8001/")
	args = parser.parse_args()

	nc=NanoClient(args.index,args.trimlen,args.readlen,args.server)
	nc.simulate(args.nreads)
