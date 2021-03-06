# Running ReadUntil in the JR Hospital

It has to be ran in level 7 Office 7724 (as you can see in the picture. It cannot run in the MiSeq room as before as that room does have only one ethernet port for the NHS network and it is used by the MiSeq machine. 

So office 7724 is suitable to run the MinION. Lou Pankhourst can help you on settig in up, but as you can see in the pictures you should be able to unplug one of the computers network cable and connect it into the MinION computer.

We need to be in the NHS network so we can communicate to the KrakenServer that will be running the nanonet server.

[Room7724](pics/11.jpg)
[Connection](pics/22.jpg)
[Connection](pics/33.jpg)

I recommend using MinION 2 as it is updated and it is set up properly (software installed putty configured etc)

## Starting Nanonet Server
Open putty and you will see already a KrakenServer setup. If not, just log into:
nick@10.134.24.7
password: compass

Once you are inside:
```bash
cd readuntil_server/server
# I like to execute it inside a screen
uwsgi --buffer-size=32768 --processes 16 --http :8001 --http-timeout 3000 -w wsgi:app
```

## Running readuntil

After libprep, Lou will start minknow and for the firt 6 minutes it will start calibrating, at this stage we still don't run any readUntil, once that has finished, open a cmd.exe and execute Readuntil

```bash
cd Desktop\RUscrips_Nanonet_master 

c:\grouper\binaries\ont-python\python.exe gReadUntil_Nanonet.py -ip 127.0.0.1 -p 8001 -ip2 10.134.24.7 -p2 8001 -procs 16 -t 20 -skip
```

We are using -p 8001, however it can change, in order to check the port execute
```bash
c:\grouper\binaries\bin\mk_manager_client.exe --list

```

The output will be a JSON, you are interested in the value: ws_event_sampler_port. That is the value that should be in -p PORT

## Transfering data into GEL

Tha quickest way to transfer data is the following:

* Ypu need to switch to university network to be able to reach GEL or WT
* Use Lou credentials or you out private key
* Execute cmd.exe
```bash
c:\cygwin\cygwin.bat
cd /cygdrive/c/data/reads
ssh -p989 -i PATH_TO_YOUR_PRIVATE_KEY compass@83.151.222.218 'mkdir /mnt/microbio/HOMES/Nanopore/RUNNAME'
tar -c . | ssh -p989 compass@83.151.222.218 'cat > /mnt/microbio/HOMES/Nanopore/RUNNAME/reads.tar' 
```

You could copy directly to WT but it is good to have a copy in JR.

## Running simulations

In order to run simulations you need to run fast5index.py and provide a folder with fast5 files, it can have subfolders it will recurse on them. Then once the index is generated, you can use the client. I recoment to not to use speeding up fature as it's a bit buggy.

You can find error while processing reads with fast5index, and that can be caused because paths within fast5 files can be changed:

Right now in the python code there are two configurations:

```python
evs=fast5['Analyses']['EventDetection_000']['Reads'][read]['Events']
#evs=fast5['Analyses']['Basecall_RNN_1D_000']['BaseCalled_template']['Events']
```


Execution example:
```bash
python fast5index.py path/to/filder/containing/fast5s /path/to/index/file.idx
# run nanonet server on http://localhost:8001/
python client.py -i /path/to/index/file.idx -t 100 -r 200 -s http://localhost:8001/ -o /path/to/index/file.simulation.out

# -t is the trimming length
# -r is the read length that is being basecalled (number ov event sent to the server)
```
