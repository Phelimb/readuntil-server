import numpy as np
import flask
from flask import Flask, request
import json
from nanonet.features import events_to_features
from nanonet.segment import segment
from nanonet import decoding, nn
from nanonet.util import kmers_to_sequence
import tempfile
import pkg_resources
import subprocess
import timeit 
import os

app = Flask(__name__)

def events_numpy_to_dict(events):
    d = {}
    try:
        d["start"] = events["start"].tolist()
        d["length"] = events["length"].tolist()
        d["mean"] = events["mean"].tolist()
        d["stdv"] = events["stdv"].tolist()
    except ValueError:
        d["start"] = events["start"].tolist()
        d["length"] = events["length"].tolist()
        d["mean"] = events["mean"].tolist()
        d["stdv"] = events["variance"].tolist()        
    return d


def events_dict_to_numpy(d):
    events = np.empty(len(d["start"]), dtype=[('start', float), ('length', float),
                                            ('mean', float), ('stdv', float)])
    events["start"] = np.array(d["start"], dtype = float)
    events["length"] = np.array(d["length"], dtype = float)
    events["mean"] = np.array(d["mean"], dtype = float)
    events["stdv"] = np.array(d["stdv"], dtype = float)
    return events



__ETA__ = 1e-300
def _basecall(events,_id,  min_prob=1e-5, trans = None, trim=10):
    modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/default_template.npy'))
    network = np.load(modelfile).item()
    features =  events_to_features(events, window=[-1, 0, 1])
    features = features[trim:-trim]
    post = network.run(features.astype(nn.dtype))
    kmers = network.meta['kmers']
    # Do we have an XXX kmer? Strip out events where XXX most likely,
    #    and XXX states entirely
    if kmers[-1] == 'X'*len(kmers[-1]):
        bad_kmer = post.shape[1] - 1
        max_call = np.argmax(post, axis=1)
        good_events = (max_call != bad_kmer)
        post = post[good_events]
        post = post[:, :-1]
        if len(post) == 0:
            return None

    weights = np.sum(post, axis=1).reshape((-1,1))
    post /= weights

    post = min_prob + (1.0 - min_prob) * post
    trans = decoding.estimate_transitions(post, trans=trans)
    score, states = decoding.decode_profile(post, trans=np.log(__ETA__ + trans), log=False)

    # Form basecall
    kmer_path = [kmers[i] for i in states]
    seq = kmers_to_sequence(kmer_path)

    return seq

def run_bwa_mem(f):
    ref = os.path.join(os.path.dirname(__file__), 'data/NC_000962.3.fasta')
    l = ['bwa','mem','-v','1', '-x', 'ont2d',ref , f]
    outsam = subprocess.check_output(l)
    return outsam

def is_tb(sam):
	sam=sam.split("\n")
	return sam[2] and not int(sam[2].split('\t')[1]) == 4

@app.route('/', methods=['GET','POST'])
def process_events():
    if request.method == 'POST':
        t0 = timeit.default_timer()
        if isinstance(request.json,dict):
            data = request.json 
        else:
            data = json.loads(request.json)
            while not  isinstance(data,dict):
               data = json.loads(data)
        _id= data.get("id", "")
        events = events_dict_to_numpy(data)
        events, _ = segment(events, section='template')   
        t1 = timeit.default_timer()          
        # print "time to convert events", t1-t0
        seq = _basecall(events, _id)
        # print (_id, seq)
        t1a = timeit.default_timer()                  
        # print "time to basecall", t1a-t1    
        channel,tmpf = tempfile.mkstemp()

        os.write(channel,">%s\n" % _id)
        os.write(channel,"%s\n" % seq)   
        ff=open("/tmp/BCALLS","a")
        ff.write(">{0}\n{1}\n".format (_id,seq))
        ff.close()
        os.close(channel)

        t1b = timeit.default_timer()  
        # print "time to mktmp", t1b-t1a
        outsam = run_bwa_mem(tmpf)
        os.unlink(tmpf)
        t2 = timeit.default_timer()
        return flask.jsonify({"id" : _id, "is_tb" : is_tb(outsam), "response_time" : t2-t0})
    else:
        sdata = os.path.join(os.path.dirname(__file__), 'data/sample_events.json')
        with open(sdata,"r") as infile:
            data = json.load(infile)         
        return """<p>Please POST event data. e.g. </p><p> 

                curl -H "Content-Type: application/json" -X POST -d '%s' http://localhost:5000/

                </p>

                 """ % (str(data))


if __name__ == "__main__":
    app.run(host='0.0.0.0')
