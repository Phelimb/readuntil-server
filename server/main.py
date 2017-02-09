import numpy as np
import flask
from flask import Flask, request
import json
from nanonet.features import events_to_features
from nanonet.segment import segment
from nanonet import decoding, nn, cl
from nanonet.util import kmers_to_sequence,kmer_overlap
import tempfile
import pkg_resources
import subprocess
import timeit 
import os

app = Flask(__name__)


def clean_post(post, kmers, min_prob):
    # Do we have an XXX kmer? Strip out events where XXX most likely,
    #    and XXX states entirely
    if kmers[-1] == 'X'*len(kmers[-1]):
        bad_kmer = post.shape[1] - 1
        max_call = np.argmax(post, axis=1)
        good_events = (max_call != bad_kmer)
        post = post[good_events]
        post = post[:, :-1]
        if len(post) == 0:
            return None, None

    weights = np.sum(post, axis=1).reshape((-1,1))
    post /= weights
    post = min_prob + (1.0 - min_prob) * post
    return post, good_events


class ProcessAttr(object):
    def __init__(self, use_opencl=False, vendor=None, device_id=0):
        self.use_opencl = use_opencl
        self.vendor = vendor
        self.device_id = device_id


def make_feature_list(eventl,trim=10):
        #for e in eventl:
        X = events_to_features(eventl, window=[-1, 0, 1])
        X = X[trim:-trim]
        return X.astype(nn.dtype)

def _contribute(scores, qdata, posmap, n_bases):
     kmerlen = len(posmap)
     for kmer_pos, seq_pos in enumerate(posmap):
         index = (kmerlen - kmer_pos - 1) * n_bases
         scores[seq_pos, :] += qdata[index:index + n_bases]
     return
def form_basecall(qdata, kmers, states, qscore_correction=None):
    bases = sorted(set(''.join(kmers)))
    kmer_len = len(kmers[0])
    n_events = len(qdata)
    n_bases = len(bases)
    kmer_path = [kmers[i] for i in states]

    moves = kmer_overlap(kmer_path)
    seq_len = np.sum(moves) + kmer_len
    scores = np.zeros((seq_len, len(bases)), dtype=np.float32)
    sequence = list(kmer_path[0])
    posmap = range(kmer_len)

    _contribute(scores, qdata[0, :], posmap, n_bases)
    for event, move in enumerate(moves):
        if move > 0:
            if move == kmer_len:
                posmap = range(posmap[-1] + 1, posmap[-1] + 1 + kmer_len)
            else:
                posmap[:-move] = posmap[move:]
                posmap[-move:] = range(posmap[-move-1] + 1, posmap[-move-1] + 1 + move)
            sequence.append(kmer_path[event][-move:])
        _contribute(scores, qdata[event, :], posmap, n_bases)
    sequence = ''.join(sequence)
    base_to_pos = {b:i for i,b in enumerate(bases)}
    scores += __ETA__
    scoresums = np.sum(scores, axis=1)
    scores /= scoresums[:, None]
    called_probs = np.fromiter(
        (scores[n, base_to_pos[base]] for n, base in enumerate(sequence)),
        dtype=float, count=len(sequence)
    )

    if qscore_correction == 'template':
        # initial scores fit to empirically observed probabilities
        #   per score using: Perror = a.10^(-bQ/10). There's a change
        #   in behaviour at Q10 so we fit two lines. (Q10 seems suspicious).
        switch_q = 10
        a, b = 0.05524, 0.70268
        c, d = 0.20938, 1.00776
        switch_p = 1.0 - np.power(10.0, - 0.1 * switch_q)
        scores = np.empty_like(called_probs)
        for x, y, indices in zip((a,c), (b,d), (called_probs < switch_p, called_probs >= switch_p)):
            scores[indices] = -(10.0 / np.log(10.0)) * (y*np.log1p(-called_probs[indices]) + np.log(x))
    elif qscore_correction in ('2d','complement'):
        # same fitting as above
        if qscore_correction == 'complement':
            x, y = 0.13120, 0.88952
        else:
            x, y = 0.02657, 0.65590
        scores = -(10.0 / np.log(10.0)) * (y*np.log1p(-called_probs) + np.log(x))
    else:
        scores = -10.0 * np.log1p(-called_probs) / np.log(10.0)

    offset = 33
    scores = np.clip(np.rint(scores, scores).astype(int) + offset, offset, 126)
    qstring = ''.join(chr(x) for x in scores)
    return sequence, qstring, kmer_path


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

def _basecall(events,_id,  min_prob=1e-5, trans = None, trim=10,cORg='CPU'):
    #modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/default_template.npy'))
    #modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/r9.4_template.npy'))
    modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/r9_template.npy'))
    network = np.load(modelfile).item()

    if cORg=='CPU':
        features =  events_to_features(events, window=[-1, 0, 1])
        features = features[trim:-trim]
        kmers = network.meta['kmers']
        post = network.run(features.astype(nn.dtype))
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
        kmer_path = [kmers[i] for i in states]
        seq = kmers_to_sequence(kmer_path)

        return seq
    if cORg=='gpu':
        print cORg
        # get features
        features_list = [make_feature_list(events,trim=trim)]

        # set up gpu/opencl
        vendor='NVIDIA'
        #vendor='APPLE' 
        device_id=0
        n_files=1
        pa = ProcessAttr(use_opencl=True, vendor=vendor, device_id=int(device_id))
        platform = [
            p for p in cl.get_platforms()
            if p.get_devices(device_type=cl.device_type.GPU)
            and pa.vendor.lower() in p.get_info(cl.platform_info.NAME).lower()
        ][0]

        device = platform.get_devices(
            device_type=cl.device_type.GPU
        )[pa.device_id]
        ctx = cl.Context([device])

        # run networks
        max_workgroup_size = device.get_info(cl.device_info.MAX_WORK_GROUP_SIZE)
        queue_list = [cl.CommandQueue(ctx)] * n_files
        post_list = network.run(features_list, ctx, queue_list)
        post_list, good_events_list = zip(*(
            clean_post(post, network.meta['kmers'], min_prob) for post in post_list
        ))
        trans_list = [np.log(__ETA__ +
            decoding.fast_estimate_transitions(post, trans=trans))
            for post in post_list]
        score, states_list = decoding.decode_profile_opencl(ctx, queue_list, post_list, trans_list=trans_list,log=False)

        # form basecall
        kmers = network.meta['kmers']
        seq_list = []
        qual_list = []
        kmer_path_list = []
        for states, post in zip(states_list, post_list):
            seq, qual, kmer_path = form_basecall(post, [x for x in network.meta['kmers'] if 'X' not in x], states)
            #seq_list.append(seq)
            #kmer_path_list.append(kmer_path)
            #qual_list.append(qual)

        return seq


#def _basecall(events,_id,  min_prob=1e-5, trans = None, trim=10):
#    #modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/default_template.npy'))
#    modelfile = os.path.abspath(pkg_resources.resource_filename('nanonet', 'data/r9_template.npy'))
#    network = np.load(modelfile).item()
#    features =  events_to_features(events, window=[-1, 0, 1])
#    features = features[trim:-trim]
#    post = network.run(features.astype(nn.dtype))
#    kmers = network.meta['kmers']
#    # Do we have an XXX kmer? Strip out events where XXX most likely,
#    #    and XXX states entirely
#    if kmers[-1] == 'X'*len(kmers[-1]):
#        bad_kmer = post.shape[1] - 1
#        max_call = np.argmax(post, axis=1)
#        good_events = (max_call != bad_kmer)
#        post = post[good_events]
#        post = post[:, :-1]
#        if len(post) == 0:
#            return None
#
#    weights = np.sum(post, axis=1).reshape((-1,1))
#    post /= weights
#
#    post = min_prob + (1.0 - min_prob) * post
#    trans = decoding.estimate_transitions(post, trans=trans)
#    score, states = decoding.decode_profile(post, trans=np.log(__ETA__ + trans), log=False)
#
#    # Form basecall
#    kmer_path = [kmers[i] for i in states]
#    seq = kmers_to_sequence(kmer_path)
#
#    return seq

def run_bwa_mem(f):
    #ref = os.path.join(os.path.dirname(__file__), 'data/NC_000962.3.fasta')
    ref = os.path.join(os.path.dirname(__file__), 'hs_ref_GRCh38.fa')
    l = ['bwa','mem','-r','6','-v','1', '-x', 'ont2d',ref , f]
    outsam = subprocess.check_output(l)
    return outsam

def is_tb(sam):
	sam=sam.split("\n")
	#return sam[2] and not int(sam[2].split('\t')[1]) == 4
        return sam[733] and not int(sam[733].split('\t')[1]) != 4

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
        seq = _basecall(events, _id,cORg='gpu')
        #seq = _basecall(events, _id)
        print (_id, seq)
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
