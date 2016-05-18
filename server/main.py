import numpy as np
from flask import Flask
app = Flask(__name__)

def events_numpy_to_dict(events):
	d = {}
	d["start"] = events["start"].tolist()
	d["length"] = events["length"].tolist()
	d["mean"] = events["mean"].tolist()
	d["stdv"] = events["stdv"].tolist()
	return d


def events_dict_to_numpy(d):
	events = np.empty(num_events, dtype=[('start', int), ('length', int),
                                            ('mean', float), ('stdv', float)])
	events["start"] = np.array(d["start"], dtype = float)
	events["length"] = np.array(d["length"], dtype = float)
	events["mean"] = np.array(d["mean"], dtype = float)
	events["stdv"] = np.array(d["stdv"], dtype = float)
	return events




@app.route("/")
def hello():
    return "Hello World!"

if __name__ == "__main__":
    app.run()
