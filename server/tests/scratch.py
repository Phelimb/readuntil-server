import sys
sys.path.append("..")
from server.main import events_numpy_to_dict
import h5py
import json
fast5 = h5py.File("tests/minion1_PC_MN16255_FAA90892_mtub_L42182_200416_1717_1_ch244_read1222_strand.fast5",'r')
read=fast5["Analyses"]["EventDetection_000"]["Reads"].keys()[0]
events=fast5["Analyses"]["EventDetection_000"]["Reads"][read]["Events"][:100]

print json.dumps(events_numpy_to_dict(events))
fast5.close()