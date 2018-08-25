from trackml.dataset import load_event, _load_event_data
from trackml.randomize import shuffle_hits
from trackml.score import score_event
import pandas
from sys import argv

if len(argv) < 2:
    print("Usage: %s <event id>"%argv[0])
    exit()
filenum = int(argv[1])
sub = pandas.read_csv('submissions/submission%d.csv'%filenum, header=0, index_col=False, dtype={'event_id':'i4','hit_id': 'i4','track_id':'i4'})

filenum = int(sub['event_id'][0])
hits, cells, particles, truth = load_event('/data/kaggle/train_100_events/event%09d'%filenum)

#sub = sub.drop(['event_id'], axis=1)
#print(sub)

print("Score: ", score_event(truth, sub))
