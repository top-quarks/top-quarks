from subprocess import call
from concurrent.futures import ThreadPoolExecutor as Pool
import os
from sys import argv


os.system("g++ src/main.cpp -g -std=c++11 -O3 -DEVAL -o predict")

parallel = 4

if len(argv) == 1:
    print("Running on events 0 - 124")
    run_list = list(range(125)) #Test
else:
    print("Running on events 1000 - 1099")
    run_list = list(range(1000, 1100)) #Validate

n = len(run_list)

def run(i):
    log = 'logs/%s.txt'%str(i)
    retval = call(['bash', '-c', '/usr/bin/time -v ./predict %s > %s 2>&1'%(str(i), log)], stderr=open(log,'ab'))
    if retval != 0:
        os.system("cat "+log)
        print('If you see "couldn\'t open" somewhere in the above, there was a problem with finding the data. Make sure you set base_path correctly in src/input.hpp')
        exit(0)
    return 1

with Pool(max_workers=parallel) as executor:
    done = 0
    for _ in executor.map(run, run_list):
        done += 1
        print("%d / %d     \r"%(done, n), end = "")
print()

body = []
for i in run_list:
    f = open("submissions/submission%d.csv"%i, "r")
    l = f.read().split("\n")
    while not l[-1]:
        l = l[:-1]
    head = [l[0]]
    body += l[1:]
    f.close()
print("Read")
f = open("submissions/submission.csv", "w")
f.write("\n".join(head+body)+'\n')
f.close()
print("Wrote")

os.system("gzip -kf submissions/submission.csv")
