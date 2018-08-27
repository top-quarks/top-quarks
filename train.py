from os import system

system("g++ src/main.cpp -g -std=c++11 -O3 -DTRAIN -o train")
system("./train 1002")

#Please ignore all the ugly output

for i in range(100):
    cmd = "python3 src/train_pair_logreg.py %d"%(i)
    print(cmd)
    system(cmd)
