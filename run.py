from os import system


for i in range(100):
    cmd = "python3 class22.py %d"%(i)
    print(cmd)
    system(cmd)
