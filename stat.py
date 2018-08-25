from glob import glob

def nice(t):
    return str(int(t//60))+':'+str(t%60)

tot, c, ma = 0,0,0
mem_tot, mem_ma = 0,0
for i in glob("logs/*"):
    f = open(i, 'r')
    t = f.read()
    f.close()
    sa = 'Elapsed (wall clock) time (h:mm:ss or m:ss): '
    a = t[t.index(sa)+len(sa):]
    m,s = a[:a.index('\n')].split(':')
    elapsed = float(m)*60+float(s)

    sb = 'Maximum resident set size (kbytes): '
    b = t[t.index(sb)+len(sb):]
    mem = float(b[:b.index('\n')])
    mem_tot += mem
    mem_ma = max(mem_ma, mem)

    tot += elapsed
    c += 1
    ma = max(ma, elapsed)
    if (ma == elapsed):print(nice(elapsed), mem/2**20, i)

print(nice(tot/c), nice(ma), c)
print(mem_tot/c/2**20, mem_ma/2**20)
