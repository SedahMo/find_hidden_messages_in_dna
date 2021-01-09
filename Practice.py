

import matplotlib.pyplot as plt
year = [1950, 1970, 1990, 2010]
pop = [2.519, 3.692, 5.263, 6.972]
plt.scatter(year, pop)
plt.show()

world = {'afg':30, 'ab':12}
print(world.items())

print(list(1:10))


import numpy as np
final_tails = []
for x in range(10000):
    tails = [0]
    for i in range(1000):
        coin = np.random.randint(0, 2)
        tails.append(tails[i] + coin)
    final_tails.append(tails[-1])


import matplotlib.pyplot as plt
plt.hist(final_tails, bins = 30)
plt.show()



import numpy as np
import matplotlib.pyplot as plt
a = open('/Users/babemomo/Dropbox/My Mac (natekiMacBook-Pro.local)/Downloads/dataset_7_6.txt')
b = a.read()
plt.plot(MinSkew(b.strip()))
plt.show()



temp = open('/Users/babemomo/Dropbox/My Mac (natekiMacBook-Pro.local)/Downloads/dataset_9_7 (1).txt').read().split('\n')
z = FrequentWordwithMismatch(temp[0], temp[1], temp[2])
print(*z)


ACT
def Neighborhood(text, d):
    
    for i in "ATCG":
        
        
        
data = 'From stephen.marquard@uct.ac.za Sat Jan  5 09:14:16 2008'
atpos = data.find('@')
print(atpos)
sppos = data.find(' ', atpos)
host = data[atpos + 1: sppos]
print(host)

text = "X-DSPAM-Confidence:    0.8475";
a = text.find(':')
b = text[a:].strip()
print(float(b))



name = input("Enter file:")
if len(name) < 1 : name = "mbox-short.txt"
handle = open(name)

hour_list = []
for line in handle:
    if line.startswith('From ') == True:
    	hour_list.append(line.strip().split()[5].split(':')[0])

hour_freq = {}
for hour in hour_list:
    hour_freq[hour] = hour_freq.get(hour, 0) + 1
    
for k, v in sorted([(k, v) for k, v in hour_freq.items()]):
    print(k, v)
    
def ImmediateNeighbors(pattern):
    neighborhood = [pattern]
    nuc = ['A', 'T', 'C', 'G']
    for i in range(len(pattern)):
        sym = pattern[i]
        for j in nuc:
            if sym != j:
                thislist = [pattern[:i], pattern[(i+1):]]
                neighborhood.append(j.join(thislist))
    return neighborhood

def IterativeNeighbors(pattern, d):
    neighbor_set = set([pattern])
    neighbor_dict = {}
    for i in range(d):
        count = iter(neighbor_set.copy())
        for j in count:
            neighbor_dict[j] = neighbor_dict.get(j, ImmediateNeighbors(j))
            neighbor_set.update(neighbor_dict[j])
    with open("print_temp.txt", "w") as f:
        print(*neighbor_set, sep = '\n', file = f)
    


for i in {1,2,3,4}:
    print(i)