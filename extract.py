import sys
import os
import time

en = []

i = 0
while i < 10:
  
    os.system("python main.py")

    time.sleep(1)
    
    a = open("quantum_energy_parts_data/double_well_quantum_energy.dat","r+")

    count = 0
    for line in a:
        count += 1
        b = line.strip().split()
        if count == 1:
            en.append( float( b[3]) )

    a.close() 

    i += 1

for j in range(len(en)):
    print 2000, en[j]


"""
print "Max = ", max(en)
print "Min = ", min(en)
print "Range = ", [max(en),  min(en)]
print "Average = ", sum(en)/len(en)
"""

