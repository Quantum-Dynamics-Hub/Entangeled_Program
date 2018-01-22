import sys
import os
import math
import time

en = []

i = 0
while i < 100:
  
    os.system("python main.py")

    time.sleep(5)
    
    a = open("energy/energy.txt","r+")

    count = 0
    for line in a:
        count += 1
        b = line.strip().split()
        if count == 1:
            en.append( float( b[3]) )

    a.close() 

    i += 1

avg = sum(en)/len(en)
print "\n"
print "##### Results from extract #####"
print "Max = ", max(en)
print "Min = ", min(en)
print "Range = ", [max(en),  min(en)]
print "Average = ", sum(en)/len(en)
#print "Percent of Cubic Barrier = ", ( (sum(en)/len(en))/0.015 ) * 100.0, " %"
print "Percent of Symmetric Double Well Barrier = ", ( ( (-0.25) - (sum(en)/len(en)) ) /(-0.25) ) * 100.0, " %"

summ = 0.0
for i in range(len(en)):
    summ = summ + (en[i] - avg)*(en[i] - avg)
std_dev = math.sqrt( summ/len(en) )

print "Standard Deviation = ", std_dev
 

