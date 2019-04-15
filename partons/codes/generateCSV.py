#!/usr/bin/python

binsX =  100
binsT =  100
xv    =  0.1
tv    = -0.1

for x in range(0,200):
    xv =  xv + 0.004
    tv = -0.1
    for t in range(0,200):
        print str(xv)+"|"+str(tv)+"|2.0|2.0|2.0"
        tv =  tv - 0.02
