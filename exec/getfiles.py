#!/usr/bin/python

import os,sys
mss = os.environ['MSS']
outpath = os.environ['VOL']+'/data'
imax = 0
run = ''
action = 'request'

args = sys.argv
if len(args)<2: sys.exit('no arguments')
for i,a in enumerate(args):
    if a in ['-r','-run']: run = args[i+1].zfill(6)
    elif a in ['-m','-max']: imax = int(args[i+1])
    elif a in ['-f','-from']: mss = args[i+1]
    elif a in ['-t','-to']: outpath = args[i+1]
    elif a in ['-c','-cancel']: action = 'cancel'
    elif a in ['-v','-view']: action = 'view'

# request files
if action=='request':
    l = [mss+'/'+x for x in os.listdir(mss) if run in x]
    l.sort()
    l.sort(key=len)
    if imax!=0: l = l[:imax+1]

    for x in l:
        print 'copy file {0} to {1}'.format(x,outpath)
        os.system('jget -n {0} {1}'.format(x,outpath))

# view queue
elif action=='view':
    if run=='': os.system('jqueue user mlevilla request jobState submit stub' )
    else: os.system('jqueue user mlevilla request jobState submit stub |  grep '+run)

# cancel jobs
elif action=='cancel':
    f= os.popen('jqueue user mlevilla request stub' )
    l = [x.split(' ') for x in f.readlines()[1:] if run in x]
    for x in l: 
        while x[0]=='': del x[0]
        os.system('jcancel '+x[0])
sys.exit(0)
