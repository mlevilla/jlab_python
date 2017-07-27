#!/usr/bin/env python

from params import *

t = catch_arg('t','all')
run_min = catch_arg('min',889)
run_max = catch_arg('max',1516)
run0 = catch_arg('r',-1)
c = catch_arg('c','','end')

runs = [run0] if run0!=-1 else get_runs_between(run_min,run_max,t)
for run in runs: _ = os.system(c+' -b -r '+str(run))


