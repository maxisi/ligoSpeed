#!/usr/bin/env python

from analysis import process
from analysis import results

import numpy as np

# PARAMETERS TO SET:
nf = 1e4
ninj = 100

detector = 'H1'
crab = 'J0534+2200'

injection_kinds = ['GR', 'G4v']
search_methods = ['GR', 'G4v']
pd = [0.0, np.pi/2, -np.pi/2]
plots = ['hinjrec', 'hinjs', 'hinjlins']

# PROCESS
for kind in injection_kinds:
    injsrch = process.InjSearch(detector, nf, kind, ninj)
    for p in pd:
        _, _, _ = injsrch.analyze(search_methods, injpdif=p)
        rh = results.Open('H1', kind, p)
        for plot in plots:
            rh.plots('J0534+2200', plot, extra_name='S5')