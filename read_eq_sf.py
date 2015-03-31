import sys, os
import Tkinter as tk
import numpy as np
import psi2rz
import dd_20140403 as dd

sf  = dd.shotfile()


def get_grid(runid):

# Parse namelist for info about equilibrium exp-diag-ed

    nshot = int(runid[:5])
    io_dir = '%s/tr_client/AUGD/%d/%s' %(os.getenv('HOME'), nshot, runid[5:])

    nml = '%s/%sTR.DAT' %(io_dir, runid)
    print('Parsing %s' %nml)

    fnml = open(nml, 'r')
    eq_d = {}
    for line in fnml.readlines():
        for eqstr in ('!exp', '!dia', '!ed'):
            if eqstr in line:
                eq_d[eqstr] = line.split('=')[1].strip()
    print(eq_d)

# Fetch the matching equilibrium shotfile

    if sf.Open(eq_d['!dia'], nshot, experiment=eq_d['!exp'], edition=int(eq_d['!ed'])):
        for sgr in ('PFM', 'Ri', 'Zj', 'PFxx'):
            eq_d[sgr] = sf.GetSignal(sgr)
        for sig in ('time', 'ixti'):
            eq_d[sig] = sf.GetSignal(sig)
        sf.Close()
    return eq_d
