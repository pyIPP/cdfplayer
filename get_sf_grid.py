import sys, os
sys.path.append('/afs/ipp/home/g/git/python/sfutils')

try:
    import Tkinter as tk
except:
    import tkinter as tk
import numpy as np
import sfread


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

    equ = sfread.SFREAD(nshot, eq_d['!dia'], exp=eq_d['!exp'], ed=int(eq_d['!ed']))
    print(equ.status)
    if equ.status:
        for obj in ('PFM', 'Ri', 'Zj', 'PFxx', 'time', 'ixti'):
            eq_d[obj] = equ.getobject(obj)

    return eq_d
