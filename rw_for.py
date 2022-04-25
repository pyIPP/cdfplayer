import numpy as np

def wr_for(arr_in, fmt='%13.6E', n_lin=6):

    arr_flat = arr_in.T.ravel()
    nx = len(arr_flat)
    out_str=''
    for jx in range(nx):
        out_str += (fmt %arr_flat[jx])
        if (jx%n_lin == n_lin - 1):
            out_str += '\n'
# Line break after writing data, but no double break
    if (nx%n_lin != 0):
        out_str += '\n'

    return out_str


def ssplit(ll):

    tmp = ll.replace('-', ' -')
    tmp = tmp.replace('e -', 'e-')
    tmp = tmp.replace('E -', 'E-')
    slist = tmp.split()
    a = [float(i) for i in slist]

    return a


def lines2fltarr(lines, dtyp=None):

    data = []
    for line in lines:
        data += ssplit(line)
    if dtyp is None:
        dtyp = np.float32

    return np.array(data, dtype=dtyp)


def fltarr_len(lines, nx, dtyp=None):

    data = []
    for jlin, line in enumerate(lines):
        data += ssplit(line)
        if len(data) >= nx:
            break
    if dtyp is None:
        dtyp = np.float32

    return jlin+1, np.array(data, dtype=dtyp)
