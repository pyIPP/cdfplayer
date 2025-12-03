import re
import numpy as np

def wr_for(arr_in, fmt='%13.6E', n_lin=6):
    arr_flat = arr_in.T.ravel()
    lines = [
        ''.join(fmt % val for val in arr_flat[i:i + n_lin])
        for i in range(0, len(arr_flat), n_lin)
    ]
    return '\n'.join(lines) + '\n'

def ssplit(line):
# This pattern matches a '-' that is immediately preceded by a digit
    line_clean = re.sub(r'(?<=\d)-', ' -', line)
    b = [float(num) for num in line_clean.split()]
    return b

def lines2fltarr(lines, dtyp=np.float32):
    return np.fromiter(
        (num for line in lines for num in ssplit(line)),
        dtype=dtyp)

def fltarr_len(lines, nx, dtyp=np.float32):
# for rabbitview
    data = []
    for jlin, line in enumerate(lines):
        data.extend(ssplit(line))
        if len(data) >= nx:
            break
    return jlin+1, np.array(data, dtype=dtyp)
