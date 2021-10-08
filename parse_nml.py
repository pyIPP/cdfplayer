import os
import numpy as np

char_list = (1794, 3842, 7938, 12034, 16130, 18178)


def parsenml(nml, parname, fmt=5):

    val = []

    for line in open(nml):
        if len(line.strip()) > 0:
            if line.strip()[0] != '!':
                tmp = line.split('=')
                vname = tmp[0].split('(')[0].strip()
                if parname.upper() == vname.upper():
                    valstr = tmp[1].split('!')[0].strip()

                    if fmt in char_list:
                        val.append( valstr.replace("'", "") )
                    if fmt in (3, 4):
                        for entry in valstr.split(','):
                            val.append(np.int32(entry))
                    if fmt == 7:
                        for entry in valstr.split(','):
                            if entry.upper().strip() in ['F', '.F.', 'FALSE', '.FALSE.']:
                                val.append( np.bool(False))
                            else:
                                val.append(np.bool(True))
                    if fmt == 5:
                        for entry in valstr.split(','):
                            val.append( np.float32(entry))

    if fmt in char_list:
        return val
    else:
        return np.array(val)


if __name__ == "__main__":
    nml = '/afs/ipp/home/g/git/tr_client/AUGD/23076/W01/23076W01TR.DAT'
    cl1 = parsenml(nml, 'gfrac', fmt=5)
    print(cl1)
    cl2 = parsenml(nml, 'frac', fmt=5)
    print(cl2)
    cl3 = parsenml(nml, 'nadvsim', fmt=3)
    print(cl3)
    cl4 = parsenml(nml, 'densim', fmt=5)
    print(cl4)
