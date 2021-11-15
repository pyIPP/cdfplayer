""" 
Write/Read data in UFILES format (version 0.2)

UFILES is an ASCII format defined in the manual:
  http://w3.pppl.gov/~pshare/help/ufiles_manual.html

How to use ufiles.py:
  
  import ufiles as uf

======
Read:
======

  u = uf.RU(fname)

  Input:
    fname [string] is the filename with the full path 

  Output:
    u.shot: shot number as appear in the first line of the ufiles (int)
    u.dim: rank of the function, i.e. rnk(f(X)) = 1, rnk(F(X,Y))=2 
    u.label['X','Y','Z']: dict of the independent variables' labels
    u.values['X','Y','Z']: dict of numpy arrays for the X, Y and Z grid
    u.flabel: label of the dependent variable
    u.fvalues: numpy array of the function data

  Example:
    plt.plot(u.values['Y'],u.fvalues[0,:])

======
Write:
======

  u = uf.WU(u_dic,udir=path)

  Input:
    path [string] is the target dir; default is $HOME/udb/$shot
    u_dic is the input dictionary:
      u_dic['pre']  = output file prefix [string]
      u_dic['ext']  = output file extension [string]
      u_dic['shot'] = shot number [string]
      u_dic['scal'] [list]
        u_dic['scal']=[[label1 [string],scal_val1 [float]], [label2, scal_val2], ...]
      u_dic['grid'] [dictionary]
        u_dic['grid']['X'] [dictionary]
          u_dic['grid']['X']['lbl'] = grid variable label + unit [string]
          u_dic['grid']['X']['arr'] = grid variable array [numpy 1D array]
      Same for 'Y' or 'Z' in case of 2D, 3D ufiles
      u_dic['data'] [dictionary]
        u_dic['data']['lbl'] = independent variable label + unit [string]
        u_dic['data']['arr'] = independent variable array [numpy 1D, 2D or 3D array]


"""

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '0.02.1'
__date__    = '03.11.2016'

import datetime, os, sys, shutil, logging
import numpy as np
import rw_for

logger = logging.getLogger('trview.ufiles')

now = datetime.datetime.now()


def WU(uf_d, udir=None):

    nshot  = int(uf_d['shot'])
    s5shot = '%05d' %nshot
    prestr = uf_d['pre']
    extstr = uf_d['ext']
    comment = 'Generated with write_u.py'
    if 'comm' in uf_d.keys():
        comment = uf_d['comm']
    if udir is None:
        udir = '%s/udb/%s' %(os.getenv('HOME'), s5shot)
    os.system('mkdir -p %s' %udir)
    uf = '%s/%s%s.%s' %(udir, prestr, s5shot, extstr)
    if os.path.isfile(uf):
        ufbak = uf
        jext = 1
        while os.path.isfile(ufbak):
            ufbak = '%s.%d' %(uf, jext)
            jext += 1
        logger.info('mv %s %s' %(uf, ufbak))
        shutil.move(uf, ufbak)

    dev = 'AUGD'
    xyz = ('X', 'Y', 'Z')

# Dimensions check
    coords = []
    dims = {}
    for coord in xyz:
        if coord in uf_d['grid'].keys():
            coords.append(coord)
            dims[coord] = len(uf_d['grid'][coord]['arr'])
    ndim = len(coords)
    if ndim != uf_d['data']['arr'].ndim:
        logger.error('Data ndim inconsistent with n of independent variables')
        logger.error('No u-file written')
        return
    for jco, coord in enumerate(coords):
        if dims[coord] != uf_d['data']['arr'].shape[jco]:
            logger.error('The %dth index has inconsistent grid vs data dimension: %d, %d', \
               jco, dims[coord], uf_d['data']['arr'].shape[jco] )
            logger.error('No u-file written')
            return
        
# Header

    ufs = '  %5d%4s %1i 0 6              ;-SHOT #- F(X) DATA -UF%1iDWR- %s\n'\
         %(nshot,  dev,  ndim,  ndim, now.strftime('%d%b%Y'))
    ufs += ''.ljust(30) + ';-SHOT DATE-  UFILES ASCII FILE SYSTEM\n'

    if 'scal' in uf_d.keys():
        nscal = len(uf_d['scal'])
        ufs += (' %2d' %nscal).ljust(30) + \
            ';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'
        for jscal in range(nscal):
            lbl, val = uf_d['scal'][jscal]
            ufs += ('%11.4E' %val).ljust(30) + ';-SCALAR,  LABEL FOLLOWS:\n'
            ufs += lbl.ljust(30) + '\n'
    else:
        ufs += '  0'.ljust(30) +';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'

    if ndim > 0:
        if (ndim == 1):
            if 'unit' in uf_d['grid']['X'].keys():
                ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL-\n' \
                %(uf_d['grid']['X']['name'], uf_d['grid']['X']['unit'])
            else:
                ufs += '%-30s;-INDEPENDENT VARIABLE LABEL-\n' %uf_d['grid']['X']['lbl']
        else:
            for coord in coords:
                if 'unit' in uf_d['grid'][coord].keys():
                    ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL: ' \
                           %(uf_d['grid'][coord]['name'], uf_d['grid'][coord]['unit'])
                else:
                    ufs += '%-30s;-INDEPENDENT VARIABLE LABEL: ' %uf_d['grid'][coord]['lbl']
                ufs += '%1s-\n' %coord

        if 'unit' in uf_d['data'].keys():
            ufs += '%-20s%-10s;-DEPENDENT VARIABLE LABEL-\n' \
            %(uf_d['data']['name'], uf_d['data']['unit'])
        else:
            ufs += '%-30s;-DEPENDENT VARIABLE LABEL-\n'   %uf_d['data']['lbl']
        ufs += '0'.ljust(30) + ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM\n'

# Grids

        for coord in coords:
            ufs += ('   %7d' %dims[coord]).ljust(30) + ';-# OF %1s PTS-\n' %coord
        for coord in coords:
            ufs += rw_for.wr_for(uf_d['grid'][coord]['arr'])

# Data

        ufs += rw_for.wr_for(uf_d['data']['arr'])

    ufs += ';----END-OF-DATA-----------------COMMENTS:-----------\n'
    ufs = ufs.replace('\n', '\n ')
    ufs += comment + '\n'

    f = open(uf, 'w')
    f.write(ufs)
    f.close()
    logger.info('Written file %s' %uf)


class RU:


    def __init__(self, fname):


# === Open the file

        f = open(fname, 'r')
        lines = f.readlines()
        self.label = {}
        self.values = {}
        self.n = {}
        self.vars = ['X', 'Y', 'Z']

# === Headers

        self.nvar = 0
        nlines = len(lines)

        for lin in lines:
            a = lin.split()

            if '-SHOT #' in lin or '; Shot #' in lin:
                self.shot = int(a[0][:5])
                self.dim = int(a[1])

            if '-INDEPENDENT VARIABLE LABEL' in lin or '; Independent variable' in lin:
                for var in self.vars:
                    svar = var + '-'
                    if svar in lin:
                        self.label[var] = lin.split(';-')[0][1: ]
                        self.nvar +=1
                if self.nvar == 0:
                    if '-INDEPENDENT' in lin:
                        self.label['X'] = lin.split(';-')[0][1:]
                    else:
                        self.label['X'] = lin.split('; Independent')[0][1:]
                    self.nvar = 1

            for var in self.vars:
                svar = '-# OF %s PTS-' %var[0]
                if svar in lin:
                    self.n[var] = int(a[0])
                svar = ';-# of radial pts  %s' %var
                if svar in lin:
                    self.n[var] = int(a[0])

            if '-DEPENDENT VARIABLE LABEL-' in lin:
                self.flabel = lin.split(';-')[0][1: ]

# ... Grid data

        if self.dim != self.nvar:
            logger.error('Inconsistency in the number of independent variables dim=%d nvar=%d', self.dim, self.nvar)
            sys.exit()

        list_var = self.vars[ :self.nvar]

        jstart = 3 + 2*self.dim + 2

        while True: # skip over lines until we find first true data line
            try:
                temp = np.array(lines[jstart].split(), dtype=float)
                break
            except:
                jstart += 1

        for var in list_var:
            vtmp = []
            for jlin, lin in enumerate(lines[jstart: ]):
                vtmp += rw_for.ssplit(lin)
                if len(vtmp) == int(self.n[var]):
                    jstart += jlin + 1
                    break
            self.values[var] = np.array(vtmp, dtype=np.float64)

# ... F

        jdata = jstart
        ftmp = []
        for lin in lines[jdata: ]:
            if 'END-OF-DATA' in lin:
                break
            ftmp += rw_for.ssplit(lin)

        farr = np.array(ftmp, dtype=np.float64)

        if self.dim == 1:
            self.fvalues = farr

        elif self.dim == 2:
            self.fvalues = farr.reshape(int(self.n['Y']), int(self.n['X'])).T

        elif self.dim == 3:
            self.fvalues = farr.reshape(int(self.n['Z']), int(self.n['Y']), int(self.n['X'])).T
