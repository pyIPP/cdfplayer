""" 
Write/Read data in UFILES format (version 0.2)

UFILES is an ASCII format defined in the manual:
  http://w3.pppl.gov/~pshare/help/ufiles_manual.html

How to use ufiles.py:
  
  import ufiles as uf
"""

__author__  = 'Giovanni Tardini (Tel. 1898)'
__version__ = '0.3.0'
__date__    = '16.02.2022'

import datetime, os, sys, shutil, logging
import numpy as np
import rw_for

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('ufiles')
logger.addHandler(hnd)
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

now = datetime.datetime.now()


def split_uname(uname):

    fbase = os.path.basename(uname)
    fname, ext = fbase.split('.')

    for jpos, letter in enumerate(fname):
        if letter.isdigit():
            break
    pre = fname[:jpos]
    shot = int(fname[jpos:])

    return pre, ext, shot


class UFILE(dict):


    def __init__(self, fin=None):

        self.__dict__ = self
        self.f = {}
        if fin is not None:
            self.read(fin)


    def read(self, fin):

# === Open the file

        f = open(fin, 'r')
        lines = f.readlines()
        dims = []
        coords = ['X', 'Y', 'Z']
        self.comment = ''

        self.pre, self.ext, nshot = split_uname(fin)

        logger.debug('Prefix: %s, extension %s', self.pre, self.ext)

# === Headers

        n_coords = 0
        jend = np.nan

        for jline, lin in enumerate(lines):

            if '-SHOT #' in lin or '; Shot #' in lin:
                a = lin.split()
                self.shot = int(a[0][:5])
                logger.debug('Shot %d %d', self.shot, nshot)
                ndim = int(a[1])

            if '-INDEPENDENT VARIABLE LABEL' in lin or '; Independent variable' in lin:
                for coord in coords:
                    svar = coord + '-'
                    if svar in lin:
                        self.__dict__[coord] = {}
                        co_d = self.__dict__[coord]
                        co_d['label'] = lin.split(';-')[0][1: ]
                        co_d['name'] = co_d['label'][:20].strip()
                        co_d['unit'] = co_d['label'][20:].strip()
                        n_coords += 1
                if n_coords == 0:
                    if '-INDEPENDENT' in lin:
                        self.X['label'] = lin.split(';-')[0][1: ]
                    else:
                        self.X['label'] = lin.split('; Independent')[0][1: ]
                    self.X['name'] = self.X['label'][:20].strip()
                    self.X['unit'] = self.X['label'][20:].strip()
                    n_coords = 1

            for var in coords:
                svar1 = '-# OF %s PTS-' %var
                svar2 = ';-# of radial pts  %s' %var
                if svar1 in lin or svar2 in lin:
                   a = lin.split()
                   dims.append(int(a[0]))

            if '-DEPENDENT VARIABLE LABEL-' in lin:
                self.f['label'] = lin.split(';-')[0][1: ]
                self.f['name'] = self.f['label'][:20].strip()
                self.f['unit'] = self.f['label'][20:].strip()
            try:
                flt_lin = rw_for.ssplit(lines[jline])
            except:
                if 'flt_lin' not in locals():
                    jstart = jline+1
                else:
                    if np.isnan(jend):
                        jend = jline

            if jline > jend:
                self.comment += lin
        self.comment = self.comment[1:]
        logger.debug(self.comment)

# ... Grid data

        if ndim != n_coords:
            logger.error('Inconsistency in the number of independent variables dim=%d nvar=%d', ndim, n_coords)
            return

        logger.debug('%d %d %s', jstart, jend, lines[jend])
        data = rw_for.lines2fltarr(lines[jstart: jend])

        ind_split = np.cumsum(dims)
        if ndim == 1:
            self.X['data'], farr = np.split(data, ind_split)
        elif ndim == 2:
            self.X['data'], self.Y['data'], farr = np.split(data, ind_split)
        elif ndim == 3:
            self.X['data'], self.Y['data'], self.Z['data'], farr = \
                np.split(data, ind_split)
        self.f['data'] = farr.reshape(dims[::-1]).T


    def average(self, axis=0):

        self.ext = self.ext + '_AVG%d' %axis
        if axis == 0:
            self.X['data'] = np.atleast_1d(np.nanmean(self.X['data']))
        elif axis == 1:
            self.Y['data'] = np.atleast_1d(np.nanmean(self.Y['data']))
        elif axis == 2:
            self.Z['data'] = np.atleast_1d(np.nanmean(self.Z['data']))
        tmp = np.nanmean(self.f['data'], axis=axis)
        shape = list(self.f['data'].shape)
        shape[axis] = 1
        self.f['data'] = tmp.reshape(shape)


    def uslice(self, val=0., axis=0):

        self.ext = self.ext + '_SLICE%d' %axis
        fun = self.f['data']
        if axis == 0:
            jclose = np.argmin(np.abs(self.X['data'] - val))
            self.X['data'] = np.atleast_1d(self.X['data'][jclose])
            self.f['data'] = fun[jclose: jclose+1]
        elif axis == 1:
            jclose = np.argmin(np.abs(self.Y['data'] - val))
            self.Y['data'] = np.atleast_1d(self.Y['data'][jclose])
            self.f['data'] = fun[:, jclose: jclose+1]
        elif axis == 2:
            jclose = np.argmin(np.abs(self.Y['data'] - val))
            self.Y['data'] = np.atleast_1d(self.Y['data'][jclose])
            self.f['data'] = fun[:, :, jclose: jclose+1]


    def write(self, udir=None, dev='AUGD'):

        self.shot = int(self.shot)
        s5shot = '%05d' %self.shot

        if hasattr(self, 'comment'):
            comment = self.comment
        else:
            comment = ''

        if udir is None:
            udir = '%s/udb/%s' %(os.getenv('HOME'), s5shot)
        os.system('mkdir -p %s' %udir)
        uf = '%s/%s%s.%s' %(udir, self.pre, s5shot, self.ext)
        logger.debug('Output file %s', uf)
        if os.path.isfile(uf):
            ufbak = uf
            jext = 1
            while os.path.isfile(ufbak):
                ufbak = '%s.%d' %(uf, jext)
                jext += 1
            logger.info('mv %s %s' %(uf, ufbak))
            shutil.move(uf, ufbak)

# Dimensions check
        dims = []
        coords = []
        if hasattr(self, 'X'):
            coords.append('X')
            dims.append(len(self.X['data']))
        if hasattr(self, 'Y'):
            coords.append('Y')
            dims.append(len(self.Y['data']))
        if hasattr(self, 'Z'):
            coords.append('Z')
            dims.append(len(self.Z['data']))

        ndim = len(dims)
        if ndim != self.f['data'].ndim:
            logger.error('Data ndim inconsistent with n of independent variables')
            logger.error('No u-file written')
            return

        for jdim, dim in enumerate(dims):
            if dim != self.f['data'].shape[jdim]:
                logger.error('The %dth index has inconsistent grid vs data dimension: %d, %d', \
                    jdim, dim, self.f['data'].shape[jdim] )
                logger.error('No u-file written')
                return

# Header

        ufs = '  %5d%4s %1i 0 6              ;-SHOT #- F(X) DATA -UF%1iDWR- %s\n'\
             %(self.shot,  dev,  ndim,  ndim, now.strftime('%d%b%Y'))
        ufs += ''.ljust(30) + ';-SHOT DATE-  UFILES ASCII FILE SYSTEM\n'

        if hasattr(self, 'scalar'):
            nscal = len(self.scalar)
            ufs += (' %2d' %nscal).ljust(30) + \
                ';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'
            for jscal in range(nscal):
                lbl = self.scalar[jscal]['label']
                val = self.scalar[jscal]['data']
                ufs += ('%11.4E' %val).ljust(30) + ';-SCALAR,  LABEL FOLLOWS:\n'
                ufs += lbl.ljust(30) + '\n'
        else:
            ufs += '  0'.ljust(30) +';-NUMBER OF ASSOCIATED SCALAR QUANTITIES-\n'
        if ndim > 0:
            if (ndim == 1):
                if 'unit' in self.X.keys() and 'name' in self.X.keys():
                    ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL-\n' \
                        %(self.X['name'], self.X['unit'])
                elif 'label' in self.X.keys():
                    ufs += '%-30s;-INDEPENDENT VARIABLE LABEL-\n' %self.X['label']
                else:
                    logger.error('X must have either "name"+"unit" or "label" key')
            else:
                for coord in coords:
                    co_d = self.__dict__[coord]
                    if 'unit' in co_d.keys() and 'name' in co_d.keys():
                        ufs += '%-20s%-10s;-INDEPENDENT VARIABLE LABEL: ' \
                            %(co_d['name'], co_d['unit'])
                    elif 'label' in co_d.keys():
                        ufs += '%-30s;-INDEPENDENT VARIABLE LABEL: ' %co_d['label']
                    else:
                        logger.error('%s must have either "name"+"unit" or "label" key', coord)
                    ufs += '%1s-\n' %coord

        if 'unit' in self.f.keys():
            ufs += '%-20s%-10s;-DEPENDENT VARIABLE LABEL-\n' \
            %(self.f['name'], self.f['unit'])
        else:
            ufs += '%-30s;-DEPENDENT VARIABLE LABEL-\n'   %self.f['label']
        ufs += '0'.ljust(30) + ';-PROC CODE- 0:RAW 1:AVG 2:SM. 3:AVG+SM\n'

# Grids

        for jcoord, coord in enumerate(coords):
            ufs += ('   %7d' %dims[jcoord]).ljust(30) + ';-# OF %1s PTS-\n' %coord
        for coord in coords:
            data = self.__dict__[coord]['data']
            ufs += rw_for.wr_for(data)

# Data

        ufs += rw_for.wr_for(self.f['data'])

        ufs += ';----END-OF-DATA-----------------COMMENTS:-----------\n'
        ufs = ufs.replace('\n', '\n ')
        ufs += comment + '\n'

        f = open(uf, 'w')
        f.write(ufs)
        f.close()
        logger.info('Written file %s' %uf)


    def plot(self):

        import matplotlib.pylab as plt

        data = self.f['data']
        if data.ndim == 1:
            plt.figure('u-file', (11, 8))
            plt.plot(self.X['data'], data)
            plt.xlabel('%s [%s]' %(self.X['name'], self.X['unit'] ) )
            plt.ylabel('%s [%s]' %(self.f['name'], self.f['unit'] ) )
        elif data.ndim == 2:
            X, Y = np.meshgrid(self.X['data'], self.Y['data'], indexing='ij')
            plt.figure('u-file', (7, 8))
            plt.subplot(1, 1, 1)
            plt.contourf(X, Y, data, levels=20)
            plt.colorbar()
            plt.xlabel('%s [%s]' %(self.X['name'], self.X['unit'] ) )
            plt.ylabel('%s [%s]' %(self.Y['name'], self.Y['unit'] ) )
            plt.title('%s [%s]' %(self.f['name'], self.f['unit'] ) )
        else:
            logger.warning('No plots for 0D nor 3D ufiles')
        plt.show()


if __name__ == '__main__':


    ufin  = 'I38335.CEZCMZ'
#    ufin  = 'M38335.MRY'
    u = UFILE(fin=ufin)
    print(u['f']['data'][20, :])
    u.plot()
    u.write()
    u.write(avg_axis=0)
