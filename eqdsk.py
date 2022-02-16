import os, logging
import numpy as np
import rw_for

fmt = logging.Formatter('%(asctime)s | %(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('eqdsk')
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.INFO)

Rlim_aug = [103.5, 105, 109, 114, 123.5, 134, 143, 150.5, 163, 177, 197, 211, \
217, 220.5, 221, 217.5, 213, 201, 189.0, 170.5, 164.5, 158, 145.8, 132.5, \
123.5,  128.6, 128,   124.5, 112.5,  106, 103.5]
Zlim_aug = [0,  25,  50,  70,  96.5,  110, 115, 119,   115, 108.5, 78, 54,  39, \
21, 0,   -15,   -28, -50, -68,   -86,   -96.6, -120.8, -106,  -106, \
-112.6, -97.6, -89.2, -82, -63.4, -30, 0]

# COCOs table I page 8

sigma = { \
    'bp'       : [ 1,  1, -1, -1,  1,  1, -1, -1], \
    'rphiz'    : [ 1, -1,  1, -1,  1, -1,  1, -1], \
    'rhothephi': [ 1,  1, -1, -1, -1, -1,  1,  1]}

explain = {'rphiz'    : {1: '(R, phi, Z) r'    , -1: '(R, Z, phi) l'}, \
           'rhothephi': {1: '(rho, the, phi) r', -1: '(rho, phi, the) l'    } }

for key, val in sigma.items():
    sigma[key] = np.array(val)


class EQDSK(dict):


    def __init__(self, GEQ=None):

        if GEQ is None:
            self.__dict__ = self
        else:
            self.__dict__ = GEQ


    def write(self, f_eqdsk):

        mixed = ('FPOL', 'PRES', 'FFPRIM', 'PPRIME')

        fmt  = '%16.9E'
        nlin = 5
        format_str = (fmt * 5 + '\n')
        xdum = 0.
        idum = 0

        f = open(f_eqdsk,  'w')
# Scalars
        f.write('%s%s%s%s%s%s %d %d %d\n' %(*self.CASE, idum, self.NW, self.NH) )
        f.write( format_str % (self.RDIM, self.ZDIM, self.RMAXIS, self.RLEFT, self.ZMID) )
        f.write( format_str % (self.RMAXIS, self.ZMAXIS, self.SIMAG, self.SIBRY, self.BCENTR) )
        f.write( format_str % (self.CURRENT, self.SIMAG, xdum, self.RMAXIS, xdum) )
        f.write( format_str % (self.ZMAXIS, xdum, self.SIBRY, xdum, xdum) )
# Profiles
        f.write( rw_for.wr_for(self.FPOL  , fmt=fmt, n_lin=nlin) )
        f.write( rw_for.wr_for(self.PRES  , fmt=fmt, n_lin=nlin) )
        f.write( rw_for.wr_for(self.FFPRIM, fmt=fmt, n_lin=nlin) )
        f.write( rw_for.wr_for(self.PPRIME, fmt=fmt, n_lin=nlin) )
        f.write( rw_for.wr_for(self.PSIRZ , fmt=fmt, n_lin=nlin) ) # Matrix
        f.write( rw_for.wr_for(self.QPSI  , fmt=fmt, n_lin=nlin) )

        if not hasattr(self, 'RLIM'):
            self.RLIM = 0.01*np.array(Rlim_aug)
            self.ZLIM = 0.01*np.array(Zlim_aug)
        n_lim = len(self.RLIM)
        n_the = len(self.RBBBS)
        f.write( '%5i%5i\n' % (n_the, n_lim) )
        tmp = np.append(self.RBBBS, self.ZBBBS)
        rz  = tmp.reshape(2,  n_the)
        f.write( rw_for.wr_for(rz, fmt=fmt, n_lin=nlin) )
        tmp = np.append(np.array(self.RLIM), np.array(self.ZLIM))
        rz2 = tmp.reshape(2,  n_lim)
        f.write( rw_for.wr_for(rz2, fmt=fmt, n_lin=nlin) )
        f.close()

        logger.info('Output stored in %s' %f_eqdsk)


    def read(self, f_eqdsk):

        if not os.path.isfile(f_eqdsk):
            logger.error('File %s not found', f_eqdsk)
        else:
            logger.info('Parsing EQDSK file %s', f_eqdsk)

        f = open(f_eqdsk, 'r')
        lines = f.readlines()

        header = lines[0].split()
        self.NW = int(header[-2])
        self.NH = int(header[-1])

        self.RDIM   , self.ZDIM  , self.RCENTR, self.RLEFT, self.ZMID   = rw_for.ssplit(lines[1])
        self.RMAXIS , self.ZMAXIS, self.SIMAG , self.SIBRY, self.BCENTR = rw_for.ssplit(lines[2])
        self.CURRENT, self.SIMAG , xdum       , self.RMAXIS, xdum       = rw_for.ssplit(lines[3])
        self.ZMAXIS , xdum       , self.SIBRY , xdum       , xdum       = rw_for.ssplit(lines[4])

        data = rw_for.lines2fltarr(lines[5:])
 
        jprof_end = self.NW * (5 + self.NH)
        ind_split = [self.NW, 2*self.NW, 3*self.NW, 4*self.NW, \
            self.NW * (4 + self.NH)]
        self.FPOL, self.FFPRIM, self.PRES, self.PPRIME, \
            self.PSIRZ, self.QPSI = np.split(data[: jprof_end], ind_split)
        self.PSIRZ = self.PSIRZ.reshape(self.NH, self.NW).T

# Read bndy, lim
        nbdy = int(data[jprof_end])
        nlim = int(data[jprof_end+1])
        jbdy = jprof_end + 2
        jlim = jbdy + 2*nbdy
        bdy = data[jbdy: jlim]
        self.RBBBS = bdy[ ::2]
        self.ZBBBS = bdy[1::2]
        lim = data[jlim: jlim + 2*nlim]
        self.RLIM = lim[ ::2]
        self.ZLIM = lim[1::2]
        print('bdy text', nbdy, nlim, jprof_end, jbdy, jlim, len(data))
        print(self.NW, self.NH, self.NW*self.NH, data[jprof_end-1])
        print(self.QPSI)
        self.get_coco()


    def plot(self):

        import matplotlib.pylab as plt

        Rgrid = np.linspace(self.RLEFT, self.RLEFT + self.RDIM, self.NW)
        Zgrid = np.linspace(self.ZMID - 0.5*self.ZDIM, self.ZMID + 0.5*self.ZDIM, self.NH)
        psi_grid = np.linspace(self.SIMAG, self.SIBRY, self.NW)

        X, Y = np.meshgrid(Rgrid, Zgrid)

        plt.figure('EQDSK pfm', (7, 8))
        plt.subplot(1, 1, 1, aspect='equal')
        plt.contourf(X, Y, self.PSIRZ.T, levels=20)
        plt.colorbar()
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('Pol. flux')
        plt.plot(self.RBBBS , self.ZBBBS , 'b-')
        plt.plot(self.RLIM  , self.ZLIM  , 'k-')
        plt.plot(self.RMAXIS, self.ZMAXIS, 'ro')

        plt.figure('EQDSK profs', (10, 8))
        plt.subplot(3, 2, 1)
        plt.plot(psi_grid[:-1], self.FPOL[:-1])
        plt.title('F')
        plt.subplot(3, 2, 2)
        plt.plot(psi_grid[:-1], self.FFPRIM[:-1])
        plt.title('dF/dpsi')
        plt.subplot(3, 2, 3)
        plt.plot(psi_grid[:-1], self.PRES[:-1])
        plt.title('Pres')
        plt.subplot(3, 2, 4)
        plt.plot(psi_grid[:-1], self.PPRIME[:-1])
        plt.title('dPres/dpsi')
        plt.subplot(3, 2, 5)
        plt.plot(psi_grid[:-1], self.QPSI[:-1])
        plt.title('Saf. factor')

        plt.show()


    def get_coco(self, ip_shot='ccw', bt_shot='cw'):
        """
        Returns the COCO number of a given eqdsk object
        """

# dpsi_sign positive if psi_sep > psi0
        dpsi_sign = np.sign(self.SIBRY - self.SIMAG)
# Known plasma discharge
        ccw_ip = 1 if (ip_shot == 'ccw') else -1 # AUG default: 1
        ccw_bt = 1 if (bt_shot == 'ccw') else -1 # AUG default: -1

# Table III

        sign_q  = np.sign(np.nanmean(self.QPSI))
        sign_ip = np.sign(self.CURRENT)
        sign_bt = np.sign(self.BCENTR)
        sigma_rphiz = sign_ip*ccw_ip
        sigma_bp    = dpsi_sign*sign_ip
# Eq 45
        sigma_rhothephi = sign_q*sign_ip*sign_bt
        logger.debug('%d, %d, %d', sigma_bp, sigma_rphiz, sigma_rhothephi)
        for jc, rhothephi in enumerate(sigma['rhothephi']):
            if( sigma['bp'   ][jc] == sigma_bp    and \
                sigma['rphiz'][jc] == sigma_rphiz and \
                rhothephi          == sigma_rhothephi):
                self.coco = jc + 1
                break

# Find out 2*pi factor for Psi

        psi1d = np.linspace(self.SIMAG, self.SIBRY, self.NW)
        dpsi = np.gradient(psi1d)

# Radial average
# It is either q_ratio ~ 1 (COCO > 10) or ~ 2*pi (COCO < 10)

        rmin = 0.5*(np.max(self.RBBBS) - np.min(self.RBBBS)) 
        q_estimate = np.abs((np.pi * self.BCENTR * rmin**2) / (self.SIBRY - self.SIMAG))
        qabs = np.abs(self.QPSI[-1])
        if np.abs(q_estimate - qabs) < np.abs(q_estimate/(2*np.pi) - qabs):
            self.coco += 10
        logger.debug('%8.4f %8.4f', np.abs(q_estimate - qabs), np.abs(q_estimate/(2*np.pi) - qabs) )
        logger.info('COCO number: %d', self.coco)


if __name__ == '__main__':

    feqdsk = '28053_1.2003s.eqdsk'
    feqdsk = '23076W01.eqdsk'
    eq = EQDSK()
    eq.read(feqdsk)
    eq.plot()
