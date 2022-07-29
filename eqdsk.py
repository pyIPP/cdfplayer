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


    def __init__(self, fin=None):

        self.__dict__ = self
        if fin is not None:
            self.read(fin)


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
        self.FPOL, self.PRES, self.FFPRIM, self.PPRIME, \
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
        self.Rgrid = np.linspace(self.RLEFT, self.RLEFT + self.RDIM, self.NW)
        self.Zgrid = np.linspace(self.ZMID - 0.5*self.ZDIM, self.ZMID + 0.5*self.ZDIM, self.NH)
        self.psi1d = np.linspace(self.SIMAG, self.SIBRY, self.NW)
        self.find_coco()


    def write(self, f_eqdsk, geq=None):

        if geq is not None:
            self.__dict__ = geq

        mixed = ('FPOL', 'PRES', 'FFPRIM', 'PPRIME')

        fmt  = '%16.9E'
        nlin = 5
        format_str = (fmt * 5 + '\n')
        xdum = 0.
        idum = 0

        f = open(f_eqdsk,  'w')
# Scalars
        if hasattr(self, 'CASE2'):
            f.write('%s %d %d %d\n' %(self.CASE2, idum, self.NW, self.NH) )
        else:
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


    def plot(self):

        import matplotlib.pylab as plt
        from matplotlib import cm

        psi_grid = np.linspace(self.SIMAG, self.SIBRY, self.NW)

        X, Y = np.meshgrid(self.Rgrid, self.Zgrid)

        plt.figure('EQDSK', (14, 8))
        plt.subplots_adjust(top=0.95, bottom=0.1, left=0.1, right=0.95, hspace=0, wspace=0.5)

        plt.subplot2grid((5, 2), (0, 0), rowspan=5, aspect='equal')
        plt.contourf(X, Y, self.PSIRZ.T, levels=20, cmap=cm.jet)
        plt.colorbar()
        plt.xlabel('R [m]')
        plt.ylabel('Z [m]')
        plt.title('Pol. flux')
        plt.plot(self.RBBBS , self.ZBBBS , 'b-')
        plt.plot(self.RLIM  , self.ZLIM  , 'k-')
        plt.plot(self.RMAXIS, self.ZMAXIS, 'ro')

        plt.subplot2grid((5, 2), (0, 1))
        plt.plot(psi_grid, self.FPOL)
        plt.ylabel('F')
        axp = plt.gca()
        axp.tick_params(axis='x', which='both', bottom='on', top='on', labelbottom='off')
        axp.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axp.yaxis.major.formatter._useMathText = True

        plt.subplot2grid((5, 2), (1, 1))
        plt.plot(psi_grid, self.FFPRIM)
        plt.ylabel('dF/dpsi')
        axp = plt.gca()
        axp.tick_params(axis='x', which='both', bottom='on', top='on', labelbottom='off')
        axp.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axp.yaxis.major.formatter._useMathText = True

        plt.subplot2grid((5, 2), (2, 1))
        plt.plot(psi_grid, self.PRES)
        plt.ylabel('Pres [Pa]')
        axp = plt.gca()
        axp.tick_params(axis='x', which='both', bottom='on', top='on', labelbottom='off')
        axp.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axp.yaxis.major.formatter._useMathText = True

        plt.subplot2grid((5, 2), (3, 1))
        plt.plot(psi_grid, self.PPRIME)
        plt.ylabel('dPres/dpsi')
        axp = plt.gca()
        axp.tick_params(axis='x', which='both', bottom='on', top='on', labelbottom='off')
        axp.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        axp.yaxis.major.formatter._useMathText = True

        if not hasattr(self, 'cocos'):
            self.find_coco()
        if self.cocos < 10:
            unit = 'Weber/rad'
        else:
            unit = 'Weber'

        plt.subplot2grid((5, 2), (4, 1))
        plt.plot(psi_grid, self.QPSI)
        plt.ylabel('Saf. factor')
        plt.xlabel(r'$\Psi$ [%s]' %unit)

        plt.show()


    def find_coco(self, ip_shot='ccw', bt_shot='cw'):
        """
        Returns the COCO number of a given eqdsk object
        """

# Known plasma discharge
        self.ip_ccw = 1 if (ip_shot == 'ccw') else -1 # AUG default: 1
        self.bt_ccw = 1 if (bt_shot == 'ccw') else -1 # AUG default: -1

# Table III
# dpsi_sign positive if psi_sep > psi0
        dpsi_sign = np.sign(self.SIBRY - self.SIMAG)
        sign_q    = np.sign(np.nanmean(self.QPSI))
        sign_ip   = np.sign(self.CURRENT)
        sign_bt   = np.sign(self.BCENTR)
        logger.debug('dpsi: %d, q: %d, Ip: %d, Bt: %d', dpsi_sign, sign_q, sign_ip, sign_bt)  
        sigma_rphiz = sign_ip*self.ip_ccw # If sigma_rphiz>0, ref frame is ccw
        if sign_bt*self.bt_ccw != sigma_rphiz:
            logger.error('Inconsistent sign of Bt')
        sigma_bp = dpsi_sign*sign_ip
# Eq 45
        sigma_rhothephi = sign_q*sign_ip*sign_bt
        logger.debug('%d, %d, %d', sigma_bp, sigma_rphiz, sigma_rhothephi)
        for jc in range(len(sigma['rhothephi'])):
            if( sigma['bp'   ][jc]     == sigma_bp    and \
                sigma['rphiz'][jc]     == sigma_rphiz and \
                sigma['rhothephi'][jc] == sigma_rhothephi):
                self.cocos = jc + 1
                break

# Find out 2*pi factor for Psi

        psi1d = np.linspace(self.SIMAG, self.SIBRY, self.NW)
        dpsi = np.gradient(psi1d)

        rmin = 0.5*(np.max(self.RBBBS) - np.min(self.RBBBS)) 
        q_estimate = np.abs((np.pi * self.BCENTR * rmin**2) / (self.SIBRY - self.SIMAG))
        qabs = np.abs(self.QPSI[-1])
        if np.abs(q_estimate - qabs) < np.abs(q_estimate/(2*np.pi) - qabs):
            self.cocos += 10

        jcoco = self.cocos%10 - 1
        ebp = self.cocos//10
        self.psi_sign = -self.ip_ccw*sigma['bp'][jcoco]*sigma['rphiz'][jcoco]
        psi_2pi  = (2.*np.pi)**ebp
        self.psi_fac = self.psi_sign*psi_2pi
        self.phi_sign = self.bt_ccw*sigma['rphiz'][jcoco]

        logger.debug('%8.4f %8.4f', np.abs(q_estimate - qabs), np.abs(q_estimate/(2*np.pi) - qabs) )
        logger.info('COCO number: %d', self.cocos)


    def to_coco(self, cocos_out=11, find=True):

        if hasattr(self, 'cococ') and not find:
            cocos_in = self.cocos
        else:
            self.find_coco()
            cocos_in = self.cocos

        logger.info('COCOS conversion from %d to %d' %(cocos_in, cocos_out))
        jc_in   = cocos_in %10 - 1
        jc_out  = cocos_out%10 - 1
        ebp_in  = cocos_in//10
        ebp_out = cocos_out//10

# Equation 9, table I, equation 39, 45
        q_sign   = sigma['rhothephi'][jc_in] * sigma['rhothephi'][jc_out]
        phi_sign = sigma['rphiz'][jc_in]*sigma['rphiz'][jc_out]
        psi_sign = sigma['rphiz'][jc_in]*sigma['bp'][jc_in] * sigma['rphiz'][jc_out]*sigma['bp'][jc_out]
        psi_2pi  = (2.*np.pi)**(ebp_out - ebp_in)
        psi_fac = psi_sign*psi_2pi

        for key, val in self.__dict__.items():
            if val is None:
                continue
            if key in ('BCENTR', 'FPOL', 'CURRENT'):
                val *= phi_sign
            elif key in ('PSIRZ', 'SIMAG', 'SIBRY'):
                val *= psi_fac
            elif key in ('PPRIME', 'FFPRIM'):
                val /= psi_fac
            elif key in ('QPSI', ):
                val *= q_sign


    def brzt(self):

        from scipy.interpolate import interp1d

        dr = self.Rgrid[1] - self.Rgrid[0]
        dz = self.Zgrid[1] - self.Zgrid[0]
        fBt = interp1d(self.psi1d, self.FPOL, kind='linear', fill_value='extrapolate')

        Bpol = np.gradient(self.PSIRZ, dr, dz)/self.Rgrid[:, None]
        self.Br = Bpol[1]
        self.Bz = -Bpol[0]
        self.Bt = fBt(self.PSIRZ)/self.Rgrid[:, None]

        self.Br = -Bpol[1]/abs(self.psi_fac)
        self.Bz =  Bpol[0]/abs(self.psi_fac)
        self.Bt *= -self.psi_sign*self.phi_sign



if __name__ == '__main__':

    feqdsk = '/afs/ipp/home/g/git/python/maingui/28053_1.2003s.eqdsk'
#    feqdsk = '23076W01.eqdsk'
    eq = EQDSK()
    eq.read(feqdsk)
    eq.plot()
    eq.to_coco(cocos_out=13)
    eq.find_coco()
    eq.plot()
