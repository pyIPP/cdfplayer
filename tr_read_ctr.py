import os
from scipy.io import netcdf_file
import numpy as np
import mom2rz


class TR_READ_CTR:


    def __init__(self, fcdf, tvec=None, rho=None, nthe=101, endpoint=False):


        if not os.path.isfile(fcdf):
            print('File %s not found' %fcdf)
            return
        else:
            print('Reading Fourier moments from file %s' %fcdf)

        cv = netcdf_file(fcdf, 'r', mmap=False).variables

        ntim, nxb = cv['RMC00'].shape

        tim = cv['TIME3'].data
        xb  = cv['XB'][0, :]

        if tvec is None:
            tind = np.arange(ntim)
        else:
            tvec = np.atleast_1d(tvec)
            tind = np.argmin(np.abs(tim[:, None] - tvec[None, :]), axis=0)

        if rho is None:
            rind = np.arange(nxb)
        else:
            rho  = np.atleast_1d(rho)
            rind = np.argmin(np.abs(xb[:, None] - rho[None, :]), axis=0)

        nt = np.size(tind)
        nrho = np.size(rind)

        self.r_mag = cv['RAXIS'][tind]
        self.z_mag = cv['YAXIS'][tind]
        self.time = tim[tind]

        nmom = 2
        rcl = 'RMC%02d' %nmom
        while rcl in cv.keys():
            nmom += 1
            rcl = 'RMC%02d' %nmom

        rc = np.zeros((nt, nrho, nmom))
        rs = np.zeros((nt, nrho, nmom))
        zc = np.zeros((nt, nrho, nmom))
        zs = np.zeros((nt, nrho, nmom))

        for jmom in range(0, nmom):
            rcl = 'RMC%02d' %jmom
            zcl = 'YMC%02d' %jmom
            rc[:, :, jmom] = cv[rcl][tind, :][:, rind]
            zc[:, :, jmom] = cv[zcl][tind, :][:, rind]
            if jmom > 0:
                rsl = 'RMS%02d' %jmom
                zsl = 'YMS%02d' %jmom
                rs[:, :, jmom] = cv[rsl][tind, :][:, rind]
                zs[:, :, jmom] = cv[zsl][tind, :][:, rind]

        rc = np.squeeze(rc)
        zc = np.squeeze(zc)
        rs = np.squeeze(rs)
        zs = np.squeeze(zs)
        self.Rsurf, self.Zsurf = mom2rz.mom2rz( \
            rc, rs, zc, zs, nthe=nthe, endpoint=endpoint)


if __name__ == "__main__":

    cdf_file = '26333A02.CDF'
    tvec = np.linspace(5, 6, 100)
    READ_EQU(cdf_file, tvec=5)
