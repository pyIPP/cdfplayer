from scipy.io import netcdf
import numpy as np


class READ_EQU:


    def __init__(self, cdf_file, tvec=None, rho=None):

        cv1 = netcdf.netcdf_file(cdf_file, 'r', mmap=False).variables
        cv = {}
        for key, val in cv1.items():
            cv[key] = np.atleast_1d(val.data)

        ntim, nr = cv['RMC00'].shape

        tim = cv['TIME3']
        xb  = cv['XB'][0, :]

        if tvec is None:
            tind = np.arange(ntim)
        else:
            tvec = np.atleast_1d(tvec)
            tind = np.argmin(np.abs(tim[:, None] - tvec[None, :]), axis=0)

        if rho is None:
            rind = np.arange(nr)
        else:
            rho  = np.atleast_1d(rho)
            rind = np.argmin(np.abs(xb[:, None] - rho[None, :]), axis=0)

        nt   = np.size(tind)
        nrho = np.size(rind)

        self.r_mag = cv['RAXIS'][tind]
        self.z_mag = cv['YAXIS'][tind]
        self.time = tim[tind]
        nmom = 0
        rcl = 'RMC00'
        while rcl in cv.keys():
            nmom += 1
            rcl = 'RMC0%d' %nmom

        rc = np.zeros((nt, nrho, nmom))
        rs = np.zeros((nt, nrho, nmom))
        zc = np.zeros((nt, nrho, nmom))
        zs = np.zeros((nt, nrho, nmom))

        for jmom in range(0, nmom):
            rcl = 'RMC0%d' %jmom
            zcl = 'YMC0%d' %jmom
            rc[..., jmom] = cv[rcl][:, rind][tind, :]
            zc[..., jmom] = cv[zcl][:, rind][tind, :]

            if jmom > 0:
                rsl = 'RMS0%d' %jmom
                zsl = 'YMS0%d' %jmom
                rs[..., jmom] = cv[rsl][:, rind][tind, :]
                zs[..., jmom] = cv[zsl][:, rind][tind, :]

        self.rc = np.squeeze(rc)
        self.zc = np.squeeze(zc)
        self.rs = np.squeeze(rs)
        self.zs = np.squeeze(zs)
        print('Eq-read done!')


if __name__ == "__main__":

    cdf_file = '26333A02.CDF'
    tvec = np.linspace(5, 6, 100)
    READ_EQU(cdf_file, tvec=5)
