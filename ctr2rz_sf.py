import os
from scipy.io import netcdf_file
import numpy as np
import get_sf_grid, tr_read_ctr


def scatter_to_rectangular(r, z, data1, data2, Rmesh, Zmesh):
    """
    Interpolate scattered data points to rectangular grid
    :param r, y: r, y coordinates of data points
    :param data1, data2: input functions, shape as r and z
    :param R, Z: meshgrid of cartesian output grid
    :return: R, Z, interpolated_data (and cache if return_cache)
    """
    import scipy
    from scipy import interpolate

    cache = scipy.spatial.Delaunay(np.vstack((r, z)).T)
    interpolant1 = scipy.interpolate.CloughTocher2DInterpolator(cache, data1)
    interpolant2 = scipy.interpolate.CloughTocher2DInterpolator(cache, data2)
    interp_data1 = np.reshape(interpolant1(np.vstack((Rmesh.flat, Zmesh.flat)).T), Rmesh.shape)
    interp_data2 = np.reshape(interpolant2(np.vstack((Rmesh.flat, Zmesh.flat)).T), Rmesh.shape)

    return interp_data1, interp_data2


class CTR2RZ:


    def __init__(self, cdf_file, tim=None, it=None, Rgrid=None, Zgrid=None, n_the=3601, cliste=True):

# Using m

        runid = cdf_file[-12:-4]
        self.shot = runid[:5]

# Read TRANSP output

        cv = netcdf_file(cdf_file, 'r', mmap=False).variables

# Time point

        if it is None:
            if tim is None:
                it = -1
            else:
                diff = np.abs(cv['TIME3'][:] - tim)
                it = np.argmin(diff)
        self.time = cv['TIME3'][it]

# rho_tor grid

        xb = cv['XB'][it]
        n_xb = len(xb)

        pf1d = 2.*np.pi*cv['PLFLX'][it] # Weber, as CLISTE
        fpol = -0.01*cv['GFUN'][it]*cv['BZXR'][it]    # vs XB; BZXR cm -> m

# Expanding contours from Fourier moments

        fm = tr_read_ctr.TR_READ_CTR(cdf_file, tvec=cv['TIME3'][it], nthe=n_the, endpoint=True)

        self.Rsurf = 0.01*fm.Rsurf
        self.Zsurf = 0.01*fm.Zsurf

        self.R0 = 0.01*fm.r_mag
        self.Z0 = 0.01*fm.z_mag

# Read CLISTE shotfile

        if cliste:
            eq = get_sf_grid.get_grid(runid)
            if eq is None:
                return
            self.Rgrid = eq['Ri'][:, 0] 
            self.Zgrid = eq['Zj'][:, 0]
            nr = len(self.Rgrid)
            nz = len(self.Zgrid)
            tdiff = np.abs(eq['time'] - self.time)
            jt_aug = np.argmin(tdiff)
            self.eq_time = eq['time'][jt_aug] # shotfile time
            print('CDF time = %8.4f, CLISTE time = %8.4f' %(self.time, self.eq_time))
            self.pf_eq_sep = eq['PFxx'][1, jt_aug]
            self.pf_shift = self.pf_eq_sep - pf1d[-1]
            pf_out = eq['PFM'][:nr, :nz, jt_aug]
        else:
            Rmin = min(self.Rsurf[-1]) - 0.02
            Rmax = max(self.Rsurf[-1]) + 0.02
            Zmin = min(self.Zsurf[-1]) - 0.02
            Zmax = max(self.Zsurf[-1]) + 0.02
            nr = len(xb) + 1
            nz = int((Zmax - Zmin)/0.02)
            self.Rgrid = np.linspace(Rmin, Rmax, nr)
            self.Zgrid = np.linspace(Zmin, Zmax, nz)
            pf_out = np.zeros((nr, nz))
            self.pf_shift = 0.

# Output cartesian grid

        self.pfm = np.zeros((nr, nz), dtype=np.float32)
        self.rbp = np.zeros((nr, nz), dtype=np.float32)

# Interpolate fpol to amgnetic axis (XB does not styart from zero)

        rb0 = fpol[1] - (fpol[1] - fpol[0])/(xb[1] - xb[0])*xb[1]
        self.rb = np.append(rb0, fpol)
        self.pf = np.append(0, pf1d) # TRANSP convention, pf=0 in the center

        Rctr = np.zeros((n_xb + 1, n_the-1))
        Zctr = np.zeros((n_xb + 1, n_the-1))
        pf_in = np.repeat(self.pf, n_the-1).T
        rb_in = np.repeat(self.rb, n_the-1).T
        Rctr[0, :] = self.R0
        Zctr[0, :] = self.Z0
        Rctr[1:, :] = self.Rsurf[:, :-1]
        Zctr[1:, :] = self.Zsurf[:, :-1]

        self.X, self.Y = np.meshgrid(self.Rgrid, self.Zgrid)
        self.pfm, self.rbp = scatter_to_rectangular(Rctr.flat, Zctr.flat, pf_in, rb_in, self.X, self.Y)
        self.pfm = self.pfm.T
        self.rbp = self.rbp.T

# Continuity with CLISTE
        self.pfm += self.pf_shift

# clean nan's outside sep
        ind_nan = np.isnan(self.pfm)

        self.pfm[ind_nan] = pf_out[ind_nan] # CLISTE outside
        self.rbp[ind_nan] = 2*self.rb[-1] - self.rb[-2]

        print('Done PFM', self.pfm.shape)

        self.pf += self.pf_shift

# B_R, B_z, B_tor

        print('Calculating B-field components')

        dr = self.Rgrid[1] - self.Rgrid[0]
        dz = self.Zgrid[1] - self.Zgrid[0]

        b_pol = np.gradient(self.pfm, dr, dz)/self.Rgrid[:, None]

        self.b_z =  b_pol[1]/(2.*np.pi)
        self.b_r = -b_pol[0]/(2.*np.pi)
        self.b_t = -self.rbp/self.Rgrid[:, None]
        self.b_z[:,  0] = self.b_z[:,  1]
        self.b_z[:, -1] = self.b_z[:, -2]
        self.b_r[ 0, :] = self.b_r[ 1, :]
        self.b_r[-1, :] = self.b_r[-2, :]


if __name__ == '__main__':


    import matplotlib.pylab as plt

    runid = '29783A01'
    shot = runid[:-3]
    tail = runid[-3:]
    fcdf = '/afs/ipp/home/g/git/tr_client/AUGD/%s/%s/%s.CDF' %(shot, tail, runid)
    rz = CTR2RZ(fcdf, tim=3., n_the=3601)
    psi = rz.pfm[:, :]
    Rax = rz.R0
    zax = rz.Z0
    X, Y = np.meshgrid(rz.Rgrid, rz.Zgrid, indexing='ij')

    b_p = np.hypot(rz.b_r, rz.b_z)

# Plots
    
    figtitle = 'Fields at time_CDF: %8.4f' %rz.time

    fig = plt.figure(1, figsize=(20, 13))
    fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    fig.canvas.set_window_title('Psi contours') 

    plt.subplot(2, 3, 1, aspect='equal')
    plt.title('Psi')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = rz.pf[::]

    plt.contour(X, Y, psi, levels)
    plt.colorbar()
    plt.plot(Rax, zax, 'go', markersize=4)
#    plt.scatter(X, Y, marker='o', color='blue', s=2)
    n_xb = rz.Rsurf.shape[0]
    for jrho in range(n_xb):
        plt.plot(rz.Rsurf[jrho, :], rz.Zsurf[jrho, :], 'g-')

    plt.subplot(2, 3, 4, aspect='equal')
    plt.title('Psi')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(psi), np.max(psi), 20)
    plt.contourf(X, Y, psi, levels)
    plt.colorbar()
#    plt.scatter(X, Y, marker='o', color='blue', s=2)
    for jrho in range(n_xb):
        plt.plot(rz.Rsurf[jrho, :], rz.Zsurf[jrho, :],'g-')

# B field

    n_levels = 21

    plt.subplot(2, 3, 2, aspect='equal')
    plt.title('BR')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(rz.b_r), np.max(rz.b_r), n_levels)
    plt.contourf(X, Y, rz.b_r, levels)
    plt.colorbar()

    plt.subplot(2, 3, 3, aspect='equal')
    plt.title('Bz')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(rz.b_z), np.max(rz.b_z), n_levels)
    plt.contourf(X, Y, rz.b_z, levels)
    plt.colorbar()

    plt.subplot(2, 3, 5, aspect='equal')
    plt.title('Bpol')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(b_p), np.max(b_p), n_levels)
    plt.contourf(X, Y, b_p, levels)
    plt.colorbar()

    plt.subplot(2, 3, 6, aspect='equal')
    plt.title('Btor')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(rz.b_t), np.max(rz.b_t), n_levels)
    plt.contourf(X, Y, rz.b_t, levels)
    plt.colorbar()

    fpdf = 'fields_tb_cliste.pdf'
    plt.savefig(fpdf)
    print('Stored %s' %fpdf)

    plt.show()
