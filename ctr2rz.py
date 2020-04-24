import os
from scipy.io import netcdf
import numpy as np
import get_sf_grid, tr_read_ctr, map_ctr


class CTR2RZ:


    def __init__(self, cdf_file, tim=None, it=None, Rgrid=None, Zgrid=None, n_the=3601):

# Using m

        runid = cdf_file[-12:-4]

# Read TRANSP output

        cv = netcdf.netcdf_file(cdf_file, 'r', mmap=False).variables

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

        pf1d = -2.*np.pi*cv['PLFLX'][it]
        fpol = -0.01*cv['GFUN'][it]*cv['BZXR'][it]    # vs XB; BZXR cm -> m

# Expanding contours from Fourier moments

        fm = tr_read_ctr.TR_READ_CTR(cdf_file, tvec=cv['TIME3'][it], nthe=n_the, endpoint=True)

        self.Rsurf = 0.01*fm.Rsurf
        self.Zsurf = 0.01*fm.Zsurf

        self.R0 = 0.01*fm.r_mag
        self.Z0 = 0.01*fm.z_mag

# Read CLISTE shotfile

        eq = get_sf_grid.get_grid(runid)
        self.Rgrid = eq['Ri'][:, 0] 
        self.Zgrid = eq['Zj'][:, 0]
        nr = len(self.Rgrid)
        nz = len(self.Zgrid)
        tdiff = np.abs(eq['time'] - self.time)
        jt_aug = np.argmin(tdiff)
        self.eq_time = eq['time'][jt_aug] # shotfile time
        print('CDF time = %8.4f, CLISTE time = %8.4f' %(self.time, self.eq_time))
        self.pf_eq_sep = eq['PFxx'][1, jt_aug]
        pf_out = eq['PFM'][:, :, jt_aug]
        self.pf_shift = self.pf_eq_sep - pf1d[-1]

# Output cartesian grid

        self.pfm = np.zeros((nr, nz), dtype=np.float32)
        self.rbp = np.zeros((nr, nz), dtype=np.float32)

# Interpolate fpol to amgnetic axis (XB does not styart from zero)

        rb0 = fpol[1] - (fpol[1] - fpol[0])/(xb[1] - xb[0])*xb[1]
        self.rb = np.append(rb0, fpol)
        self.pf = np.append(0, pf1d) # TRANSP convention, pf=0 in the center

# Reference frame centerd in Rmag, zmag

        Rctr = np.zeros((n_xb + 1, n_the-1))
        Zctr = np.zeros((n_xb + 1, n_the-1))
        Rctr[1:, :] = self.Rsurf[:, :-1] - self.R0
        Zctr[1:, :] = self.Zsurf[:, :-1] - self.Z0
        Rrect = self.Rgrid - self.R0
        Zrect = self.Zgrid - self.Z0

        Rmesh, Zmesh = np.meshgrid(Rrect, Zrect)
        th_mesh = np.arctan2(Zmesh.T, Rmesh.T)
        norm_mesh = np.hypot(Rmesh.T, Zmesh.T)
        norm3 = np.zeros(3)
        iclose = map_ctr.map_norm(Rrect, Zrect, Rctr, Zctr) # bottle-neck

        print('Mapping PFM')

        i_left = np.maximum(iclose-1, 0)
        i_left = np.minimum(i_left, n_xb -2)
        for jr in range(nr):
            for jz in range(nz):
                inext = iclose[jr, jz]
                norm_sur = map_ctr.interp_norm(Rctr[inext], Zctr[inext], th_mesh[jr, jz])
                if (inext == n_xb) and (norm_mesh[jr, jz] >= norm_sur): # outside separatrix
                    self.pfm[jr, jz] = pf_out[jr, jz]  # Take from CLISTE
                    self.rbp[jr, jz] = fpol[0]
                else:
                    ileft = i_left[jr, jz]
                    for jrho in range(3):
                        if (ileft + jrho == 0):
                            norm3[jrho] = 0
                        elif (ileft + jrho == inext):
                            norm3[jrho] = norm_sur # already computed
                        else:
                            norm3[jrho] = map_ctr.interp_norm(Rctr[ileft+jrho], Zctr[ileft+jrho], th_mesh[jr, jz])

                    self.pfm[jr, jz] = map_ctr.quad_int(norm_mesh[jr, jz], norm3, self.pf[ileft: ileft+3])
                    self.rbp[jr, jz] = map_ctr.quad_int(norm_mesh[jr, jz], norm3, self.rb[ileft: ileft+3])
                    self.pfm[jr, jz] += self.pf_shift

        self.pf += self.pf_shift

# B_R, B_z, B_tor

        print('Calculating B-field components')

        dR  = np.gradient(self.Rgrid)
        dz  = np.gradient(self.Zgrid)

        d_dr = np.gradient(self.pfm, axis=0)/dR[:, None]
        d_dz = np.gradient(self.pfm, axis=1)/dz[None, :]
        self.b_z =  .5/np.pi*d_dz/self.Rgrid[:, None]
        self.b_r = -.5/np.pi*d_dr/self.Rgrid[:, None]
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
    psi = rz.pfm[:, :].T
    Rax = rz.R0
    zax = rz.Z0
    X, Y = np.meshgrid(rz.Rgrid, rz.Zgrid)

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
    levels = rz.pf[::-1]

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
    plt.contourf(X, Y, rz.b_r.T, levels)
    plt.colorbar()

    plt.subplot(2, 3, 3, aspect='equal')
    plt.title('Bz')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(rz.b_z), np.max(rz.b_z), n_levels)
    plt.contourf(X, Y, rz.b_z.T, levels)
    plt.colorbar()

    plt.subplot(2, 3, 5, aspect='equal')
    plt.title('Bpol')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(b_p), np.max(b_p), n_levels)
    plt.contourf(X, Y, b_p.T, levels)
    plt.colorbar()

    plt.subplot(2, 3, 6, aspect='equal')
    plt.title('Btor')
    plt.xlabel('R [cm]')
    plt.ylabel('z [cm]')
    levels = np.linspace(np.min(rz.b_t), np.max(rz.b_t), n_levels)
    plt.contourf(X, Y, rz.b_t.T, levels)
    plt.colorbar()

    fpdf = 'fields_tb_cliste.pdf'
    plt.savefig(fpdf)
    print('Stored %s' %fpdf)

    plt.show()
