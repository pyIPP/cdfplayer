from netCDF4 import Dataset
import matplotlib.pylab as plt
import numpy as np
import sys, os
import read_eq_sf, read_equ, mom2rz


def quad_int(x_out, x_in, y_in):

    if len(x_in) != 3:
        print('Error, need exactly 3 points, not %d' %len(x_in) )
        return None
    if x_in[0] == x_in[1] or x_in[0] == x_in[2] or x_in[1] == x_in[2]:
        print('Need 3 different x points', x_in)
        return None

    ratio = (y_in[2] - y_in[1])/(x_in[2] - x_in[1])
    A     = (y_in[0] - y_in[1])/(x_in[0] - x_in[1]) - ratio
    A *= 1./(x_in[0] - x_in[2])
    B = ratio - A*(x_in[2] + x_in[1])
    C = y_in[0] - A*x_in[0]**2 - B*x_in[0]

    return A*x_out**2 + B*x_out + C


def close_theta(x_in, y_in, th_ref):

    nthe = len(x_in)
    f = np.arctan(y_in/x_in)
    ind = np.where(x_in < 0)
    f[ind] += np.pi
    ind = np.where(f < 0)
    f[ind] += 2*np.pi

# Discontinuity at 2pi -> 0
    f += -th_ref
    i0 = np.argmin(np.abs(f))
    i1 = i0 - 1
    i2 = i0 + 1
    if i0 == 0:
        i1 = nthe - 1
        f[i1] += -2.*np.pi
    if i0 == nthe-1:
        i2 = 0
        f[i2] += 2.*np.pi

    ind = [i1, i0, i2]
    g = f[ind]
    rho_zero = np.hypot(x_in[ind], y_in[ind])

    return quad_int(0., g, rho_zero)


class RZ2PSI:


    def __init__(self, cdf_file, it=-1, vplot=True, n_the=3601):

        runid = cdf_file[-12:-4]
        eq = read_eq_sf.get_grid(runid)

# CLISTE shotfile, R(t) and z(t) grid are fixed in time
        Rgrid = 100*eq['Ri'][:, 0] 
        zgrid = 100*eq['Zj'][:, 0]
        pf_sep = eq['PFxx'][1]
        nr = len(Rgrid)
        nz = len(zgrid)

#==========
# Profiles
#==========

        print(cdf_file)
        cdf = Dataset(cdf_file, 'r', format='NETCDF4')
        cv = cdf.variables

        n_xb = cv['XB'].shape[1]

        data = {}
        sgr = ('PLFLX2PI', 'BPOL', 'CUR', 'GFUN')
        sig = ('RAXIS', 'YAXIS', 'BZXR')
        for lbl in sig + sgr:
            data[lbl] = cv[lbl][it]
        data['PLFLX2PI'] *= -1.
        data['fpol'] = data['GFUN']         # vs XB
        data['BZXR'] *= -0.01
        data['fpol'] *= data['BZXR']

#======================
# Poloidal flux matrix 
#======================

        fm = read_equ.READ_EQU(cdf_file, tvec=cv['TIME3'][it])

        r_plot, z_plot = mom2rz.mom2rz(fm.rc, fm.rs, \
                                       fm.zc, fm.zs, nthe=n_the)

        self.pfm = np.zeros((nr, nz), dtype=np.float32)
        rbp = np.zeros((nr, nz), dtype=np.float32)

        th_step = 2*int(float(n_the)/float(nz))
        rh_step = int(3*float(n_xb)/float(nz)) + 2

        print('Time = %8.4f' %cv['TIME3'][it])
        diff = np.abs(eq['time'] - cv['TIME3'][it])
        jt_aug = np.argmin(diff)
        R0 = data['RAXIS']
        z0 = data['YAXIS']
        pf = data['PLFLX2PI']
        pf_out = eq['PFM'][:, :, jt_aug]
        self.pf_shift = pf_sep[jt_aug] - pf[-1]
        rb = data['fpol']
        for jr in range(nr):
            Rgrd = Rgrid[jr]
            for jz in range(nz):
                zgrd = zgrid[jz]
                thgrd = np.arctan((zgrd - z0)/(Rgrd - R0))
                if Rgrd < R0:
                    thgrd += np.pi
                if thgrd < 0:
                    thgrd += 2.*np.pi
                if jz == 0:
                    rh1 = n_xb - 1
                    rh2 = n_xb
                    th1 = 0
                    th2 = n_the
                else:
                    rh1 = max(0, irho - rh_step)
                    rh2 = min(n_xb, irho + rh_step)
                    th1 = max(0, ithe - th_step)
                    th2 = min(n_the, ithe + th_step)
                    if (th2 >= n_the-th_step) or (th1 <= th_step) or (irho <= 5*rh_step):
                        th1 = 0
                        th2 = n_the
                rdist = (r_plot[rh1:rh2, th1:th2] - Rgrd)**2
                zdist = (z_plot[rh1:rh2, th1:th2] - zgrd)**2
                dist = rdist + zdist
                rz_min = np.where(dist == np.min(dist))
                irho = rh1 + rz_min[0][0]
                ithe = th1 + rz_min[1][0]

                norm_grd = np.hypot(Rgrd - R0, zgrd - z0) # scalar
                rtmp = r_plot[irho, th1:th2] - R0
                ztmp = z_plot[irho, th1:th2] - z0
                norm_sur = close_theta(rtmp, ztmp, thgrd)
                if (irho == n_xb-1) and (norm_grd >= norm_sur): #out of separatrix
                    self.pfm[jr, jz] = pf_out[jr, jz]
                    rbp[jr, jz] = data['fpol'][0]
                else:
                    norm = np.zeros(3)
                    if irho == 0:
                        norm[0] = 0
                        j = 1
                        for jrho in range(0,2):
                            rtmp = r_plot[jrho, th1:th2] - R0
                            ztmp = z_plot[jrho, th1:th2] - z0
                            norm[j] = close_theta(rtmp, ztmp, thgrd)
                            j += 1
                        pfa = 0.
                        rba = data['fpol'][0]
                        self.pfm[jr, jz] = quad_int(norm_grd, norm, [pfa, pf[0], pf[1]])
# Ensure continuity with EQI at separatrix
                        self.pfm[jr, jz] += self.pf_shift 
                        rbp[jr, jz] = quad_int(norm_grd, norm, [rba, rb[0], rb[1]])
                    else:
                        if irho == n_xb - 1:
                            irho1 = n_xb - 3
                            irho2 = n_xb - 1
                        else:
                            irho1 = irho - 1
                            irho2 = irho + 1

                    j=0
                    for jrho in range(irho1, irho2 + 1):
                        rtmp = r_plot[jrho, th1:th2] - R0
                        ztmp = z_plot[jrho, th1:th2] - z0
                        norm[j] = close_theta(rtmp, ztmp, thgrd)
                        j += 1
                    self.pfm[jr, jz] = quad_int(norm_grd, norm, pf[irho1:irho2+1])
                    self.pfm[jr, jz] += pf_sep[jt_aug] - pf[-1] 
                    rbp[jr, jz] = quad_int(norm_grd, norm, rb[irho1:irho2+1])

#=================
# B_R, B_z, B_tor
#=================

        b_z = np.zeros((nr, nz))
        b_r = np.zeros((nr, nz))
        b_p = np.zeros((nr, nz))
        b_t = np.zeros((nr, nz))
        dz  = np.gradient(zgrid)
        dR  = np.gradient(Rgrid)

        for jr in range(nr):
            tmpz = self.pfm[jr, :]
            der = np.gradient(tmpz)/dz
            der *= 50./np.pi
            b_z[jr, :] = 100.*der/Rgrid[jr]
        for jz in range(nz):
            tmpR = self.pfm[:, jz]
            der = np.gradient(tmpR)/dR
            der *= 50./np.pi
            b_r[:, jz] = - 100.*der/Rgrid
            b_t[:, jz] = - 100.*rbp[:, jz]/Rgrid
        b_z[:, 0] = b_z[:, 1]
        b_z[:,-1] = b_z[:,-2]
        b_r[0, :] = b_r[1, :]
        b_r[-1,:] = b_r[-2,:]
        b_p = np.hypot(b_z, b_r)

#=======
# Plots
#=======

        print('%8.4f' %cv['TIME3'][it])

# View

        if vplot:
            psi = self.pfm[:, :].T
            Rax = data['RAXIS']
            zax = data['YAXIS']
            X, Y = np.meshgrid(Rgrid, zgrid)

            figtitle = 'Pol. flux matrix  time_SF: %8.4f, time_CDF: %8.4f' \
                       %(eq['time'][jt_aug], cv['TIME3'][it])

            plt.ion()

            fig11 = plt.figure(11, figsize=(16,11))
            fig11.clf()
            fig11.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
            fig11.canvas.set_window_title('Psi contours') 

            ax1 = fig11.add_subplot(1, 2, 1, aspect='equal')
            levels = data['PLFLX2PI'][::-1]
            levels += pf_sep[jt_aug] - pf[-1] 
            ctr1 = ax1.contour(X, Y, psi, levels)
            fig11.colorbar(ctr1, aspect=10, shrink=0.9)
            ax1.plot(Rax, zax, 'go', markersize=4)
            ax1.scatter(X, Y, marker='o', color='blue', s=2)
            ax1.set_title(figtitle)
            for jrho in range(n_xb):
                ax1.plot(r_plot[jrho, :], z_plot[jrho, :], 'g-')

            ax2 = fig11.add_subplot(1, 2, 2, aspect='equal')
            levels = np.linspace(np.min(psi), np.max(psi), 20)
            ctr2 = ax2.contourf(X, Y, psi, levels)
            fig11.colorbar(ctr2, aspect=10, shrink=0.9)
            ax2.scatter(X, Y, marker='o', color='blue', s=2)
            ax2.set_title(figtitle)
            for jrho in range(n_xb):
                ax2.plot(r_plot[jrho, :], z_plot[jrho, :],'g-')

# Fields
            jzplot = nz/2 + 6

            fig12 = plt.figure(12, figsize=(19,11))
            fig12.clf()
            fig12.canvas.set_window_title('B-fields') 
            fig12.text( .37, .96, 'Runid %s  Time=%8.4f [s]  at z = %8.4f [m]' \
                       %(runid, cv['TIME3'][it], 0.01*zgrid[jzplot]) )

            fig12.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.93)

            ax21 = fig12.add_subplot(2, 2 ,1)
            ax21.plot(Rgrid, b_r[:, jzplot])
            ax21.set_xlabel('R [m]')
            ax21.set_ylabel('B_R [T]')

            ax22 = fig12.add_subplot(2, 2, 2)
            ax22.plot(Rgrid, b_z[:, jzplot])
            ax22.set_xlabel('R [m]')
            ax22.set_ylabel('B_z [T]')

            ax23 = fig12.add_subplot(2, 2, 3)
            ax23.plot(Rgrid, b_p[:, jzplot])
            ax23.set_xlabel('R [m]')
            ax23.set_ylabel('B_pol [T]')

            ax24 = fig12.add_subplot(2, 2, 4)
            ax24.plot(Rgrid, b_t[:,jzplot])
            ax24.set_xlabel('R [m]')
            ax24.set_ylabel('B_tor [T]')

            plt.draw()

        cdf.close()

