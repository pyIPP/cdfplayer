import os
from scipy.io import netcdf
import numpy as np


class READ_FBM:


    def __init__(self, f_fbm):


        if not os.path.isfile(f_fbm):
            print('File %s not found' %f_fbm)
            return

        cv = netcdf.netcdf_file(f_fbm, 'r', mmap=False).variables

        spc_lbl = "".join(cv['SPECIES_1'].data)
        print spc_lbl

# Read data from cdf

        self.fbm  = 0.5*cv['F_'+spc_lbl].data
        self.a_d  = cv['A_'+spc_lbl].data
        self.e_d  = cv['E_'+spc_lbl].data
        self.r2d  = cv['R2D'].data
        self.z2d  = cv['Z2D'].data
        self.x2d  = cv['X2D'].data
        self.th2d = cv['TH2D'].data
        self.time = cv['TIME'].data
        bmvol = cv['BMVOL'].data

        vol = 1.e-6*np.sum(bmvol)
        print('Volume is %8.4f m^-3' %vol)
        if cv['SYMFLAG'].data == 0:
            n_sym = 2 # asymmetric equilibrium
        else:
            n_sym = 1

        n_cells, n_pit, n_E = self.fbm.shape
        n_zones = cv['N_ZONEROWS'].data

        rho_lab = np.zeros(n_cells, dtype=int)      # rho index, takes values 0:n_zones-1

        self.rho_grid = np.unique(self.x2d)
        rmaj_min = np.zeros(n_zones)                # min(R(rho))

# Boundary grids for theta and rho (both are equispaced in TRANSP)

        max_jrho = 0
        thb_grid2d = []
        for jrho in range(n_zones):
            min_jrho = max_jrho
            max_jrho = min_jrho + 2*n_sym*(jrho + 1) # amount of poloidal cells at jrho
            rho_lab[min_jrho: max_jrho] = jrho
            rmaj_min[jrho] = min(self.r2d[min_jrho: max_jrho])
            the_loc = self.th2d[min_jrho: max_jrho]
            the_step = the_loc[1] - the_loc[0]
            thb_grid2d.append(the_loc - 0.5*the_step)

        rho_b = np.zeros(n_zones + 1)
        rho_step = self.rho_grid[1] - self.rho_grid[0]
        rho_b[:n_zones] = self.rho_grid     - 0.5*rho_step
        rho_b[n_zones]  = self.rho_grid[-1] + 0.5*rho_step

        dE = np.diff(cv['EB_%s' %spc_lbl].data)
        dpa = np.diff(cv['AB_%s' %spc_lbl].data)
        dpa_dE = np.outer(dpa, dE)

# Passing vs trapped

        fbm_pass = np.zeros((n_cells, n_E))
        fbm_trap = np.zeros((n_cells, n_E))
        self.trap_pit = 1. - rmaj_min[rho_lab[:]]/self.r2d[:]
        for jcell in range(n_cells):
            ind_trap = (self.a_d**2 <= self.trap_pit[jcell])
            if ind_trap.any():
                fbm_trap[jcell, :] = np.tensordot(self.fbm[jcell, ind_trap , :], dpa[ind_trap] , axes=(0, 0))
            if (~ind_trap).any():
                fbm_pass[jcell, :] = np.tensordot(self.fbm[jcell, ~ind_trap, :], dpa[~ind_trap], axes=(0, 0))

# Integrals

        self.int_en_pit = np.tensordot(self.fbm, dpa_dE, axes=((1, 2), (0, 1)))
        int_en_pass = np.tensordot(fbm_pass, dE, axes=(1, 0))
        int_en_trap = np.tensordot(fbm_trap, dE, axes=(1, 0))

        self.dens_zone = np.zeros((n_zones, n_pit, n_E))
        self.dens_vol  = np.tensordot(self.fbm, bmvol, axes=(0,0))/np.sum(bmvol)
        vol_zone = np.zeros(n_zones)
        for jrho in range(n_zones):
            indrho = (rho_lab == jrho)
            vol_zone[jrho] = np.sum(bmvol[indrho])
            self.dens_zone[jrho, :, :] = np.tensordot(self.fbm[indrho, :, :], bmvol[indrho], axes = (0, 0))
            self.dens_zone[jrho, :, :] *= 1/vol_zone[jrho]

        self.bdens = np.tensordot(self.dens_zone, dpa_dE, axes=((1, 2), (0, 1)))

        part_tot = np.sum(self.bdens*vol_zone)
        pass_tot = np.sum(int_en_pass*bmvol)
        trap_tot = np.sum(int_en_trap*bmvol)
        print('Passing #%12.4e  Trapped #%12.4e    Total #%12.4e' %(pass_tot, trap_tot, part_tot))
        print('Volume averaged fast ion density #12.4e m^-3' %(part_tot/vol))
        self.frac_trap = int_en_trap/(int_en_pass + int_en_trap)

# Rsurf, zsurf only of NUBEAM zones

        self.rsurf = []
        self.zsurf = []
        for rhob in rho_b:
            x_dist = (cv['XSURF'].data - rhob)**2
            irho = np.argmin(x_dist)
            self.rsurf.append(cv['RSURF'][irho, :])
            self.zsurf.append(cv['ZSURF'][irho, :])
        self.rsurf = np.array(self.rsurf)
        self.zsurf = np.array(self.zsurf)

# MC cells: grid bars

        self.rbar = []
        self.zbar = []
        for jrho in range(n_zones):
            for ythe in thb_grid2d[jrho]:
                ithe = np.min(np.where(cv['THSURF'].data > ythe))
                ind = [ithe - 1, ithe ]
                th_ref = cv['THSURF'][ind]
                rbar = np.zeros(2)
                zbar = np.zeros(2)
                rbar[0] = np.interp(ythe, th_ref, self.rsurf[jrho  , ind])
                zbar[0] = np.interp(ythe, th_ref, self.zsurf[jrho  , ind])
                rbar[1] = np.interp(ythe, th_ref, self.rsurf[jrho+1, ind])
                zbar[1] = np.interp(ythe, th_ref, self.zsurf[jrho+1, ind])
                self.rbar.append(rbar)
                self.zbar.append(zbar)
