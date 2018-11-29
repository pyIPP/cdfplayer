import os
from scipy.io import netcdf
import numpy as np


class READ_FBM:


    def __init__(self, f_fbm):


        if not os.path.isfile(f_fbm):
            print('File %s not found' %f_fbm)
            return

        cv = netcdf.netcdf_file(f_fbm, 'r', mmap=False).variables

#-------------------------------
# Read space-time grids from cdf
#-------------------------------

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

        n_cells = len(bmvol)
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

#-----------------------------------------------
# Read distribution, energy, pitch data from cdf
#-----------------------------------------------

        spec_lab = [] 
        for j in range(10):
            lbl = 'SPECIES_%d' %j
            if lbl in cv.keys():
                slbl = "".join(cv[lbl].data)
                spec_lab.append(slbl)

        self.fdist  = {}
        self.a_d    = {}
        self.e_d    = {}
        self.bdens  = {}
        self.btrap  = {}
        self.trap_pit   = {}
        self.int_en_pit = {}
        self.int_en_pit_frac_trap = {}
        self.dens_zone  = {}
        self.trap_zone  = {}
        self.dens_vol   = {}
        self.trap_vol   = {}

        for spc_lbl in spec_lab:

            self.fdist[spc_lbl] = 0.5*cv['F_%s' %spc_lbl].data
            print(cv['F_%s' %spc_lbl].units)
            self.a_d[spc_lbl]   =     cv['A_%s' %spc_lbl].data
            self.e_d[spc_lbl]   =     cv['E_%s' %spc_lbl].data
            fbm = self.fdist[spc_lbl]
            dE = np.diff(cv['EB_%s' %spc_lbl].data)
            dpa = np.diff(cv['AB_%s' %spc_lbl].data)

            n_cells, n_pit, n_E = fbm.shape

            dpa_dE = np.outer(dpa, dE)

# Trapped particles distribution

            fbm_trap = np.zeros((n_cells, n_pit, n_E))
            fbm_trap[:, :, :] = fbm[:, :, :]
            self.trap_pit[spc_lbl] = 1. - rmaj_min[rho_lab[:]]/self.r2d[:]
            for jcell in range(n_cells):
                ind_trap = (self.a_d[spc_lbl]**2 <= self.trap_pit[spc_lbl][jcell])
                fbm_trap[jcell, ~ind_trap, :] = 0

# Integrals

            self.int_en_pit[spc_lbl] = np.tensordot(fbm, dpa_dE, axes=((1, 2), (0, 1)))
            int_en_pit_trap = np.tensordot(fbm_trap, dpa_dE, axes=((1, 2), (0, 1)))

            self.dens_zone[spc_lbl] = np.zeros((n_zones, n_pit, n_E))
            self.trap_zone[spc_lbl] = np.zeros((n_zones, n_pit, n_E))

            self.dens_vol[spc_lbl]  = np.tensordot(fbm     , bmvol, axes=(0, 0))/np.sum(bmvol)
            self.trap_vol[spc_lbl]  = np.tensordot(fbm_trap, bmvol, axes=(0, 0))/np.sum(bmvol)
            vol_zone = np.zeros(n_zones)
            for jrho in range(n_zones):
                indrho = (rho_lab == jrho)
                vol_zone[jrho] = np.sum(bmvol[indrho])
                self.dens_zone[spc_lbl][jrho, :, :] = np.tensordot(fbm[     indrho, :, :], bmvol[indrho], axes = (0, 0))
                self.trap_zone[spc_lbl][jrho, :, :] = np.tensordot(fbm_trap[indrho, :, :], bmvol[indrho], axes = (0, 0))
                self.dens_zone[spc_lbl][jrho, :, :] *= 1/vol_zone[jrho]
                self.trap_zone[spc_lbl][jrho, :, :] *= 1/vol_zone[jrho]

            self.bdens[spc_lbl] = np.tensordot(self.dens_zone[spc_lbl], dpa_dE, axes=((1, 2), (0, 1)))
            self.btrap[spc_lbl] = np.tensordot(self.trap_zone[spc_lbl], dpa_dE, axes=((1, 2), (0, 1)))

            part_tot = np.sum(self.bdens[spc_lbl]*vol_zone)
            trap_tot = np.sum(self.btrap[spc_lbl]*vol_zone)
            print('Trapped #%12.4e    Total #%12.4e    Fraction %12.4e' %(trap_tot, part_tot, trap_tot/part_tot))
            print('Volume averaged fast ion density #12.4e m^-3' %(part_tot/vol))
            self.int_en_pit_frac_trap[spc_lbl] = int_en_pit_trap/self.int_en_pit[spc_lbl]
