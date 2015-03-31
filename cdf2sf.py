# Used in cdfplayer

import sys, os, shutil
from netCDF4 import Dataset
import matplotlib.pylab as plt
import numpy as np
import ww_20150119 as ww
import sfh_20140403 as sfh
import rz2psi, read_eq_sf, mom2rz, read_equ

ww = ww.shotfile()
sfh = sfh.SFH()


class CDF2TRE:


    def __init__(self, cdf_file, vplot=True, time=-1, sfh_tre='TRE00000.sfh'):

        cdf = Dataset(cdf_file, 'r', format='NETCDF4')
        cv = cdf.variables

# Variables (user defined)

        runid = cdf_file.split('/')[-1][:8]

# Parameters
        if not os.path.isfile(sfh_tre):
            print('%s not found' %sfh_tre)
            return

        sfhdir = sfh_tre[:-13]
        os.chdir(sfhdir)

        coords = ('RMAJM', 'XB', 'X')
        cdf2d = ['DVOL', 'DAREA', 'ELONG', 'TRFLX', 'BPOL', 'PLFLX2PI', \
                 'Q', 'PMHD_IN', 'PPLAS', 'CUR', 'PLCURTOT', 'GFUN', 'PLFMP']
        cdf1d = ['ASHAF', 'BETAT', 'LIO2', 'PCUR', 'BZ', 'BZXR', 'RAXIS', 'YAXIS']

        prof1  = ('Ri', 'Zj', 'Lpf', 'IpiPSI')
        prof2  = ('PFL', 'TFLx', 'Qpsi', 'FFP')
        specs  = ('LPFx', 'PFxx', 'RPFx', 'zPFx')
        mixed  = ('Vol', 'Area', 'Pres', 'Jpol')
        sig3d  = ('PFM',) #'CDM')
        sig1d  = ('ixti', 'time')
        switch = ('ikCAT',)

        ssqlist = ('Rsquad', 'zsquad', 'Rgeo', 'Zgeo', 'Rxpu', 'Zxpu', 'ahor', 'k', \
                   'Rin', 'Raus', 'Zunt', 'Rzunt', 'Zoben', 'Rzoben', 'Rmag', 'Zmag', \
                   'Vol', 'IVSF', 'Slunt', 'Srunt', 'fXP1fIL', 'fXP1fAL', 'XPfdif', \
                   'lenH-1', 'lenH-2' ,'lenH-3', 'Rxpo', 'Zxpo', 'Slobn', 'Srobn', \
                   'bpli2', 'betpol', 'li', 'Wmhd', 'Circumf', 'bver', 'q95/q0', \
                   'Finstab', 'Fstab', 'Rax-Rgeo', 'fax-bnd', 'fbnd-ref', \
                   'q0', 'q95', 'fbnd-f12', 'fbnd-f17', 'fbnd-f22', 'fbnd-f25', \
                   'q25', 'q50', 'q75', 'Zgeri2b', 'Zgera2b', 'Zskewi2b', 'Zskewa2b', \
                   'Rrrofi2b', 'Rroofa2b', 'lenH-4', 'delRoben', 'delRuntn', \
                   'GapAbot', 'GamAmid', 'GapAtop', 'GapAmin', \
                   'lenH-5', 'lenV-1', 'koben', 'kuntn', 'lenV-2', 'lenV-3', 'dRXP')

#==========
# Profiles
#==========

        grid = {}
        for coord in coords:
            grid[coord] = cv[coord][0, :]

        nt = len(cv['TIME3'][:])
        n_xb = len(grid['XB'])
        n_xb1 = n_xb + 1
        n_x  = len(grid['X'])
        n_ssq = len(ssqlist)

        if time == -1:
            time = cv['TIME3'][-1]

#   grid['XB'] = 0.025, 0.05, ...,1.0
#   grid['X']  = 0.0125,0.0375,...,0.9875

        gfun = cv['GFUN'][:].T
        Jpol = (-5.e4*gfun*cv['BZXR']).T
        xb_prof = {}
        psi_prof = {}
        ddpsi = {}

        for prof in mixed + prof2:
            xb_prof [prof] = np.zeros((nt, n_xb1), dtype=np.float32)
            psi_prof[prof] = np.zeros((nt, n_xb1), dtype=np.float32)
            ddpsi[prof]    = np.zeros((nt, n_xb1), dtype=np.float32)

        xb_prof['Vol'] [:, 1:n_xb1] = 1.e-6*np.cumsum(cv['DVOL'] [:, :], axis=1)
        xb_prof['Area'][:, 1:n_xb1] = 1.e-4*np.cumsum(cv['DAREA'][:, :], axis=1)

# Interpolate on X and XB grids, to have profiles and derivatives on XB grid

        xb_grid = np.append(0, grid['XB'])
        for jt in range(nt):
            xb_prof['PFL'] [jt, :] = np.append(0, -cv['PLFLX2PI'][jt, :])
            xb_prof['TFLx'][jt, :] = np.append(0, -cv['TRFLX']   [jt, :])
            xb_prof['Qpsi'][jt, :] = np.interp(xb_grid, grid['XB'], cv['Q'][jt, :])
            xb_prof['Jpol'][jt, :] = np.interp(xb_grid, grid['XB'], Jpol[jt, :])
            xb_prof['Pres'][jt, :] = np.interp(xb_grid, grid['X'] , cv['PMHD_IN'][jt, :])

# Regular psi grid

        for jt in range(nt):
            pfl = xb_prof['PFL'][jt, :]
            psi_prof['PFL'][jt, :] = np.linspace(pfl[0], pfl[-1], n_xb1)
        dpsi = psi_prof['PFL'][:, 1] - psi_prof['PFL'][:, 0] # time dependent
        for lbl in  ('TFLx', 'Qpsi') + mixed:
            print 'TEST PFL monotonic?', psi_prof['PFL'][0, :]
            for jt in range(nt):
                tmp = np.interp(psi_prof['PFL'][jt, ::-1], xb_prof['PFL'][jt, ::-1], xb_prof[lbl][jt, ::-1])
                psi_prof[lbl][jt, :] = tmp[::-1]
# Derivatives d/dPsi
                if lbl in mixed:
                    ddpsi[lbl][jt, :] = np.gradient(psi_prof[lbl][jt, :])/dpsi[jt]
        psi_prof['FFP'] = 4.e-14*psi_prof['Jpol']*ddpsi['Jpol']

#==================
# Psi(R, z), j(R, z)
#==================

        eq = read_eq_sf.get_grid(runid)
        Rgrid = eq['Ri'][:, 0]
        zgrid = eq['Zj'][:, 0]
        nr = len(Rgrid)
        nz = len(zgrid)

#================
# shotfile output
#================

        nshot = int(runid[:5])

        tre = {}
        tre['ixti'] = np.arange(nt) + 1
        tre['time'] = np.float32(cv['TIME3'])
# R, z grids
        tre['Ri'] = np.repeat(np.float32(Rgrid), nt).reshape(nr, nt)
        tre['Zj'] = np.repeat(np.float32(zgrid), nt).reshape(nz, nt)
# Lpf profiles, from plasma edge to magnetic axis
        tre['Lpf']  = np.repeat( np.int32(n_xb), nt ).reshape(1, nt)

        for lbl in prof2:
            tre[lbl] = psi_prof[lbl].T   # (t, r) -> (r, t)
            tre[lbl] = tre[lbl][::-1, :] # [0]<-> separatrix

# Profiles with d/dpsi mixed
        for lbl in mixed:
            tre[lbl] = np.zeros((2*n_xb1, nt), dtype=np.float32)
            for jt in range(nt):
                tre[lbl][:, jt] = np.array([psi_prof[lbl][jt, ::-1], ddpsi[lbl][jt, ::-1]]).T.ravel()

        tre['IpiPSI'] = np.float32(cv['PCUR'][:].reshape(1, nt))
        tre['ikCAT']  = np.repeat( np.int32(4), nt ).reshape(1, nt)

# Matrices and separatrix
        tre['PFM'] = np.zeros((nr, nz, nt), dtype=np.float32)

        nthe_eqd = 101

# Special points: axis, sepx, lim=sepx
        lpfx = 5
        tre['LPFx'] = np.repeat( np.int32(lpfx-1), nt ).reshape(1, nt)
        tre['PFxx'] = np.zeros((lpfx, nt), dtype=np.float32)
        tre['RPFx'] = np.zeros((lpfx, nt), dtype=np.float32)
        tre['zPFx'] = np.zeros((lpfx, nt), dtype=np.float32)

# SSQ

        ssq_d = {}
        for lbl in ssqlist:
            ssq_d[lbl] = np.zeros(nt)

        fm = read_equ.READ_EQU(cdf_file)
        r_plot, z_plot = mom2rz.mom2rz(fm.rc, fm.rs, fm.zc, fm.zs, nthe=nthe_eqd)
        r_sep = 0.01*r_plot[:, -1, :]
        z_sep = 0.01*z_plot[:, -1, :]
        j_xpoint = np.argmin(z_sep, axis=1)
        j_up     = np.argmax(z_sep, axis=1)

        tre['zPFx'][1, :] = np.min(z_sep, axis=1)  # z lower X-point
        ssq_d['Zoben'] = np.max(z_sep, axis=1)
        ssq_d['Rin']   = np.min(r_sep, axis=1)
        ssq_d['Raus']  = np.max(r_sep, axis=1)

#==============
# Profile plots
#==============

        jtplot = np.argmin(np.abs(cv['TIME3'] - time))
        print('TIME:', jtplot, cv['TIME3'][jtplot], time)

        tplot = [jtplot]

# Fig 1: profiles

        n_col =3
        n_row = (len(cdf2d) + 3)/n_col
        if vplot:
            plt.ion()
            plt.figure(101, figsize=(6*n_col, 2.2*n_row))
            plt.subplots_adjust(left=0.07, bottom=0.05, right=0.98, top=0.97, wspace=0.15, hspace=0.3)

        for jt in tplot:
            if vplot:
                plt.clf()
            j_plot=1
            for lbl in cdf2d:
                drho  = cv[lbl].dimensions[1]
                dunit = cv[lbl].units
                dname = cv[lbl].long_name
                if vplot:
                    plt.subplot(n_row, n_col, j_plot)
                    plt.xlabel(drho, fontsize=12)
                    if lbl == 'GFUN':
                        plt.plot(grid[drho], Jpol[jt])
                        plt.ylabel('fpol TESLA * M', fontsize=12)
                    else:
                        plt.plot(grid[drho], cv[lbl][jt])
                        plt.ylabel(lbl.strip()+' '+dunit.strip(), fontsize=12)
                    j_plot += 1

            if vplot:
                plt.subplot(n_row, n_col, j_plot)
                plt.xlabel('PFL', fontsize=12)
                plt.ylabel('VOLUME M**3', fontsize=12)
                plt.plot(psi_prof['PFL'][jt, :], psi_prof['Vol'][jt, :])
                j_plot += 1
                plt.subplot(n_row, n_col, j_plot)
                plt.xlabel('PFL', fontsize=12)
                plt.ylabel('AREA M**2', fontsize=12)
                plt.plot(psi_prof['PFL'][jt, :], psi_prof['Area'][jt, :])
                plt.draw()

# Fig 2 (derivatives)

        n_row = len(mixed)
        plt.figure(102, figsize=(7, 1.8*n_row))
        plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.95)
        for jt in tplot:
            plt.clf()
            j_plot=1
            for lbl in mixed:
                plt.subplot(n_row, 1, j_plot)
                plt.plot(psi_prof['PFL'][jt, :], ddpsi[lbl][jt, :])
                plt.xlabel('Pol.flux', fontsize=12)
                plt.ylabel('d%s/dPsi' %lbl.strip(), fontsize=12)
                j_plot += 1
            plt.draw()

#-----------------
# Bottle neck: PFM
#-----------------

        for jt in range(nt):
            rz = rz2psi.RZ2PSI(cdf_file, it=jt)
            tre['PFM'][:, :, jt] = rz.pfm           # Pol. flux matrix
            tre['PFL'][:, jt] += rz.pf_shift
            tre['RPFx'][1, jt] = r_sep[jt, j_xpoint[jt]]  # R lower X-point
            ssq_d['Rzoben'][jt] = r_sep[jt, j_up[jt]]
            for rho in (0, 0.25, 0.50, 0.75, 0.95):
                lbl = 'q%d' %int(100*rho)
                ssq_d[lbl][jt] = np.interp(rho, xb_grid**2, xb_prof['Qpsi'][jt])

        ssq_d['Rzunt'] = tre['RPFx'][1, :] 
        ssq_d['Zunt']  = tre['zPFx'][1, :]
        print 'Test PFL ', tre['PFL'][:, -1]

# Magnetic axis
        tre['PFxx'][0, :] = tre['PFL'][-1, :]
        tre['RPFx'][0, :] = 0.01*cv['RAXIS'][:]
        tre['zPFx'][0, :] = 0.01*cv['YAXIS'][:]
# X point, 2nd X-point
        tre['PFxx'][1, :] = tre['PFL'][0, :]
        tre['PFxx'][3, :] = tre['PFL'][0, :]
        tre['PFxx'][4, :] = tre['PFL'][0, :]

# Parameters
        tre['PARMV'] = {}
        tre['PARMV']['M']     = np.int32(nr-1)
        tre['PARMV']['N']     = np.int32(nz-1)
        tre['PARMV']['NCL']   = np.int32(13)
        tre['PARMV']['NTIME'] = np.int32(nt)

        ssq_d['Rgeo']     = 0.01*(cv['RAXIS'][:] - cv['ASHAF'][:])
        ssq_d['Rax-Rgeo'] = 0.01*cv['ASHAF'][:]
        ssq_d['Zgeo']     = np.zeros(nt)
        ssq_d['Rmag']     = 0.01*cv['RAXIS'][:]
        ssq_d['Zmag']     = 0.01*cv['YAXIS'][:]
        ssq_d['betpol']   = cv['BETAT'][:]
        ssq_d['li']       = cv['LIO2'][:]
        ssq_d['Vol']      = psi_prof['Vol'][:, -1]
        ssq_d['k']        = cv['ELONG'][:, -1]

# Create SSQ array in right sequence for shotfile writing
        ssql = []
        for lbl in ssqlist:
            ssql.append(ssq_d[lbl])
        ssq = np.array(ssql, dtype=np.float32)

# Modify sfh

        sig_gr = sig3d + prof1 + prof2 + mixed + specs + switch

        err = sfh.Open(sfh_tre)
        print 'Error=', err
        if err == 0:
            for lbl in sig1d:
                print('SFHmod %s' %lbl)
                sfh.Modtim(lbl, nt)
            for lbl in sig_gr:
                print('SFHmod %s' %lbl)
                dims = tre[lbl].shape
                sfh.Modsgr(lbl, dims)
            sfh.Modsgr('SSQ', (n_ssq, nt))
            sfh.Close()

# Write TRE shotfile

        exp  = os.getenv('USER')

        diag = 'TRE'
        if ww.Open(exp, diag, nshot):
            for lbl in sig1d:
                print('Writing %s' %lbl)
                status = ww.SetSignal(lbl, tre[lbl])
            for lbl in sig_gr:
                print 'Writing ' + lbl, tre[lbl].shape
                status = ww.SetSignalGroup(lbl, tre[lbl])

# Parameter Sets
            ps = 'PARMV'
            for pn in tre[ps].iterkeys():
                val = tre[ps][pn]
                print 'Writing ' + pn, ' value=', val
                status = ww.SetParameter(ps, pn, val)

            status = ww.SetParameter('SSQnames', 'SSQn', ssqlist, typ=6, cjust=8)

            status = ww.SetParameter('RunDocu', 'runid', runid, typ=6, cjust=8)
            try:
                status = ww.SetParameter('RunDocu', 'CLISTexp', eq['!exp'], typ=6, cjust=8)
            except:
                print('No CLISTE-exp found to store in ParameterSet RunDocu')
            try:
                status = ww.SetParameter('RunDocu', 'CLISTdia', eq['!dia'], typ=6, cjust=8)
            except:
                print('No CLISTE-diag found to store in ParameterSet RunDocu')
            try:
                cl_ed = np.int32(eq['!ed'])
                status = ww.SetParameter('RunDocu', 'CLISTed', cl_ed, typ=1)
            except:
                print('No CLISTE-edition found to store in ParameterSet RunDocu')

            
# Write SSQ, SSQnam
            status = ww.SetSignalGroup('SSQ', ssq)

            ww.Close()


class CDF2TRA:


    def __init__(self, runid):

        homdir = os.getenv('HOME')
        os.system('mkdir -p %s/%s' %(homdir, 'shotfiles/TRA'))
        source = '/afs/ipp/home/t/transp/pub/TRA00000.sfh.temp'
        sfhdir = '%s/tr_client/AUGD' %homdir
        nshot = int(runid[0:5])
        tail = runid[5:8]
        cdf_file = '%s/%d/%s/%s.CDF' %(sfhdir, nshot, tail, runid)
        fsfh = '%s/TRA00000.sfh' %sfhdir
        nml = '%sTR.DAT' %cdf_file[: -4]
        if not os.path.isfile(cdf_file):
            print ('%s not found' %cdf_file)
            return
        if not os.path.isfile(nml):
            print ('%s not found' %nml)
            return
        try:
            shutil.copy2(source, fsfh)
        except IOError, e:
            print('Unable to copy file %s\n%s' %(source, e))
        sfh_dic = sfh.sfhmod(cdf_file, nml=nml, fsfh=fsfh)
        ww.write_sf(nshot, sfh_dic, sfhdir, diag='TRA')


if __name__ == "__main__":

    if len(sys.argv) > 1:
        runid = sys.argv[1]
    CDF2TRA(runid)
