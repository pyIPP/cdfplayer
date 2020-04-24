# Used in cdfplayer

import os, shutil
from scipy.io import netcdf
import numpy as np
import ww_20180130
import sfh_20180130
import ctr2rz, get_sf_grid, tr_read_ctr

ww = ww_20180130.shotfile()
sfh = sfh_20180130.SFH()


def cdf2tre(runid, time=-1):

    homdir = os.getenv('HOME')
    os.system('mkdir -p %s/shotfiles/TRE' %homdir)
    source = '/afs/ipp/home/t/transp/pub/TRE00000.sfh.temp'
    sfhdir = '%s/tr_client/AUGD' %homdir
    nshot = int(runid[0:5])
    tail = runid[5:8]
    fcdf = '%s/%d/%s/%s.CDF' %(sfhdir, nshot, tail, runid)
    fsfh = '%s/tr_client/AUGD/TRE00000.sfh' %homdir
    if not os.path.isfile(fcdf):
        print('%s not found' %fcdf)
        return
    try:
        shutil.copy2(source, fsfh)
    except:
        print('Unable to copy file %s to %s' %(source, fsfh))
        return
    cv = netcdf.netcdf_file(fcdf, 'r', mmap=False).variables

# Parameters

    os.chdir(sfhdir)

    coords = ('RMAJM', 'XB', 'X')
    cdf2d = ['DVOL', 'DAREA', 'ELONG', 'TRFLX', 'BPOL', 'PLFLX', \
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

    nt = len(cv['TIME3'].data)
    n_xb = len(grid['XB'])
    n_xb1 = n_xb + 1
    n_x  = len(grid['X'])
    n_ssq = len(ssqlist)

    if time == -1:
        time = cv['TIME3'][-1]

#   grid['XB'] = 0.025, 0.05, ...,1.0
#   grid['X']  = 0.0125,0.0375,...,0.9875

    gfun = cv['GFUN'].data.T
    Jpol = (-5.e4*gfun*cv['BZXR'].data).T
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
        xb_prof['PFL'] [jt, :] = np.append(0, -2*np.pi*cv['PLFLX'][jt, :])
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

    eq = get_sf_grid.get_grid(runid)
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
    tre['time'] = np.float32(cv['TIME3'].data)
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

    tre['IpiPSI'] = np.float32(cv['PCUR'].data.reshape(1, nt))
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

    fm = tr_read_ctr.TR_READ_CTR(fcdf, rho=1, nthe=nthe_eqd)
    r_sep = 0.01*fm.Rsurf
    z_sep = 0.01*fm.Zsurf
    j_xpoint = np.argmin(z_sep, axis=1)
    j_up     = np.argmax(z_sep, axis=1)

    tre['zPFx'][1, :] = np.min(z_sep, axis=1)  # z lower X-point
    ssq_d['Zoben'] = np.max(z_sep, axis=1)
    ssq_d['Rin']   = np.min(r_sep, axis=1)
    ssq_d['Raus']  = np.max(r_sep, axis=1)

#-----------------
# Bottle neck: PFM
#-----------------

    for jt in range(nt):
        rz = ctr2rz.CTR2RZ(fcdf, it=jt, vplot=False)
        tre['PFM'][:, :, jt] = rz.pfm           # Pol. flux matrix
        tre['PFL'][:, jt] += rz.pf_shift
        tre['RPFx'][1, jt] = r_sep[jt, j_xpoint[jt]]  # R lower X-point
        ssq_d['Rzoben'][jt] = r_sep[jt, j_up[jt]]
        for rho in (0, 0.25, 0.50, 0.75, 0.95):
            lbl = 'q%d' %int(100*rho)
            ssq_d[lbl][jt] = np.interp(rho, xb_grid**2, xb_prof['Qpsi'][jt])

    ssq_d['Rzunt'] = tre['RPFx'][1, :] 
    ssq_d['Zunt']  = tre['zPFx'][1, :]

# Magnetic axis
    tre['PFxx'][0, :] = tre['PFL'][-1, :]
    tre['RPFx'][0, :] = 0.01*cv['RAXIS'].data
    tre['zPFx'][0, :] = 0.01*cv['YAXIS'].data
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

    ssq_d['Rgeo']     = 0.01*(cv['RAXIS'].data - cv['ASHAF'].data)
    ssq_d['Rax-Rgeo'] = 0.01*cv['ASHAF'].data
    ssq_d['Zgeo']     = np.zeros(nt)
    ssq_d['Rmag']     = 0.01*cv['RAXIS'].data
    ssq_d['Zmag']     = 0.01*cv['YAXIS'].data
    ssq_d['betpol']   = cv['BETAT'].data
    ssq_d['li']       = cv['LIO2'].data
    ssq_d['Vol']      = psi_prof['Vol'][:, -1]
    ssq_d['k']        = cv['ELONG'][:, -1]

# Create SSQ array in right sequence for shotfile writing
    ssq = np.zeros((n_ssq, nt), dtype=np.float32)
    for jssq, lbl in enumerate(ssqlist):
        ssq[jssq, :] = ssq_d[lbl]

# Modify sfh

    sig_gr = sig3d + prof1 + prof2 + mixed + specs + switch

    err = sfh.Open(fsfh)
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
    else:
        print('Error=', err)

# Write TRE shotfile

    exp  = os.getenv('USER')

    diag = 'TRE'
    if ww.Open(exp, diag, nshot):
        for lbl in sig1d + sig_gr:
            status = ww.SetObject(lbl, tre[lbl])

# Parameter Sets
        ps = 'PARMV'
        for pn in tre[ps].keys():
            val = tre[ps][pn]
            print('Writing %s  value = %s' %(pn, val))
            status = ww.SetParameter(ps, pn, val)

        status = ww.SetParameter('RunDocu', 'runid', runid)
        try:
            status = ww.SetParameter('RunDocu', 'CLISTexp', eq['!exp'])
        except:
            print('No CLISTE-exp found to store in ParameterSet RunDocu')
        try:
            status = ww.SetParameter('RunDocu', 'CLISTdia', eq['!dia'])
        except:
            print('No CLISTE-diag found to store in ParameterSet RunDocu')
        try:
            cl_ed = np.int32(eq['!ed'])
            status = ww.SetParameter('RunDocu', 'CLISTed', cl_ed)
        except:
            print('No CLISTE-edition found to store in ParameterSet RunDocu')

# Write SSQ, SSQnam
        status = ww.SetObject('SSQ', ssq)
        status = ww.SetObject('SSQnam', ssqlist)

        ww.Close()


def cdf2tra(runid):

    homdir = os.getenv('HOME')
    os.system('mkdir -p %s/shotfiles/TRA' %homdir)
    source = '/afs/ipp/home/t/transp/pub/TRA00000.sfh.temp'
    sfhdir = '%s/tr_client/AUGD' %homdir
    nshot = int(runid[0:5])
    tail = runid[5:8]
    fcdf = '%s/%d/%s/%s.CDF' %(sfhdir, nshot, tail, runid)
    fnml = '%sTR.DAT' %fcdf[: -4]
    fsfh = '%s/TRA00000.sfh' %sfhdir
    for fname in fcdf, fnml:
        if not os.path.isfile(fname):
            print('%s not found' %fname)
            return
    try:
        shutil.copy2(source, fsfh)
    except:
        print('Unable to copy file %s to %s' %(source, fsfh))
        return
    sfh_dic = sfh_20180130.sfhmod(fcdf, nml=fnml, fsfh=fsfh)
    ww_20180130.write_sf(nshot, sfh_dic, sfhdir, 'TRA', exp=os.getenv('USER'))


def cdf2nub(runid):

    homdir = os.getenv('HOME')
    os.system('mkdir -p %s/shotfiles/NUB' %homdir)
    source = '/afs/ipp/home/n/nubeam/pub/NUB00000.sfh.temp'
    sfhdir = '%s/nb_client/AUGD' %homdir
    nshot = int(runid[0:5])
    tail = runid[5:8]
    fcdf = '%s/%d/%s/%s.cdf' %(sfhdir, nshot, tail, runid)
    fnml = '%s/%d/%s.nml'    %(sfhdir, nshot, runid)
    fsfh = '%s/NUB00000.sfh' %sfhdir
    for fname in fcdf, fnml:
        if not os.path.isfile(fname):
            print('%s not found' %fname)
            return
    try:
        shutil.copy2(source, fsfh)
    except:
        print('Unable to copy file %s to %s' %(source, fsfh))
        return
    sfh_dic = sfh_20180130.sfhmod(fcdf, nml=fnml, fsfh=fsfh)
    ww_20180130.write_sf(nshot, sfh_dic, sfhdir, 'NUB', exp=os.getenv('USER'))


if __name__ == "__main__":

    import sys
    if len(sys.argv) > 1:
        runid = sys.argv[1]
        cdftra(runid)
