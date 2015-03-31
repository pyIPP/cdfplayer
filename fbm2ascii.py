import os
import numpy as np
from netCDF4 import Dataset


class FBM2ASCII:


    def __init__(self, runid, t_id=1):

        sshot = runid[0:5]
        tail  = runid[5:]
        cdf_dir  = '%s/tr_client/AUGD/%s/%s' %(os.getenv('HOME'), sshot, tail)
        cdf_file = '%s/%s.CDF' %(cdf_dir, runid) 
        fbm_file = '%s/%s_fi_%d.cdf' %(cdf_dir, runid, int(t_id))

        print cdf_file
        print fbm_file

        fbm = Dataset(fbm_file, 'r', format='NETCDF4')
        cdf = Dataset(cdf_file, 'r', format='NETCDF4')
        cv = cdf.variables
        fv = fbm.variables

        for key, val in fv.iteritems():
            print key, val
        t_fbm = fv['TIME'][0]
        t_cdf = cv['TIME3'][:]
        it = np.argmin(np.abs(t_cdf-t_fbm))
        print('FBM time = %8.4f, CDF closest time = %8.4f' %(t_fbm, t_cdf[it]))

# CDF

        cdf_str  = ('R0 = %8.4f\n' %cv['RAXIS'][it])
        cdf_str += ('z0 = %8.4f\n' %cv['YAXIS'][it])
        cdf_str += ('Shafranov shift = %8.4f\n' %cv['ASHAF'][it])

        n_xb = len(cv['XB'][it])

        for lbl in ('XB', 'PLFLX', 'Q'):
            cdf_str += '\n%s%s\n' %(lbl, cv[lbl].units)
            for j_xb in range(n_xb):
                cdf_str += ('%12.4f' %cv[lbl][it, j_xb])
                if j_xb%6 == 5:
                    cdf_str += '\n'

        cdf_out = '%s/%scdf.txt' %(cdf_dir, runid)
        print('Writing %s' %cdf_out)
        f = open(cdf_out, 'w')
        f.write(cdf_str)
        f.close()

# Fourier moments

        fou_str = ''

        nmom = 0
        rc = 'RMC00'
        while rc in cv.iterkeys():
            nmom += 1
            rc = 'RMC0%d' %nmom
        print nmom

# Radial coordinate of Fourier moments is XB

        fou_out=''

        for lbl in ('XB', 'RMC00', 'YMC00'):
            fou_str += '\n\n'+lbl
            for j_xb in range(n_xb):
                fou_str += (' %8.4f' %cv[lbl][it, j_xb])
                if j_xb%8 == 7:
                    fou_str += '\n'
    
        for jmom in range(1, nmom):
            for lbl in ('RMC0', 'RMS0', 'YMC0', 'YMS0'):
                mom_str = '%s%d' %(lbl, jmom)
                fou_str += '\n\n%s' %mom_str
                for j_xb in range(n_xb):
                    fou_str += (' %8.4f' %cv[mom_str][it, j_xb])
                    if j_xb%8 == 7:
                        fou_str += '\n'

        fou_out = '%s/%sfou.txt' %(cdf_dir, runid)
        print('Writing %s' %fou_out)
        f = open(fou_out, 'w')
        f.write(fou_str)
        f.close()

# FBM

        n_cells, n_pit, n_E = fv['F_D_NBI'].shape
        print('#cells: %d    #p.a.: %d    #Energy: %d' %(n_cells, n_pit, n_E))

        fbm_str = ''
        for lbl in ('R2D', 'Z2D', 'X2D', 'TH2D'):
            fbm_str += lbl
            for jcell in range(n_cells):
                if jcell%8 == 0:
                    fbm_str += '\n'
                fbm_str += (' %8.4f' %fv[lbl][jcell])
            fbm_str += '\n'

        fbm_str += '\n\nDistribution function F_D_NBI\n\n'

        fbm_out = '%s/%sfbm.txt' %(cdf_dir, runid)

        f = open(fbm_out, 'w')
        f.write(fbm_str)

        for jE in range(n_E):
            print('jE = %d' %jE)
            for jpit in range(n_pit):
                f.write('\n')
                f.write('En   = %8.4f\n' %fv['E_D_NBI'][jE])
                f.write('p.a. = %8.4f\n' %fv['A_D_NBI'][jpit])
                for jcell in range(0, n_cells, 6):
                    dat = fv['F_D_NBI'][jcell:jcell+6, jpit, jE]
                    np.savetxt(f, dat, fmt=' %10.4e', newline='')
                    f.write('\n')

        cdf.close()
        fbm.close()
        print('Written %s' %fbm_out)
        f.close()


if __name__ == "__main__":

    runid = '23076E06'
    t_id  = 1
    FBM2ASCII(runid, t_id=t_id)

