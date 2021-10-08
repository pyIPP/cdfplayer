import os
import numpy as np
import aug_sfutils as sfu
import parse_nml

sfh = sfu.SFH()

def sfhmod(cdf_file, nml=None, fsfh='TRA00000.sfh'):

    from scipy.io import netcdf

    if not os.path.isfile(cdf_file):
        print('%s not found' %cdf_file)
        return
    fsfh = os.path.abspath(fsfh)
    print(fsfh)
    if not os.path.isfile(fsfh):
        print('%s not found' %fsfh)
        return

    runid = os.path.split(cdf_file)[1][:8]
    cv = netcdf.netcdf_file(cdf_file, 'r', mmap=False).variables
    sfh_d = {}

    sfo = sfu.SFREAD(sfh=fsfh)
    err = sfh.Open(fsfh)
    if err != 0:
        return

    tbases  = [ lbl for lbl in sfo.sfh['obj_nam'] if sfo.sfh[lbl].obj_type == 8 ]

    for obj in tbases:
        sfh_d[obj] = cv[obj].data
        nt_cdf = sfh_d[obj].shape[0]
        sfh.Modtim(obj, nt_cdf)

    for obj in sfo.parsets:
        if nml is None:
            print('Parameter sets not read from any namelist')
        else: # Read values from namelist
            parset = sfo.getparset(obj)
            sfh_d[obj] = {}
            for pn, val in parset.items():
                if pn == 'runid':
                    sfh_d[obj][pn] = [runid]
                else:
                    sfh_d[obj][pn] = parse_nml.parsenml(nml, pn, fmt=val.data_format)

    for obj in sfo.objects:

        name_flg = False
        if obj in cv.keys():
            name_flg = True
            cdfobj = obj
        else:
            if len(obj) == 8:
                name_flg = False
                for cdfobj in cv.keys():
                    if obj in cdfobj:
                        print('CDF signal %s stored as shotfile %s' %(cdfobj, obj))
                        name_flg = True
                        break

        if name_flg:
            sfh_d[obj] = cv[cdfobj].data
            nt_cdf = sfh_d[obj].shape[0]
            devtyp = sfo.sfh[obj].obj_type
            if sfh_d[obj].ndim == 2:
                nx_cdf = sfh_d[obj].shape[1]
            if devtyp == 13:
                if sfh_d[obj].ndim == 1:
                    nx_cdf = sfh_d[obj].shape[0]
                    sfh.Mdarea(obj, 1, nx_cdf, 0, 0)
                else:
                    sfh.Mdarea(obj, nt_cdf, nx_cdf, 0, 0)
            if devtyp == 6:
                sfh.Mdindex(obj, nx_cdf, 0, 0)
            if devtyp in (6, 7):
                sfh.Modtim(obj, nt_cdf)
        else:
            user_str = '/%s/' %os.getenv('USER')
            if user_str in cdf_file:
                cdff = '~/%s' %cdf_file.split(user_str)[1]
            else:
                cdff = cdf_file
            print('Signal %s not found in %s' %(obj, cdff))

    sfh.Close()

    return sfh_d
