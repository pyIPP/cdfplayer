import numpy as np


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


def interp_norm(r_in, z_in, th_ref):

    nthe = len(r_in)
    theta = np.arctan2(z_in, r_in) - th_ref

# Discontinuity at 2pi -> 0
    jmid = np.argmin(np.abs(theta))
    j1 = jmid - 1
    j2 = jmid + 1
    if jmid == 0:
        j1 = nthe - 1
        theta[j1] += -2.*np.pi
    elif jmid == nthe-1:
        j2 = 0
        theta[j2] += 2.*np.pi

    ind3 = [j1, jmid, j2]
    theta3 = theta[ind3]
    norm3 = np.hypot(r_in[ind3], z_in[ind3])

    return quad_int(0., theta3, norm3)


def map_norm(Rrect, Zrect, Rctr, Zctr):

    nr = len(Rrect)
    nz = len(Zrect)
    n_xb = Rctr.shape[0] - 1
    iclose = np.zeros((nr, nz), dtype=np.int32)
    for jr in range(nr):
        rdist = (Rctr - Rrect[jr])**2
        for jz in range(nz):
            zdist = (Zctr - Zrect[jz])**2
            jmin = np.argmin(rdist + zdist)
            iclose[jr, jz], _ = np.unravel_index(jmin, Rctr.shape, order='C')

    return iclose
