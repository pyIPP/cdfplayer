import numpy as np

def mom2rz(rcos, rsin, zcos, zsin, nthe=101, endpoint=False):

# Fourier moments -> {R,z} of magnetic surfaces
# rcos [jt, jrho, jmom] or [jrho, jmom] or [jt, jmom] or [jmom]

    theta = np.linspace(0, 2*np.pi, nthe, endpoint=endpoint) 

    nmom  = np.size(rcos, -1)
    angle = np.outer(np.arange(nmom), theta)
    cos   = np.cos(angle)
    sin   = np.sin(angle)
    r_plot  = np.tensordot(rcos, cos, axes=([-1, 0]))
    r_plot += np.tensordot(rsin, sin, axes=([-1, 0]))
    z_plot  = np.tensordot(zcos, cos, axes=([-1, 0]))
    z_plot += np.tensordot(zsin, sin, axes=([-1, 0]))

    return r_plot, z_plot
