import numpy as np

class PHILOS:

    def __init__(self, par_d, npoints=5000):

        theta  = par_d['theta']
        phi_in = par_d['phi'] 
        x0 = par_d['x0']
        y0 = par_d['y0']
        z0 = par_d['z0']
        xend = par_d['xend']
        phi0 = np.arctan(x0/y0)
        phi = phi_in + phi0
        dx = xend - x0
        dr = dx/np.sin(phi)
        self.dr = np.abs(dr)
        yend = y0 + dr*np.cos(phi)
        zend = z0 + self.dr*np.tan(theta)
        self.xline = np.linspace(x0, xend, npoints)
        self.yline = np.linspace(y0, yend, npoints)
        self.zline = np.linspace(z0, zend, npoints)
        self.rline = np.hypot(self.xline, self.yline)


class XYLOS:

    def __init__(self, par_d, npoints=5000):

        theta  = par_d['theta']
        phi_in = par_d['phi'] 
        x0 = par_d['x0']
        y0 = par_d['y0']
        z0 = par_d['z0']
        xend = par_d['xend']
        yend = par_d['yend']
        dx = xend - x0
        dy = yend - y0
        self.dr = np.hypot(dx, dy)
        zend = z0 + self.dr*np.tan(theta)
        self.xline = np.linspace(x0, xend, npoints)
        self.yline = np.linspace(y0, yend, npoints)
        self.zline = np.linspace(z0, zend, npoints)
        self.rline = np.hypot(self.xline, self.yline)
