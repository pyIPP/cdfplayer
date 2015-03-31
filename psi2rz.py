import numpy as np
import matplotlib.pylab as plt


class PSI2RZ:


    def plot_psi(self, plot_d, fill=True):

        psi = plot_d['psi'].T
        fig11 = plt.figure(11, figsize=(10.5, 12))
        plt.clf()

        nw = len(plot_d['Rgrid'])
        X, Y = np.meshgrid(plot_d['Rgrid'], plot_d['zgrid'])
        ax = fig11.add_subplot(1, 1, 1, aspect='equal')

#    levels=np.linspace(0, 1.6, 10)
        levels = np.linspace(np.min(psi), np.max(psi), 10)

        if fill:
            ctr = ax.contourf(X, Y, psi, levels)
        else:
            ctr = ax.contour(X, Y, psi, levels)
        fig11.colorbar(ctr, aspect=10, shrink=0.9)
        ax.scatter(X, Y, color='b', marker='o', s=2)
        if 'r_bdy' in plot_d.iterkeys():
            ax.plot(plot_d['r_bdy'], plot_d['z_bdy'], 'k-')
        for pos in ('axis', 'x'):
            Rc = 'R%s' %pos
            if Rc in plot_d.iterkeys():
                ax.plot(plot_d[Rc], plot_d['z%s' %pos], 'go')
        if 'eqdsk' in plot_d.iterkeys():
            str_title = 'PFM from eqdsk file %s' %plot_d['eqdsk']
        else:
            str_title = '# %d: PFM from shotfile %s %s,  edition %d at t = %6.3f' \
                        %( plot_d['shot'], plot_d['exp'], plot_d['diag'], \
                           plot_d['edit'], plot_d['time'] )
        ax.set_title(str_title)

# Profiles

        names = ('pres', 'fprim', 'pprim', 'qpsi')
        ncols = 2
        nrows = (len(names)+1)/ncols
        plt.figure(12, figsize = (6*ncols, 2.2*nrows))
#    plt.clf()
        plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.93, hspace=0.3)
        jplot = 0
        nx = len(plot_d['fprim'])
        xgrid = np.linspace(plot_d['psi_ax'], plot_d['psi_bdy'], nx)
        for lbl in names:
            plt.subplot(nrows, ncols, jplot)
            plt.title(lbl)
            plt.plot(xgrid, plot_d[lbl])
            jplot += 1

        plt.show()

if __name__ == "__main__":

    eqfile = '17870Q02.eqdsk'
    eqfile = 'demo.eqdsk'
    PSI2RZ(eqfile)
