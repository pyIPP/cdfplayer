import numpy as np
import datetime
import rw_for

now   = datetime.datetime.now()
w_for = rw_for.RW_FOR().wr_for
sspl  = rw_for.RW_FOR().ssplit


class EQDSK:


    def write_eqdsk(self, plot_d, f_eqdsk, time, eqlbl='TRANSP:', f_cdf=''):

        mixed = ('fpol', 'pres', 'fprim', 'pprim')

        rleft   = plot_d['Rgrid'][0]
        rdim    = plot_d['Rgrid'][-1]-rleft
        zdim    = plot_d['zgrid'][-1]-plot_d['zgrid'][0]
        rmaxis  = plot_d['Raxis']
        zmaxis  = plot_d['zaxis']
        simag   = plot_d['psi_ax']
        sibry   = plot_d['psi_bdy']
        current = plot_d['Ipl']
        bcentr  = plot_d['bcentr']
        rcentr  = plot_d['rcentr']
        psi_in  = plot_d['PFL']
        xdum  = 0.
        zmid  = 0.5*(plot_d['zgrid'][-1] + plot_d['zgrid'][0])
        nw    = len(plot_d['Rgrid'])
        nh    = len(plot_d['zgrid'])
        n_the = len(plot_d['r_bdy'])

        psi_grid = np.linspace(psi_in[0], psi_in[-1], nw)

        ind = np.argsort(psi_in)

# Interpolate onto an equidistant Psi-grid
        geq = {}
        for lbl in mixed + ('qpsi', ):
            geq[lbl] = np.interp(psi_grid, psi_in[ind], plot_d[lbl][ind])#[::-1]
        geq['fprim'] *= -2*np.pi

        fmt  = '%16.9E'
        nlin = 5
        format_str = (fmt * 5 + '\n')
        f = open(f_eqdsk,  'w')

# Scalars
        f.write( 'RZ2PSI %s  %s %s  t=%7.3f    0%4i%4i\n' \
               %(now.strftime('%d%b%Y'),  eqlbl,  f_cdf,  time,  nw,  nh) )
        f.write( format_str % (rdim, zdim, rcentr, rleft, zmid) )
        f.write( format_str % (rmaxis, zmaxis, simag, sibry, bcentr) )
        f.write( format_str % (current, simag, xdum, rmaxis, xdum) )
        f.write( format_str % (zmaxis, xdum, sibry, xdum, xdum) )

# Profiles
        for prof in mixed:
            f.write( w_for(geq[prof], fmt=fmt, n_lin=nlin) )
        f.write( w_for(plot_d['psi'], fmt=fmt, n_lin=nlin) ) # Matrix
        f.write( w_for(geq['qpsi']  , fmt=fmt, n_lin=nlin) )
        f.write( '%5i%5i\n' % (n_the, n_the) )
        tmp = np.append(plot_d['r_bdy'], plot_d['z_bdy'])
        rz  = tmp.reshape(2,  n_the)
        f.write( w_for(rz, fmt=fmt, n_lin=nlin) )
        f.write( w_for(rz, fmt=fmt, n_lin=nlin) )
        f.close()

        print('Output stored in %s' %f_eqdsk)


    def read_eqdsk(self, eqfile):

        f = open(eqfile, 'r')
        lines = f.readlines()

        header = lines[0].split()
        nw = int(header[-2])
        nh = int(header[-1])

        rdim, zdim, rcentr, rleft, zmid      = sspl(lines[1])
        rmaxis, zmaxis, simag, sibry, bcentr = sspl(lines[2])
        current, simag, xdum, rmaxis, xdum   = sspl(lines[3])
        zmaxis, xdum, sibry, xdum, xdum      = sspl(lines[4])

        jline = 5
        len_d = {'fpol': nw, 'fprim': nw, 'pres': nw, 'pprim': nw, \
                 'psi': nw*nh, 'qpsi':nw}

        plot_d = {}
        n_lin  = 5
        geq_sig = ('fpol', 'pres', 'fprim', 'pprim', 'psi', 'qpsi')
        for lbl in geq_sig:
          nx = len_d[lbl]
          plot_d[lbl] = np.zeros(nx)
          jl = 0
          while jl < nx:
            jr = min([jl+n_lin, nx])
            plot_d[lbl][jl: jr] = sspl(lines[jline])
            jline += 1
            jl += n_lin
        psi = plot_d['psi'].reshape(nh, nw).T

# Read bndy, lim

        nbbbs, limitr = lines[jline].split()
        n_bdy = int(nbbbs)
        n_lim = int(limitr)
        jline += 1
        len_d = {'bdy': n_bdy, 'lim': n_lim}
        for lbl in ('bdy', 'lim'):
          nx = 2*len_d[lbl]
          plot_d[lbl] = np.zeros(nx)
          jl = 0
          while jl < nx:
            jr = min([jl+5, nx])
            plot_d[lbl][jl:jr] = sspl(lines[jline])
            jline += 1
            jl += 5

        bdy = plot_d['bdy'].reshape(n_bdy, 2)
        rbbbs = bdy[0, :]
        zbbbs = bdy[1, :]
        lim = plot_d['lim'].reshape(n_lim, 2)

        plot_d['r_bdy']   = bdy[:, 0]
        plot_d['z_bdy']   = bdy[:, 1]
        plot_d['Raxis']   = rmaxis
        plot_d['zaxis']   = zmaxis
        plot_d['Rgrid']   = np.linspace(rleft, rleft + rdim, nw)
        plot_d['zgrid']   = np.linspace(zmid - 0.5*zdim, zmid + 0.5*zdim, nh)
        plot_d['psi']     = psi
        plot_d['psi_ax']  = simag
        plot_d['psi_bdy'] = sibry

        return plot_d
