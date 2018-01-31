def pol_sect(ax):

    Rlbl = 'R [cm]'
    zlbl = 'z [cm]'
    fsize = 12
    ax.set_xlabel(Rlbl, fontsize=fsize)
    ax.set_ylabel(zlbl, fontsize=fsize)
    try:
        import plot_aug
        plot_aug.vessel_pol(ax, fac=100.)
        xpol = (90, 230)
        ypol = (-125, 125)
        ax.set_xlim(xpol)
        ax.set_ylim(ypol)
    except:
        pass


def tor_sect(ax, rsurf):

    import numpy as np

    ntheta = 101
    Rmaj = rsurf[0, 0]
    Rtor_in  = np.min(rsurf)
    Rtor_out = np.max(rsurf)
    phi_tor = np.linspace(0, 2*np.pi, ntheta)
    cosp = np.cos(phi_tor)
    sinp = np.sin(phi_tor)

    ax.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
    ax.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
    ax.plot(Rmaj*cosp, Rmaj*sinp, 'r--')
    try:
        import plot_aug
        plot_aug.vessel_tor(ax, fac=100.)
        ax.set_xlim((-600, 400))
    except:
        pass

# Overlay MC grid cells
def pol_cells(ax, fbm):

    if fbm is None:
        return
    for irho in range(fbm.rsurf.shape[0]):
        ax.plot(fbm.rsurf[irho, :], fbm.zsurf[irho, :], 'r-')
    for jbar, myr in enumerate(fbm.rbar):
        ax.plot(myr, fbm.zbar[jbar], 'r-')
