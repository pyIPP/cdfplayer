import os, sys
sys.path.append('/afs/ipp/home/g/git/python/repository')
import numpy as np
from scipy.io import netcdf
import Tkinter as tk
import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import plot_sections, read_fbm
from tkMessageBox import *


def read_birth(birth_file, topframe=None):


    print('Reading %s' %birth_file)

    if not os.path.isfile(birth_file):
        showerror("Error", '%s not found' %birth_file)
        return

    birthfile  = birth_file.split('/')[-1]
    runid    = birthfile[:8]
    t_id     = birth_file[-1]

    fbm_file = birth_file.replace('birth.cdf%s' %t_id, 'fi_%s.cdf' %t_id)
    if not os.path.isfile(fbm_file):
        print('%s not found' %fbm_file)
        fbm = None
    else:
        fbm = read_fbm.READ_FBM(fbm_file)

    cv = netcdf.netcdf_file(birth_file, 'r', mmap=False).variables

    print('PLOT_BIRTH')
    print(birth_file)

    comp = ('full', 'half', '1/3')

    mcl = cv['mclabel'].data
    mc_label = "".join(mcl[0]).strip()

# Read data from cdf

# r, z : deposition location for each particle
# time
# einj = injection energy
# xksid = pitch angle
# zeta = toroidal angle
# ib = #NBI source associated to each MC marker

    Rj      = cv['bs_rgc_%s'   %mc_label].data
    zj      = cv['bs_zgc_%s'   %mc_label].data
    Einj    = cv['bs_einj_%s'  %mc_label].data
    pitch   = cv['bs_xksid_%s' %mc_label].data
    weight  = cv['bs_wght_%s'  %mc_label].data
    tor_ang = cv['bs_zeta_%s'  %mc_label].data
    t_birth = cv['bs_time_%s'  %mc_label][0]
    j_nbi   = cv['bs_ib_%s'    %mc_label].data

    n_birth = len(Rj)
    print('# of MC particles: %d' %n_birth)
    src_arr = np.unique(j_nbi)
    print('Sources: ', src_arr)

# Vessel compoments for plot

#    try:
    if True:
        import plot_aug
        xlin, ylin, rlin, zlin = plot_aug.nbi_plot(nbis=src_arr, runid=runid, raug=False)
#    except:
#        pass

    Efull = {}
    for jsrc in src_arr:
        index = np.where(j_nbi == jsrc )
        Efull[jsrc] = np.max(Einj[index])

# Determine 1/2 and 1/3 fractions

    mix_lbl = np.zeros(n_birth, dtype='|S4')
    mix_lbl[:] = comp[0]
    for jpart in range(n_birth):
        if Einj[jpart] <= 0.51*Efull[j_nbi[jpart]]:
            mix_lbl[jpart] = comp[1]
        if Einj[jpart] <= 0.34*Efull[j_nbi[jpart]]:
            mix_lbl[jpart] = comp[2]

    nr = 141
    nz = 101
    n_comp = 3
    Rmin = 100
    Rmax = 240
    zmin = -100
    zmax = 100
    R_grid = np.linspace(Rmin, Rmax, nr)
    z_grid = np.linspace(zmin, zmax, nz)
    deltaR = R_grid[1] - R_grid[0]
    deltaz = z_grid[1] - z_grid[0]
    dep_matrix={}
    for jsrc in src_arr:
        dep_matrix[jsrc] = {}
        for jcomp in comp:
            dep_matrix[jsrc][jcomp] = np.zeros((nr, nz))

    for jpart in range(n_birth):
        jsrc = j_nbi[jpart]
        jcomp = mix_lbl[jpart]
        Rpart = Rj[jpart]
        zpart = zj[jpart]
        if (Rpart > Rmin) and (Rpart < Rmax) and (zpart > zmin) and (zpart < zmax):
            jr = int((Rpart - Rmin)/deltaR)
            jz = int((zpart - zmin)/deltaz)
            dep_matrix[jsrc][jcomp][jr, jz] += weight[jpart]

    res_R = {}
    for jsrc in src_arr:
        res_R[jsrc] = {}
        for jcomp in comp:
            dep_R = np.sum(dep_matrix[jsrc][jcomp], axis=1) # z-sum
            res_R[jsrc][jcomp] = np.cumsum(dep_R)
            print('Deposited particles for source %d, component %s: %10.3e/s' %(jsrc, jcomp, res_R[jsrc][jcomp][-1]) )

# NBI geometry for plots
    print('RUNID = %s' %runid)

#------
# Plots
#------

# Above view

    phi_dep = np.radians(tor_ang)
    cols = ('r', 'b', 'g', 'm', 'y')

# Component by component: v_par/v_perp

    nrows = 2
    ncols = n_comp + 1

    delta_pit = 0.2

    jfig = 12

    xgrid, ygrid = np.meshgrid(R_grid, z_grid)

#------
# Plots
#------

    fsize = 12

# One tab for each source

    if topframe is None:
        topframe = tk.Toplevel()
        topframe.title('Birth location')
        topframe.geometry('1500x940')
    nbsource = ttk.Notebook(topframe, name='nb source')
    nbsource.pack(side=tk.TOP, fill=tk.X)

    for j, jsrc in enumerate(src_arr):
        frame_source = tk.Frame(nbsource)
        lbl = ' NBI #%d ' %jsrc

        nbbirth = ttk.Notebook(frame_source, name='nbbirth')
        nbbirth.pack(side=tk.TOP, fill=tk.X)

        frame_birth = tk.Frame(nbbirth)

        fig_birth = Figure(figsize=(14., 8.45), dpi=100)
        can_birth = FigureCanvasTkAgg(fig_birth, master=frame_birth)
        can_birth._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_birth.clf()
        fig_birth.subplots_adjust(left=0.05, bottom=0.08, right=0.97, top=0.97,  \
                                  wspace=0.15, hspace=0.)

        fig_birth.text(0.5, 0.95, '%s, t =%6.3f s, top view' %(runid, t_birth), ha='center')
        fig_birth.text(0.5, 0.55, 'Poloidal section'  , ha='center')

        jsplot = 1

# Overplot 3 species mix
# Above view
        axtop = fig_birth.add_subplot(nrows, ncols, jsplot, aspect='equal')
        axtop.set_title('All energy components', fontsize=fsize)

# Poloidal section

        axpol = fig_birth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')
        axpol.set_title('All energy components', fontsize=fsize)

        jcol=0
        for jcomp in comp:
            ind = np.where((mix_lbl == jcomp) & (j_nbi == jsrc))
            axtop.plot(Rj[ind]*np.cos(phi_dep[ind]), Rj[ind]*np.sin(phi_dep[ind]), cols[jcol]+'o', label=jcomp)
            axpol.plot(Rj[ind], zj[ind], '%so' %cols[jcol], label=jcomp)
            jcol += 1
        axtop.legend(loc=2, numpoints=1, prop={'size': 8})
        axpol.legend(loc=2, numpoints=1, prop={'size': 8})

        if fbm is not None:
            plot_sections.tor_sect(axtop, fbm.rsurf)
        axtop.plot(xlin[j], ylin[j], 'g-', linewidth=2.5)
        plot_sections.pol_sect(axpol)
        plot_sections.pol_cells(axpol, fbm)
        try:
            axpol.plot(rlin[j], zlin[j], 'g-', linewidth=2.5)
        except:
            pass

# For each species overplot 5 pitch angle range

        jcol = 0
        jsplot = 2

        for jcomp in comp:
            axtop = fig_birth.add_subplot(nrows, ncols, jsplot      , aspect='equal')
            axpol = fig_birth.add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')

            axtop.set_title('%s energy' %jcomp, fontsize=fsize)
            axpol.set_title('%s energy' %jcomp, fontsize=fsize)

            p1 = 0
            jcol = 0
            while p1 <= 1-delta_pit:
                p2 = p1 + delta_pit
                ind = np.where((mix_lbl == jcomp) & (j_nbi == jsrc) & \
                               (pitch > p1) & (pitch < p2)) 
                axtop.plot(Rj[ind]*np.cos(phi_dep[ind]), \
                           Rj[ind]*np.sin(phi_dep[ind]), '%so' %cols[jcol], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))
                axpol.plot(Rj[ind], zj[ind], '%so' %cols[jcol], \
                           label='%3.1f < p.a. < %3.1f' %(p1, p2))
                p1 += delta_pit
                jcol += 1

            if fbm is not None:
                plot_sections.tor_sect(axtop, fbm.rsurf)
            try:
                axtop.plot(xlin[j], ylin[j], 'g-', linewidth=2.5)
            except:
                pass
            plot_sections.pol_sect(axpol)
            plot_sections.pol_cells(axpol, fbm)
            try:
                axpol.plot(rlin[j], zlin[j], 'g-', linewidth=2.5)
            except:
                pass

            axtop.legend(loc=2, numpoints=1, prop={'size': 8})
            axpol.legend(loc=2, numpoints=1, prop={'size': 8})
            jsplot += 1

        toolbar = NavigationToolbar2TkAgg(can_birth, frame_birth)
        toolbar.update()

#-------------------------        
# Deposition & attenuation
#-------------------------

        frame_dep = tk.Frame(nbbirth)

        fig_pol = Figure(figsize=(6., 3.9), dpi=100)
        fig_pol.clf()
        frame_pol = tk.Frame(frame_dep, height=600)
        frame_pol.pack(side=tk.TOP, fill=tk.BOTH)
        frame_pol.pack_propagate(0)
        can_pol = FigureCanvasTkAgg(fig_pol, master=frame_pol)
        can_pol._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_pol.subplots_adjust(left=0.05, bottom=0.1, right=0.98, \
                                top=0.92)

        fig_pol.text(0.33, 0.95, '%s, t =%6.3f s' %(runid, t_birth), ha='center')
# 2D deposition, poloidal section
        jsplot = 1
        for jcomp in comp:
            zgrid = dep_matrix[jsrc][jcomp].T
            ind = np.where(zgrid == 0)
            zgrid[ind] = None
            axpol = fig_pol.add_subplot(1, n_comp, jsplot, aspect='equal')

            axpol.set_title('%s energy' %jcomp, fontsize=fsize)
            ctr = axpol.contourf(xgrid, ygrid, zgrid)
            fig_pol.colorbar(ctr, shrink=0.9, aspect=10)

            plot_sections.pol_sect(axpol)
            plot_sections.pol_cells(axpol, fbm)
            try:
                axpol.plot(rlin[j], zlin[j], 'g-', linewidth=2.5)
            except:
                pass

            jsplot += 1

        toolbar = NavigationToolbar2TkAgg(can_pol, frame_pol)
        toolbar.update()

# Attenutation

        fig_att = Figure(figsize=(3., 2.4), dpi=100)
        fig_att.clf()
        frame_att = tk.Frame(frame_dep)
        frame_att.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        frame_att.pack_propagate(0)
        can_att = FigureCanvasTkAgg(fig_att, master=frame_att)
        can_att._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig_att.subplots_adjust(left=0.05, bottom=0.2, right=0.98, \
                                top=0.9)

        jsplot = 1
        for jcomp in comp:
            axatt = fig_att.add_subplot(1, n_comp, jsplot)
            axatt.set_title('%s energy' %jcomp, fontsize=fsize)
            axatt.set_xlabel('R [cm]', fontsize=fsize)
            axatt.set_ylabel('NBI attenuation', fontsize=fsize)
            axatt.plot(R_grid, res_R[jsrc][jcomp])
            jsplot += 1

        toolbar = NavigationToolbar2TkAgg(can_att, frame_att)
        toolbar.update()

# Add tabs

        nbbirth.add(frame_birth, text='Birth location')
        nbbirth.add(frame_dep  , text='NBI penetration')

        nbsource.add(frame_source, text=lbl)


if __name__ == "__main__":

    fbirth = '/afs/ipp/home/g/git/tr_client/AUGD/29795/A05/29795A05_birth.cdf1'
    birth_tk = tk.Tk()
    birth_tk.title('Birth location')
    birth_tk.geometry('1500x940')
    read_birth(fbirth, topframe=birth_tk)
    birth_tk.mainloop()
