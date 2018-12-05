import sys, os
from scipy.io import netcdf
try:
    import Tkinter as tk
    import ttk
    import tkFileDialog as tkfd
    import tkMessageBox as tkmb
except:
    import tkinter as tk
    from tkinter import ttk
    from tkinter import filedialog as tkfd
    from tkinter import messagebox as tkmb

import read_fbm, plot_birth, plot_sections
import tkhyper
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.mlab import griddata
import matplotlib as mpl
import numpy as np


bdens_d = {'D_NBI': 'BDENS', 'He3_FUSN': 'FDENS_3', \
           'H_FUSN': 'FDENS_P', 'T_FUSN': 'FDENS_T'}

lframe_wid = 630
rframe_wid = 750
fsize = 12


class FBM:


    def __init__(self):

# Widget frame

        viewer = tk.Tk()
        viewer.title('FBM viewer')
        viewer.geometry('%dx960' %(lframe_wid + rframe_wid))
        viewer.option_add("*Dialog.msg.wrapLength", "10i")
        viewer.option_add("*Font", "Helvetica")
# Menubar

        menubar = tk.Menu(viewer)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Read FBM..." , command=self.load_fbm)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.about)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Help", menu=helpmenu)

        viewer.config(menu = menubar)

        nb = ttk.Notebook(viewer, name='nb')
        nb.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.fbmframe = ttk.Frame(nb)
        self.birthframe = ttk.Frame(nb)
        self.neutframe  = ttk.Frame(nb)
        self.trapframe  = ttk.Frame(nb)

        nb.add(self.fbmframe  , text='2D dist')
        nb.add(self.trapframe , text='Trap. frac.')
        nb.add(self.neutframe , text='Bt BB neut')
        nb.add(self.birthframe, text='Birth profile')

#----------
# FBM plots
#----------

        polframe = ttk.Frame(self.fbmframe, width=rframe_wid)
        polframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        polframe.pack_propagate(0)

        r_fbmframe = ttk.Frame(self.fbmframe, width=rframe_wid)
        r_fbmframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        r_fbmframe.pack_propagate(0)

        comm_fbmframe = ttk.Frame(r_fbmframe)
        comm_fbmframe.pack(side=tk.TOP, fill=tk.BOTH)
        bdens_frame = ttk.Frame(r_fbmframe)
        bdens_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        cell_fbmframe = ttk.Frame(r_fbmframe)
        cell_fbmframe.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.butt_d = {}
        self.buttons(comm_fbmframe, 'D_NBI')

# Poloidal canvas 

        fig_pol = Figure()
        can_pol = FigureCanvasTkAgg(fig_pol, master=polframe)
        can_pol._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        axpol = fig_pol.add_subplot(1, 1, 1, aspect='equal')

        try:
            plot_sections.pol_sect(axpol)
        except:
            pass

        self.can_fbm   = {}
        self.cell_mark = {}
        self.fig_cell = {}

# Bdens plot

        self.plot_bdens(bdens_frame, 'D_NBI')

# Phase-space plot

        fig_cell  = Figure(figsize=(3.5, 3.), dpi=100)
        axcell = fig_cell.add_subplot(1, 1, 1)
        fig_cell.subplots_adjust(left=0.14, bottom=0.15, right=0.82, top=0.92)
        can_cell = FigureCanvasTkAgg(fig_cell, master=cell_fbmframe)
        can_cell._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        axcell.set_xlabel('Energy [keV]', fontsize=fsize)
        axcell.set_ylabel('Pitch angle' , fontsize=fsize)
        axcell.set_ylim([-1,1])

        viewer.mainloop()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/aug/manuals/transp/fbm/plot_fbm.html">FBM homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "300x60")


    def buttons(self, frame, spc_lbl):
# Buttons

        box = {}
        n_hframes = 6
        for jframe in range(n_hframes):
            box[jframe] = ttk.Frame(frame)
            box[jframe].pack(side=tk.TOP, fill=tk.X, pady=2)

        self.butt_d[spc_lbl] = {}
        self.butt_d[spc_lbl]['th_int']    = tk.BooleanVar()
        self.butt_d[spc_lbl]['vol_int']   = tk.BooleanVar()
        self.butt_d[spc_lbl]['log_scale'] = tk.BooleanVar()
        self.butt_d[spc_lbl]['log_scale'].set(True)

        lb1 = ttk.Label(box[0], text='Right-mouse click on a cell for local FBM plot')

        cb1 = ttk.Checkbutton(box[2], text='Theta averaged FBM', variable=self.butt_d[spc_lbl]['th_int'])
        cb2 = ttk.Checkbutton(box[3], text='Volume averaged FBM', variable=self.butt_d[spc_lbl]['vol_int'])

        for but in (cb1, cb2, lb1):
            but.pack(side=tk.TOP, anchor=tk.W, padx=10, pady=2)

        Elbl = ttk.Label(box[4], text='Emax [keV]')
        self.butt_d[spc_lbl]['Emax'] = ttk.Entry(box[4], width=12)
        self.butt_d[spc_lbl]['Emax'].insert(0, 100)
        for but in (Elbl, self.butt_d[spc_lbl]['Emax']):
            but.pack(side=tk.LEFT, anchor=tk.W, padx=10, pady=2)

        cb3 = ttk.Checkbutton(box[5], text='Log scale', variable=self.butt_d[spc_lbl]['log_scale'])
        fminlbl = ttk.Label(box[5], text='f_log_min')
        fmaxlbl = ttk.Label(box[5], text='f_log_max')
        self.butt_d[spc_lbl]['fmin'] = ttk.Entry(box[5], width=12)
        self.butt_d[spc_lbl]['fmin'].insert(0, 4.5)
        self.butt_d[spc_lbl]['fmax'] = ttk.Entry(box[5], width=12)
        self.butt_d[spc_lbl]['fmax'].insert(0, 8.5)
        for but in (cb3, fminlbl, self.butt_d[spc_lbl]['fmin'], fmaxlbl, self.butt_d[spc_lbl]['fmax']):
            but.pack(side=tk.LEFT, anchor=tk.W, padx=10, pady=5)



    def plot_bdens(self, frame, spc_lbl):

        for child in frame.winfo_children():
            child.destroy()
        
        fig_bdens = Figure(figsize=(3.5, 3.), dpi=100)
        axbdens = fig_bdens.add_subplot(1, 1, 1)
        fig_bdens.subplots_adjust(left=0.13, bottom=0.15, right=0.95, top=0.94)
        can_bdens = FigureCanvasTkAgg(fig_bdens, master=frame)
        can_bdens._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        axbdens.set_xlabel(r'$\rho_{tor}$',fontsize=fsize)
        axbdens.set_ylabel('%s [1/cm**3]' %bdens_d[spc_lbl],fontsize=fsize)
        axbdens.set_xlim([0, 1])

        if not hasattr(self, 'fbmr'):
            axbdens.plot([], [], 'r-', label='From FBM')
            axbdens.plot([], [], 'g-', label='From CDF')
        else:
            axbdens.plot(self.fbmr.rho_grid, self.fbmr.bdens[spc_lbl])
            axbdens.plot(self.cv['X'], self.cv[bdens_d[spc_lbl]])

        axbdens.legend()
        can_bdens.draw()


    def plot_trapdens(self, frame, spc_lbl):

        for child in frame.winfo_children():
            child.destroy()

        fig_trapdens = Figure(figsize=(5.5, 4.), dpi=100)
        axtrapdens = fig_trapdens.add_subplot(1, 1, 1)
        fig_trapdens.subplots_adjust(left=0.13, bottom=0.15, right=0.95, top=0.94)
        can_trapdens = FigureCanvasTkAgg(fig_trapdens, master=frame)
        can_trapdens._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH)
        axtrapdens.set_xlabel(r'$\rho_{tor}$',fontsize=fsize)
        axtrapdens.set_ylabel('Trapped fraction',fontsize=fsize)
        axtrapdens.set_xlim([0, 1])
        axtrapdens.set_ylim([0, 1])

        if hasattr(self, 'fbmr'):
            frac_trap = self.fbmr.btrap[spc_lbl]/self.fbmr.bdens[spc_lbl]
            axtrapdens.plot(self.fbmr.rho_grid, frac_trap, 'r-')
        else:
            axtrapdens.plot([], [])
        toolbar = NavigationToolbar2TkAgg(can_trapdens, frame)
        toolbar.update()

        can_trapdens.draw()


    def fig_plot(self, frame, title, zdata, zmin=None, zmax=None, spc_lbl=None):

        fbmfile  = self.ffbm.split('/')[-1]
        runid    = fbmfile[:8]

        print('Plot %s' %title)
        title += ', Run %s, t =%6.3f s' %(runid, self.fbmr.time)

        for child in frame.winfo_children():
            child.destroy()

        fig = Figure()
        can = FigureCanvasTkAgg(fig, master=frame)
        can._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        fig.subplots_adjust(left=0.08, bottom=0.08, right=0.8, top=0.92)
        fig.text(0.5, 0.95, title, ha='center')

        if not hasattr(self, 'x_grid'):
            xgrid = np.linspace(self.fbmr.r2d.min(), self.fbmr.r2d.max(), 1000)
            ygrid = np.linspace(self.fbmr.z2d.min(), self.fbmr.z2d.max(), 1000)
            self.x_grid, self.y_grid = np.meshgrid(xgrid, ygrid)

        ax = fig.add_subplot(1, 1, 1, aspect='equal')

        plot_sections.pol_sect(ax)
        plot_sections.pol_cells(ax, self.fbmr)

# Selected point

        ctr_f = griddata(self.fbmr.r2d, self.fbmr.z2d, zdata, \
             self.x_grid, self.y_grid, interp='linear')
        if zmin is None:
            zmin = np.nanmin(ctr_f)
        if zmax is None:
            zmax = np.nanmax(ctr_f)
        n_levels = 11
        bounds = np.linspace(zmin, zmax, n_levels)
        ax.contourf(self.x_grid, self.y_grid, ctr_f, bounds)
        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
        cb_ax = fig.add_axes([0.82, 0.08, 0.09, 0.84])
        mpl.colorbar.ColorbarBase(cb_ax,norm=norm,boundaries=bounds,ticks=bounds)

        if spc_lbl is not None:
            self.can_fbm[spc_lbl]    = can
            self.cell_mark[spc_lbl], = ax.plot([], [], 'ro')

# Toolbar
        toolbar = NavigationToolbar2TkAgg(can, frame)
        toolbar.update()


    def load_fbm(self):

        dir_init = os.getenv('HOME')+'/tr_client/AUGD'
        self.ffbm = tkfd.askopenfilename(  initialdir=dir_init, filetypes= \
                                        (("FBM", "*.cdf"),)    )

        print(self.ffbm)

        if self.ffbm.strip() == '':
            tkmb.showerror("Error", 'Select a FBM file: run File->read_fbm')
        else:
            self.read_all()


    def read_all(self):

        fbmdir, fbmfile  = os.path.split(self.ffbm)

        tmp = fbmfile.split('_')
        runid = tmp[0]
        t_id  = tmp[2].split('.')[0]

# Read FBM 
        self.fbmr = read_fbm.READ_FBM(self.ffbm)

        self.species = []
        for spc_lbl in self.fbmr.int_en_pit_frac_trap.keys():
            self.species.append(spc_lbl)

# Read CDF
        tr_file = '%s/%s.CDF' %(fbmdir, runid)
        if not os.path.isfile(tr_file):
            print('%s not found' %tr_file)
            tkmb.showerror("Error", '%s not found')
        cv_all = netcdf.netcdf_file(tr_file, 'r', mmap=False).variables
        sigs = ['BTNT2_DD', 'BBNT2_DD', 'TIME3', 'X'] + bdens_d.values()
        tdist = (cv_all['TIME3'].data - self.fbmr.time)**2
        jtclose = np.argmin(tdist)
        self.cv = {}
        for lbl in sigs:
            if lbl in cv_all.keys():
                self.cv[lbl] = cv_all[lbl][jtclose]

# Plots

        if hasattr(self, 'x_grid'):
            del self.x_grid, self.y_grid

        self.plot_neut(self.neutframe)
        self.plot_dist(self.fbmframe)
        self.plot_trap(self.trapframe)
        birth_file =  '%s/%s_birth.cdf%s' %(fbmdir, runid, t_id)
        plot_birth.read_birth(birth_file, topframe=self.birthframe)


    def plot_neut(self, frame):

        for child in frame.winfo_children():
            child.destroy()

        btneut = self.cv['BTNT2_DD']
        bbneut = self.cv['BBNT2_DD']
        indbt = (btneut == 0)
        indbb = (btneut == 0)
        btneut[indbt] = np.nan
        bbneut[indbb] = np.nan
        print('Time in CDF and cdf: %6.3f, %6.3f s' %(self.cv['TIME3'], self.fbmr.time) )

        btframe = ttk.Frame(frame)
        bbframe = ttk.Frame(frame)
        for frame in btframe, bbframe:
            frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        self.fig_plot(btframe, 'Beam-target neutrons', btneut)
        self.fig_plot(bbframe, 'Beam-beam neutrons'  , bbneut)


    def plot_trap(self, frame):

        for child in frame.winfo_children():
            child.destroy()

        nbtrap = ttk.Notebook(frame, name='nb')
        nbtrap.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        for spc_lbl in self.species:
            trap_frame = ttk.Frame(nbtrap)
            pol_trapframe = ttk.Frame(trap_frame, width=rframe_wid)
            prof_trapframe = ttk.Frame(trap_frame, width=rframe_wid)
            pol_trapframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            pol_trapframe.pack_propagate(0)
            prof_trapframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            prof_trapframe.pack_propagate(0)
            nbtrap.add(trap_frame, text=spc_lbl)

            title = 'Trapped fast ion fraction %s' %spc_lbl
            self.fig_plot(pol_trapframe, title, self.fbmr.int_en_pit_frac_trap[spc_lbl], zmin=0)
            self.plot_trapdens(prof_trapframe, spc_lbl)


    def plot_dist(self, frame):

        for child in frame.winfo_children():
            child.destroy()

        self.nbdist = ttk.Notebook(frame, name='nb')
        self.nbdist.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        for spc_lbl in self.species:
            nbframe = ttk.Frame(self.nbdist)
            pol_frame = ttk.Frame(nbframe, width=rframe_wid)
            pol_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            pol_frame.pack_propagate(0)
            r_fbmframe = ttk.Frame(nbframe, width=rframe_wid)
            r_fbmframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
            r_fbmframe.pack_propagate(0)

            comm_fbmframe = ttk.Frame(r_fbmframe)
            comm_fbmframe.pack(side=tk.TOP, fill=tk.BOTH)
            bdens_frame = ttk.Frame(r_fbmframe, height=200)
            bdens_frame.pack(side=tk.TOP, fill=tk.BOTH)
            cell_frame = ttk.Frame(r_fbmframe)
            cell_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

            self.nbdist.add(nbframe, text=spc_lbl)

            title = r'2D distribution %s, $\int\int$ dE dp.a.' %spc_lbl
            self.fig_plot(pol_frame, title, self.fbmr.int_en_pit[spc_lbl], zmin=0, spc_lbl=spc_lbl)

            self.fig_cell[spc_lbl]  = Figure(figsize=(3.5, 3.), dpi=100)
            self.fig_cell[spc_lbl].subplots_adjust(left=0.14, bottom=0.15, right=0.82, top=0.92)
            can_cell = FigureCanvasTkAgg(self.fig_cell[spc_lbl], master=cell_frame)
            can_cell._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(can_cell, cell_frame)
            toolbar.update()

# Plots on the right
            self.can_fbm[spc_lbl].mpl_connect('button_press_event', self.my_call)

            self.buttons(comm_fbmframe, spc_lbl)
            self.plot_bdens(bdens_frame, spc_lbl)
            self.plot_fbm_cell(0, spc_lbl)


    def my_call(self, event):

        if event.button in (2,3):
            dist2 = (self.fbmr.r2d - event.xdata)**2 + (self.fbmr.z2d - event.ydata)**2
            jcell = np.argmin(dist2)
            self.plot_fbm_cell(jcell)


    def plot_fbm_cell(self, jcell, spc_lbl=None):

        if spc_lbl is None:
            spc_lbl = self.species[self.nbdist.index('current')]

        self.fig_cell[spc_lbl].clf()
        ax = self.fig_cell[spc_lbl].gca()
        ax.set_ylim((-1, 1))
        ax.set_xlabel('Energy [keV]', fontsize=fsize)
        ax.set_ylabel('Pitch angle' , fontsize=fsize)

        thint    = self.butt_d[spc_lbl]['th_int'].get()
        volint   = self.butt_d[spc_lbl]['vol_int'].get()
        logscale = self.butt_d[spc_lbl]['log_scale'].get()

        jrho = np.where(self.fbmr.rho_grid == self.fbmr.x2d[jcell])[0][0]

        if volint:
            self.cell_mark[spc_lbl].set_data(self.fbmr.r2d, self.fbmr.z2d)
            tit_lbl = 'Volume averaged, t=%6.3f' %self.fbmr.time
            zarr_lin = self.fbmr.dens_vol[spc_lbl]
        else:
            if thint:
                ind=np.where(self.fbmr.x2d == self.fbmr.x2d[jcell])
                self.cell_mark[spc_lbl].set_data(self.fbmr.r2d[ind], self.fbmr.z2d[ind])
                tit_lbl = r'$\rho_{tor}$' + \
                          r' = %8.4f, $\theta$ averaged, t=%6.3f' %(self.fbmr.x2d[jcell], self.fbmr.time)
                zarr_lin = self.fbmr.dens_zone[spc_lbl][jrho]
            else:
                self.cell_mark[spc_lbl].set_data(self.fbmr.r2d[jcell], self.fbmr.z2d[jcell])
                tit_lbl = r'$\rho_{tor} = $ %8.4f $\theta = $%8.4f deg, t=%6.3f s' \
                          %(self.fbmr.x2d[jcell], np.degrees(self.fbmr.th2d[jcell]), self.fbmr.time)
                zarr_lin = self.fbmr.fdist[spc_lbl][jcell]
        self.can_fbm[spc_lbl].draw() # Update red marker

        zmax_lin = np.max(zarr_lin[~np.isnan(zarr_lin)])
        zmin_lin = 1e-8*zmax_lin
#        flog_max = np.log10(zmax_lin)
#        flog_min = np.log10(zmin_lin)
        flog_min = float(self.butt_d[spc_lbl]['fmin'].get())
        flog_max = float(self.butt_d[spc_lbl]['fmax'].get())
        indzero = np.where(zarr_lin <= zmin_lin)
        zarr_lin[indzero] = np.nan
        zarr_log = np.log10(zarr_lin)
        n_levels = 15
        if logscale:
            zmin = flog_min
            zmax = flog_max
            zarr = zarr_log
        else:
            zmin = 0
            zmax = zmax_lin
            zarr = zarr_lin

        bounds = np.linspace(zmin, zmax, n_levels)

        Emax = float(self.butt_d[spc_lbl]['Emax'].get())
        ax.set_xlim((0, Emax))
        ax.set_title(tit_lbl, fontsize=fsize)

        ctr = ax.contourf(1e-3*self.fbmr.e_d[spc_lbl], self.fbmr.a_d[spc_lbl], zarr, bounds)
        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
        cb_ax = self.fig_cell[spc_lbl].add_axes([0.86, 0.05, 0.05, 0.91])
        mpl.colorbar.ColorbarBase(cb_ax,norm=norm,boundaries=bounds,ticks=bounds)
        ax.plot([0, Emax], [0, 0], 'k-')
        if not volint and not thint:
            ax.plot([0, Emax], [ self.fbmr.trap_pit[spc_lbl][jcell],  self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')
            ax.plot([0, Emax], [-self.fbmr.trap_pit[spc_lbl][jcell], -self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')

        ax.figure.canvas.draw()


if __name__ == "__main__":

    FBM()
