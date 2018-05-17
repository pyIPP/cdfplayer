import sys, os
from scipy.io import netcdf
import Tkinter as tk
import ttk
from tkFileDialog import *
from tkMessageBox import *
import read_fbm, plot_birth, plot_sections
import tkhyper
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.mlab import griddata
import matplotlib as mpl
import numpy as np


class FBM:


    fsize = 12


    def __init__(self):

# Widget frame

        polframe = tk.Tk()
        polframe.title('FBM viewer')
        polframe.geometry('1500x960')
        polframe.option_add("*Dialog.msg.wrapLength", "10i")

# Menubar

        menubar = tk.Menu(polframe)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Read FBM..." , command=self.load_fbm)
        filemenu.add_command(label="Plot dist...", command=self.plot_dist)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.about)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Help", menu=helpmenu)

        polframe.config(menu = menubar)

        nb = ttk.Notebook(polframe, name='nb')
        nb.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        fbmplotsframe   = tk.Frame(nb)
        self.birthframe = tk.Frame(nb)
        self.neutframe  = tk.Frame(nb)
        self.trapframe  = tk.Frame(nb)

        nb.add(fbmplotsframe  , text='2D dist')
        nb.add(self.trapframe , text='Trap. frac.')
        nb.add(self.neutframe , text='Bt BB neut')
        nb.add(self.birthframe, text='Birth profile')

        self.distframe = tk.Frame(fbmplotsframe, width=630)
        self.distframe.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        self.distframe.pack_propagate(0)

        right_frame = tk.Frame(fbmplotsframe, width=870)
        right_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        right_frame.pack_propagate(0)

        comm_frame = tk.Frame(right_frame)
        comm_frame.pack(side=tk.TOP, fill=tk.BOTH)
        bdens_frame = tk.Frame(right_frame)
        bdens_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        cell_frame = tk.Frame(right_frame)
        cell_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Buttons

        box = {}
        n_hframes = 6
        for jframe in range(n_hframes):
            box[jframe] = tk.Frame(comm_frame)
            box[jframe].pack(side=tk.TOP, fill=tk.X, pady=2)

        self.th_int    = tk.BooleanVar()
        self.vol_int   = tk.BooleanVar()
        self.log_scale = tk.BooleanVar()
        self.log_scale.set(True)

        lb1 = tk.Label(box[0], text='Right mouse click on a cell for FBM', fg="#900000")

        cb1 = tk.Checkbutton(box[2], text='Theta averaged FBM', variable=self.th_int)
        cb2 = tk.Checkbutton(box[3], text='Volume averaged FBM', variable=self.vol_int)

        for but in (cb1, cb2, lb1):
            but.pack(side=tk.TOP, anchor=tk.W, padx=10, pady=2)

        Elbl = tk.Label(box[4], text='Emax [keV]')
        self.Emax = tk.Entry(box[4], width=12, bg='#ffffff')
        self.Emax.insert(0, 100)
        for but in (Elbl, self.Emax):
            but.pack(side=tk.LEFT, anchor=tk.W, padx=10, pady=2)

        cb3 = tk.Checkbutton(box[5], text='Log scale', variable=self.log_scale)
        fminlbl = tk.Label(box[5], text='f_log_min', pady=5)
        fmaxlbl = tk.Label(box[5], text='f_log_max', pady=5)
        self.fmin = tk.Entry(box[5], width=12, bg='#ffffff')
        self.fmin.insert(0, 4.5)
        self.fmax = tk.Entry(box[5], width=12, bg='#ffffff')
        self.fmax.insert(0, 8.5)
        for but in (cb3, fminlbl, self.fmin, fmaxlbl, self.fmax):
            but.pack(side=tk.LEFT, anchor=tk.W, padx=10, pady=2)

# Bdens plot

        fig_bdens = Figure(figsize=(3.5, 3.), dpi=100)
        axbdens = fig_bdens.add_subplot(1, 1, 1)
        fig_bdens.subplots_adjust(left=0.13, bottom=0.15, right=0.95, top=0.94)
        self.can_bdens = FigureCanvasTkAgg(fig_bdens, master=bdens_frame)
        self.can_bdens._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        axbdens.set_xlabel(r'$\rho_{tor}$',fontsize=self.fsize)
        axbdens.set_ylabel('BDENS [1/cm**3]',fontsize=self.fsize)
        axbdens.set_xlim([0, 1])
        axbdens.set_ylim([0, 1e13])

        self.bdplot1, = axbdens.plot([], [], 'r-', label='From FBM')
        self.bdplot2, = axbdens.plot([], [], 'g-', label='From CDF')
        axbdens.legend()

# Cell plot

        self.fig_cell  = Figure(figsize=(3.5, 3.), dpi=100)
        self.axcell = self.fig_cell.add_subplot(1, 1, 1)
        self.fig_cell.subplots_adjust(left=0.14, bottom=0.15, right=0.82, top=0.92)
        self.can_cell = FigureCanvasTkAgg(self.fig_cell, master=cell_frame)
        self.can_cell._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.axcell.set_xlabel('Energy [keV]', fontsize=self.fsize)
        self.axcell.set_ylabel('Pitch angle' , fontsize=self.fsize)
        self.axcell.set_ylim([-1,1])
        toolbar = NavigationToolbar2TkAgg(self.can_cell, cell_frame)
        toolbar.update()


        self.can_fbm   = {}
        self.cell_mark = {}

# Poloidal canvas 

        fig_pol = Figure()
        self.can_pol = FigureCanvasTkAgg(fig_pol, master=self.distframe)
        self.can_pol._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        axpol = fig_pol.add_subplot(1, 1, 1, aspect='equal')

        try:
            plot_sections.pol_sect(axpol)
        except:
            pass

        polframe.mainloop()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/aug/manuals/transp/fbm/plot_fbm.html">FBM homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "300x60")


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
        self.ffbm = askopenfilename(  initialdir=dir_init, filetypes= \
                                    (("FBM", "*.cdf"),)    )

        print(self.ffbm)

        if self.ffbm.strip() == '':
            showerror("Error", 'Select a FBM file: run File->read_fbm')
        else:
            self.read_fbm()


    def read_fbm(self):

        fbmfile  = self.ffbm.split('/')[-1]
        runid    = fbmfile[:8]
        t_id     = fbmfile[12:13]
        self.fbmr = read_fbm.READ_FBM(self.ffbm)
        self.plot_dist()
        birth_file =  '%s_birth.cdf%s' %(self.ffbm[:-9], t_id)
        plot_birth.read_birth(birth_file, topframe=self.birthframe)


    def plot_dist(self):
# Neutrons

        tr_file = self.ffbm[:-9] + '.CDF'
        if not os.path.isfile(tr_file):
            print('%s not found' %tr_file)
            showerror("Error", '%s not found')
            
        cv = netcdf.netcdf_file(tr_file, 'r', mmap=False).variables

# Select time
        tdist = (cv['TIME3'].data - self.fbmr.time)**2
        jtclose = np.argmin(tdist)
        btneut = cv['BTNT2_DD'][jtclose, :]
        bbneut = cv['BBNT2_DD'][jtclose, :]
        indbt = (btneut == 0)
        indbb = (btneut == 0)
        btneut[indbt] = np.nan
        bbneut[indbb] = np.nan
        print('Time in CDF and cdf: %6.3f, %6.3f s' %(cv['TIME3'][jtclose], self.fbmr.time) )

#------
# Plots
#------

# Right frame: BDENS

# BDENS canvas

        self.bdplot1.set_data(self.fbmr.rho_grid, self.fbmr.bdens['D_NBI'])
        self.bdplot2.set_data(cv['X'][jtclose, :], cv['BDENS'][jtclose, :])
        self.can_bdens.draw()

# Left frame: tabs

        if hasattr(self, 'can_pol'):
            self.can_pol.get_tk_widget().destroy()

        if hasattr(self, 'x_grid'):
            del self.x_grid, self.y_grid
        
# Neutrons
        for child in self.neutframe.winfo_children():
            child.destroy()
        btframe = tk.Frame(self.neutframe)
        bbframe = tk.Frame(self.neutframe)
        for frame in btframe, bbframe:
            frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)

        self.fig_plot(btframe, 'Beam-target neutrons', btneut)
        self.fig_plot(bbframe, 'Beam-beam neutrons'  , bbneut)

# FBM

        for child in self.distframe.winfo_children():
            child.destroy()
        for child in self.trapframe.winfo_children():
            child.destroy()

        self.nbdist = ttk.Notebook(self.distframe, name='nb')
        self.nbdist.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        nbtrap = ttk.Notebook(self.trapframe, name='nb')
        nbtrap.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        self.species = []
        for spc_lbl in self.fbmr.frac_trap.keys():
            self.species.append(spc_lbl)
            dframe = tk.Frame(self.nbdist)
            self.nbdist.add(dframe, text=spc_lbl)
            trapframe = tk.Frame(nbtrap)
            nbtrap.add(trapframe, text=spc_lbl)

            title = r'2D distribution %s, $\int\int$ dE dp.a.' %spc_lbl
            self.fig_plot(dframe, title, self.fbmr.int_en_pit[spc_lbl], zmin=0, spc_lbl=spc_lbl)
            self.can_fbm[spc_lbl].mpl_connect('button_press_event', self.my_call)

# Trapped fraction
            title = 'Trapped fast ion fraction %s' %spc_lbl
            self.fig_plot(trapframe, title, self.fbmr.frac_trap[spc_lbl], zmin=0)

# Plot 1st cell by default
            self.plot_cell(0)


    def my_call(self, event):

        if event.button in (2,3):
            dist2 = (self.fbmr.r2d - event.xdata)**2 + (self.fbmr.z2d - event.ydata)**2
            jcell = np.argmin(dist2)
            self.plot_cell(jcell)


    def plot_cell(self, jcell):

        thint    = self.th_int.get()
        volint   = self.vol_int.get()
        logscale = self.log_scale.get()

# "None" better then "0" for clarity in contour plot

        jrho = np.where(self.fbmr.rho_grid == self.fbmr.x2d[jcell])[0][0]
        spc_lbl = self.species[self.nbdist.index('current')]
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
        self.can_fbm[spc_lbl].draw()

        zmax_lin = np.max(zarr_lin[~np.isnan(zarr_lin)])
        zmin_lin = 1e-6*zmax_lin
#        flog_max = np.log10(zmax_lin)
#        flog_min = np.log10(zmin_lin)
        flog_min = float(self.fmin.get())
        flog_max = float(self.fmax.get())
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

        Emax = float(self.Emax.get())

        self.axcell.cla()
        self.axcell.set_xlim((0, Emax))
        self.axcell.set_ylim((-1, 1))
        self.axcell.set_xlabel('Energy [keV]', fontsize=self.fsize)
        self.axcell.set_ylabel('Pitch angle' , fontsize=self.fsize)
        self.axcell.set_title(tit_lbl, fontsize=self.fsize)

        ctr = self.axcell.contourf(1e-3*self.fbmr.e_d[spc_lbl], self.fbmr.a_d[spc_lbl], zarr, bounds)
        norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
        cb_ax = self.fig_cell.add_axes([0.86, 0.05, 0.05, 0.91])
        mpl.colorbar.ColorbarBase(cb_ax,norm=norm,boundaries=bounds,ticks=bounds)
        self.axcell.plot([0, Emax], [0, 0], 'k-')
        if not volint and not thint:
            self.axcell.plot([0, Emax], [ self.fbmr.trap_pit[spc_lbl][jcell],  self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')
            self.axcell.plot([0, Emax], [-self.fbmr.trap_pit[spc_lbl][jcell], -self.fbmr.trap_pit[spc_lbl][jcell]], 'g-')

        self.can_cell.draw()


if __name__ == "__main__":

    FBM()
