import sys,os
from netCDF4 import Dataset
import Tkinter as tk
from tkFileDialog import *
from tkMessageBox import *
import fconf, nbi_geo, los, tkhyper, fbm2ascii
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.mlab import griddata
import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np


class FBM:

    fsize = 12
    Rlbl = 'R [cm]'
    zlbl = 'z [cm]'
    MCgrid_col = 'red'

    err_msg = "Please browse a FBM file from the 'File' menu"


    def __init__(self, tok='AUGD'):

        if tok == 'AUGD':
            sys.path.append('/afs/ipp/aug/ads-diags/common/python/lib/')
            import plot_aug
            import kk_20140416 as kk
            self.kgc = kk.kkGCd0()
            self.tor_d = plot_aug.STRUCT().tor_old

# Widget frame

        if __name__ == '__main__':
            fbmframe = tk.Tk()
        else:
            fbmframe = tk.Toplevel()
        fbmframe.title('FBM viewer')
        fbmframe.geometry('1050x950')
        fbmframe.option_add("*Dialog.msg.wrapLength", "10i")

        canv_frame = tk.Frame(fbmframe, width=630)
        canv_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        canv_frame.pack_propagate(0)
        comm_frame = tk.Frame(fbmframe, width=320)
        comm_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=1)
        comm_frame.pack_propagate(0)

# Buttons

        box = {}
        n_hframes = 6
        for jframe in range(n_hframes):
            box[jframe] = tk.Frame(comm_frame)
            box[jframe].pack(side=tk.TOP, fill=tk.X)

        self.th_int    = tk.BooleanVar()
        self.rho_int   = tk.BooleanVar()
        self.log_scale = tk.BooleanVar()

        lb1 = tk.Label(box[0], text='Right mouse click on a cell for FBM', fg="#000000", pady=5)
        lb2 = tk.Label(box[1], text='No check for single cell FBM', fg="#900000", pady=5)
        cb1 = tk.Checkbutton(box[2], text='Theta averaged FBM', variable=self.th_int)
        cb2 = tk.Checkbutton(box[3], text='Rho+theta averaged FBM', variable=self.rho_int)
        cb3 = tk.Checkbutton(box[4], text='Log scale', variable=self.log_scale)

        self.fig_pol = Figure()
        self.canvas = FigureCanvasTkAgg(self.fig_pol, master=canv_frame)
        self.fig_pol.subplots_adjust(left=0.08, bottom=0.08, right=0.8, top=0.92)
        self.tit = self.fig_pol.text(0.5, 0.95,'Fast ion density, integrated over E and p.a.', \
                                     ha='center')
        self.ax3 = self.fig_pol.add_subplot(1, 1, 1, aspect='equal')
        self.ax3.set_xlabel(self.Rlbl, fontsize=self.fsize)
        self.ax3.set_ylabel(self.zlbl, fontsize=self.fsize)
        self.ax3.set_xlim((90, 230))
        self.ax3.set_ylim((-125, 125))
        try:
            for key in self.kgc.gc_x.iterkeys():
                self.ax3.plot(100*self.kgc.gc_x[key], 100*self.kgc.gc_y[key], 'b-')
        except:
            print('No plot of vessel components available')
        self.canvas.draw()

        for but in (cb1, cb2, cb3, lb1, lb2):
            but.pack(side=tk.TOP, anchor=tk.W, padx=10, pady=2)
        toolbar = NavigationToolbar2TkAgg(self.canvas, comm_frame)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Menubar

        menubar = tk.Menu(fbmframe)
        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Read FBM...", command=self.load_file)
        filemenu.add_command(label="Store FBM as ASCII", command=self.fbm2ascii)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)
        plotmenu = tk.Menu(menubar, tearoff=0)
        plotmenu.add_command(label="Plot FBM"               , command=self.plot_fbm)
        plotmenu.add_command(label="Plot trapped fraction"  , command=self.plot_trap)
        plotmenu.add_command(label="Plot GC birth locations", command=self.plot_birth)
        plotmenu.add_command(label="Plot fast ion density"  , command=self.plot_bdens)
        plotmenu.add_command(label="Plot b-th neutrons"     , command=self.plot_btneut)
        plotmenu.add_command(label="Plot b-b  neutrons"     , command=self.plot_bbneut)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.about)
        menubar.add_cascade(label="File", menu=filemenu)
        menubar.add_cascade(label="Plot", menu=plotmenu)
        menubar.add_cascade(label="Help", menu=helpmenu)

        fbmframe.config(menu = menubar)

        fbmframe.mainloop()


    def about(self):
        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/aug/manuals/transp/fbm/plot_fbm.html">FBM homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "300x60")


    def load_file(self):
        dir_init = os.getenv('HOME')+'/tr_client/AUGD'
        self.fname = askopenfilename(  initialdir=dir_init, filetypes= \
                                    (("FBM", "*.cdf"),)    )
        self.read_fbm()


    def fbm2ascii(self):

        try:
            fbm_file = self.fname
        except:
            showerror("Error", self.err_msg)
            return
            
        fbmfile = fbm_file.split('/')[-1]
        runid   = fbmfile[0:8]
        t_id    = fbmfile[12:13]

        if not os.path.isfile(fbm_file):
            print(fbm_file+' not found')
            sys.exit()
        fbm2ascii.FBM2ASCII(runid, t_id=t_id)


    def reset_plot(self):

        self.ax3.cla()
        self.ax3.set_xlabel(self.Rlbl, fontsize=self.fsize)
        self.ax3.set_ylabel(self.zlbl, fontsize=self.fsize)
        self.ax3.set_xlim((90, 230))
        self.ax3.set_ylim((-125, 125))

        try:
            for key in self.kgc.gc_x.iterkeys():
                self.ax3.plot(100*self.kgc.gc_x[key], 100*self.kgc.gc_y[key], 'b-')
        except:
            print('')
# Selected point
        self.line1, = self.ax3.plot([], [], 'ro')

# Cells

        for irho in self.ind_zone:
            self.ax3.plot(self.rsurf[irho, :], self.zsurf[irho, :], color=self.MCgrid_col)
        for jbar, myr in enumerate(self.rbar):
            self.ax3.plot(myr, self.zbar[jbar], color=self.MCgrid_col)
        self.canvas.draw()


    def read_fbm(self):

        fbm_file = self.fname
        fbmfile  = fbm_file.split('/')[-1]
        runid    = fbmfile[:8]
        t_id     = fbmfile[12:13]

        if not os.path.isfile(fbm_file):
            print('%s not found' %fbm_file)
            sys.exit()
        cdf = Dataset(fbm_file, 'r', format='NETCDF4')
        cv = cdf.variables

        print('PLOT_FBM')
        print(fbm_file)

        spc_lbl = "".join(cv['SPECIES_1'][:])

# Read data from CDF

        self.fbm  = cv['F_'+spc_lbl][:]
        self.a_d  = cv['A_'+spc_lbl][:]
        self.e_d  = cv['E_'+spc_lbl][:]*1.e-3
        self.r2d  = cv['R2D'][:]
        self.z2d  = cv['Z2D'][:]
        self.x2d  = cv['X2D'][:]
        self.th2d = cv['TH2D'][:]
        self.time = cv['TIME'][:]
        eb_d  = cv['EB_%s' %spc_lbl][:]
        bmvol = cv['BMVOL'][:]

        vol = 1.e-6*np.sum(bmvol)
        print('Volume is %8.4f m^-3' %vol)
        if cv['SYMFLAG'][0] == 0:
            n_sym = 2 # asymmetric equilibrium
        else:
            n_sym = 1

        n_cells, n_pit, n_E = self.fbm.shape
        n_zones = cv['N_ZONEROWS'][0]

        rho_lab = np.zeros(n_cells, dtype=int)    # rho index, takes values 0:n_zones-1
        the_lab = np.zeros(n_cells, dtype=int)    # theta index at different rho's

        self.rho_grid = np.zeros(n_zones) # rho_grid: rho array without redundances
        self.the_size = np.zeros(n_zones, dtype=int) # n of cells at given rho
        rmaj_min   = np.zeros(n_zones) # min(R(rho))
        the_grid2d = np.zeros(2*n_sym*n_zones) #TH2D reshaped
        self.thb_grid2d = np.zeros((n_zones, 2*n_sym*n_zones))

        max_jrho = 0
        for jrho in range(n_zones):
            self.the_size[jrho] = 2*n_sym*(jrho + 1) # amount of poloidal cells at jrho
            min_jrho = max_jrho
            max_jrho += self.the_size[jrho]
            rho_lab[min_jrho: max_jrho] = jrho
            self.rho_grid[jrho] = self.x2d[min_jrho]
            rmaj_min[jrho] = min(self.r2d[min_jrho: max_jrho])
            the_lab[min_jrho: max_jrho] = range(max_jrho - min_jrho)

            for i_the in range(self.the_size[jrho]):
                cell_lab = n_sym*jrho*(jrho + 1) + i_the # absolute cell index
                the_grid2d[i_the] = self.th2d[cell_lab]
            the_step = the_grid2d[1] - the_grid2d[0]
            self.thb_grid2d[jrho, :] = the_grid2d[:] - 0.5*the_step

# Boundary grids for theta and rho (both are equispaced in TRANSP)

        self.rho_b = np.zeros(n_zones + 1)
        rho_step = 1./n_zones
        self.rho_b[:n_zones] = self.rho_grid     - 0.5*rho_step
        self.rho_b[n_zones]  = self.rho_grid[-1] + 0.5*rho_step

        dE = np.diff(eb_d)
        dpit = 1./float(n_pit)

# Passing vs trapped

        fbm_pass = np.zeros((n_cells, n_E))
        fbm_trap = np.zeros((n_cells, n_E))
        self.trap_pit = 1. - rmaj_min[rho_lab[:]]/self.r2d[:]
        for jcell in range(n_cells):
            ind_trap = (self.a_d**2 <= self.trap_pit[jcell])
            ind_pass = (self.a_d**2 >  self.trap_pit[jcell])
            fbm_trap[jcell, :] = np.sum(self.fbm[jcell, ind_trap, :], axis=0)
            fbm_pass[jcell, :] = np.sum(self.fbm[jcell, ind_pass, :], axis=0)
        fbm_trap *= dpit
        fbm_pass *= dpit

# Integrals

        int_en = np.tensordot(self.fbm, dE, axes=(2, 0))
        self.int_en_pass = np.tensordot(fbm_pass, dE, axes=(1, 0))
        self.int_en_trap = np.tensordot(fbm_trap, dE, axes=(1, 0))
        int_en_pit = np.average(int_en, axis=1) # Would be np.sum(int_en)*dpit, but dpit = 1/npit

        self.dens_zone = np.zeros((n_zones, n_pit, n_E))
        self.dens_vol  = np.tensordot(self.fbm, bmvol, axes=(0,0))/np.sum(bmvol)
        vol_zone = np.zeros(n_zones)
        for jcell in range(n_cells):
            jrho = rho_lab[jcell]
            self.dens_zone[jrho, :, :] += self.fbm[jcell, :, :]*bmvol[jcell]
            vol_zone[jrho] += bmvol[jcell]

        for jrho in range(n_zones):
            self.dens_zone[jrho] *= 1/vol_zone[jrho]

        bdtmp = np.tensordot(self.dens_zone, dE, axes=(2,0))
        self.bdens = dpit*np.sum(bdtmp, axis=1)

        part_tot = np.sum(self.bdens*vol_zone)
        pass_tot = np.sum(self.int_en_pass*bmvol)
        trap_tot = np.sum(self.int_en_trap*bmvol)
        print('Passing #%12.4e  Trapped #%12.4e    Total #%12.4e' %(pass_tot, trap_tot, part_tot))
        print('Volume averaged fast ion density #12.4e m^-3' %(part_tot/vol))

        frac_trap = self.int_en_trap/(self.int_en_pass + self.int_en_trap)

        self.idlbl = '%s  Time = %8.4f s' %(runid, self.time)
        print self.idlbl

        self.thsurf = cv['THSURF'][:]
        self.xsurf  = cv['XSURF'][:]
        self.rsurf  = cv['RSURF'][:]
        self.zsurf  = cv['ZSURF'][:]

# MC cells: grid bars

        self.rbar = []
        self.zbar = []
        for jrho in range(n_zones):
            x_dist = (self.xsurf - self.rho_b[jrho])**2
            irho = np.argmin(x_dist)
            for jthe in range(self.the_size[jrho]):
                ythe = self.thb_grid2d[jrho,jthe]
                ithe = np.min(np.where(self.thsurf > ythe))
                r = np.zeros(2)
                z = np.zeros(2)
                if (jthe == 0):
                    r = self.rsurf[[irho, irho+2], ithe]
                    z = self.zsurf[[irho, irho+2], ithe]
                else:
                    th_ref = self.thsurf[ithe-1: ithe+1]
                    r[0] = np.interp(ythe, th_ref, self.rsurf[irho, ithe-1: ithe+1])
                    z[0] = np.interp(ythe, th_ref, self.zsurf[irho, ithe-1: ithe+1])
                    krho = irho + 2
                    r[1] = np.interp(ythe, th_ref, self.rsurf[krho, ithe-1: ithe+1])
                    z[1] = np.interp(ythe, th_ref, self.zsurf[krho, ithe-1: ithe+1])
                self.rbar.append(r)
                self.zbar.append(z)

# Select magnetic surfaces to plot

        self.ind_zone = []
        for rhob in self.rho_b:
            x_dist = (self.xsurf - rhob)**2
            irho = np.argmin(x_dist)
            self.ind_zone.append(irho)

# FBM plot
        print('Plot 2D fast ions density')
        xgrid = np.linspace(self.r2d.min(), self.r2d.max(), 1000)
        ygrid = np.linspace(self.z2d.min(), self.z2d.max(), 1000)
        self.x_grid, self.y_grid = np.meshgrid(xgrid, ygrid)
        self.zfbm  = griddata(self.r2d, self.z2d, int_en_pit, self.x_grid, self.y_grid)
        self.ztrap = griddata(self.r2d, self.z2d, frac_trap , self.x_grid, self.y_grid)
        self.cb_ax = self.fig_pol.add_axes([0.82, 0.08, 0.09, 0.84])

        self.cid = self.canvas.mpl_connect('button_press_event', self.my_call)

        self.reset_plot()
        self.plot_fbm()


    def plot_btneut(self):
        self.plot_neut()


    def plot_bbneut(self):
        self.plot_neut(type='bb')


    def plot_neut(self, type='bt'):

        try:
            tr_file = self.fname[:-9] + '.CDF'
        except:
            showerror("Error", self.err_msg)
            return

        try:
            self.tit.set_text('Beam-target neutron emissivity 1/(s cm**3) %s' %self.idlbl)
            print('%s' %tr_file)
            cdf = Dataset(tr_file, 'r', format='NETCDF4')
            cv = cdf.variables
            if type == 'bt':
                btneut = cv['BTNT2_DD'][:]
            else:
                btneut = cv['BBNT2_DD'][:]
            cdftime = cv['TIME3'][:]
# Select time
            tdist = (cdftime - self.time)**2
            jtclose = np.argmin(tdist)
            print('Time in CDF and cdf: %8.4f, %8.4f s' %(cdftime[jtclose], self.time) )
            self.bt    = griddata(self.r2d, self.z2d, btneut[jtclose,:], self.x_grid, self.y_grid)
            zmin = np.min(self.bt)
            zmax = np.max(self.bt)
            n_levels = 11
            bounds = np.linspace(zmin, zmax, n_levels)
            ctr3 = self.ax3.contourf(self.x_grid, self.y_grid, self.bt, bounds)
            norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
            mpl.colorbar.ColorbarBase(self.cb_ax,norm=norm,boundaries=bounds,ticks=bounds)
            self.canvas.draw()
        except:
            showerror("Error", 'CDF file\n%s\nnot found' %tr_file)


    def plot_fbm(self):

        try:
            self.tit.set_text('Fast ion density, integrated over E and p.a., %s' %self.idlbl)
            zmin = 0
            zmax = np.max(self.zfbm)
            n_levels = 11
            bounds = np.linspace(zmin, zmax, n_levels)
            ctr3 = self.ax3.contourf(self.x_grid, self.y_grid, self.zfbm, bounds)
            norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
            mpl.colorbar.ColorbarBase(self.cb_ax,norm=norm,boundaries=bounds,ticks=bounds)
            self.canvas.draw()
        except:
            showerror("Error", self.err_msg)


    def plot_trap(self):

        try:
            self.tit.set_text('Trapped fast ion fraction %s' %self.idlbl)
            zmin = 0
            zmax = 1
            n_levels = 11
            bounds = np.linspace(zmin, zmax, n_levels)
            ctr3 = self.ax3.contourf(self.x_grid, self.y_grid, self.ztrap)
            norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
            mpl.colorbar.ColorbarBase(self.cb_ax,norm=norm,boundaries=bounds,ticks=bounds)
            self.canvas.draw()
        except:
            showerror("Error", self.err_msg)


    def my_call(self, event):

        if event.button in (2,3):
            thint=self.th_int.get()
            rhoint=self.rho_int.get()
            logscale=self.log_scale.get()
            dist2 = (self.r2d-event.xdata)**2 + (self.z2d-event.ydata)**2
            jcell = np.argmin(dist2)
            if rhoint:
                self.line1.set_data(self.r2d, self.z2d)
            else:
                if thint:
                    ind=np.where(self.x2d == self.x2d[jcell])
                    self.line1.set_data(self.r2d[ind], self.z2d[ind])
                else:
                    self.line1.set_data(self.r2d[jcell], self.z2d[jcell])

            self.canvas.draw()
            self.plot_cell(jcell,thint=thint,rhoint=rhoint,logscale=logscale)


    def plot_bdens(self):

        try:
            plt.figure(2,figsize=(8,6))
            plt.clf()
            plt.plot(self.rho_grid,self.bdens)
            plt.xlabel(r'$\rho_{tor}$',fontsize=self.fsize)
            plt.ylabel('BDENS [1/cm**3]',fontsize=self.fsize)
            plt.title('Summary '+self.idlbl,fontsize=self.fsize)
            plt.show()
        except:
            showerror("Error", self.err_msg)


    def plot_birth(self):
#        try:
            self.plt_birth()
#        except:
#            showerror("Error", "Please browse a FBM file from the 'File' menu")


    def plt_birth(self):


        fbm_file = self.fname
        fbmfile  = fbm_file.split('/')[-1]
        runid    = fbmfile[:8]
        t_id     = fbmfile[12:13]
        birth_file =  '%s_birth.cdf%s' %(fbm_file[:-9], t_id)
        print('Reading %s' %birth_file)

        if not os.path.isfile(birth_file):
            print('%s not found' %birth_file)
            sys.exit()
        cdf = Dataset(birth_file, 'r', format='NETCDF4')
        cv = cdf.variables

        for key in cv.iterkeys():
            print key

        print('PLOT_BIRTH')
        print(birth_file)

        comp = ('full', 'half', '1/3')

        mcl = cv['mclabel'][:]
        mc_label = "".join(mcl[0]).strip()

# Read data from CDF

        cdf_d = {}
        sig = ('r', 'z', 'einj', 'xksid', 'zeta', 'time', 'ib')
        for lbl in sig:
            sigid = 'bs_%s_%s' %(lbl, mc_label)
            cdf_d[lbl]=cdf.variables[sigid][:]

        Rj = cdf_d['r']
        zj = cdf_d['z']

        l = cdf_d['ib'].tolist()
        n_birth = len(l)
        print('# of MC particles: %d' %n_birth)
        l.sort()
        unique = [x for i, x in enumerate(l) if not i or x != l[i-1]]
        src_arr = np.array(unique)
        print('Sources: ', src_arr)
        n_src = len(src_arr)

        Efull = {}
        for jsrc in src_arr:
            index = np.where(cdf_d['ib'] == jsrc )
            Efull[jsrc] = np.max(cdf_d['einj'][index])

# Determine 1/2 and 1/3 fractions

        mix_lbl = np.zeros(n_birth, dtype='|S4')
        mix_lbl[:] = comp[0]
        for jpart in range(n_birth):
            if cdf_d['einj'][jpart] <= 0.51*Efull[cdf_d['ib'][jpart]]:
                mix_lbl[jpart] = comp[1]
            if cdf_d['einj'][jpart] <= 0.34*Efull[cdf_d['ib'][jpart]]:
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
        deltaR = R_grid[1]-R_grid[0]
        deltaz = z_grid[1]-z_grid[0]
        dep_matrix={}
        for jsrc in src_arr:
            dep_matrix[jsrc] = {}
            for jcomp in comp:
                dep_matrix[jsrc][jcomp] = np.zeros((nr, nz))

        for jpart in range(n_birth):
            jsrc = cdf_d['ib'][jpart]
            jcomp = mix_lbl[jpart]
            Rpart = Rj[jpart]
            zpart = zj[jpart]
            if (Rpart > Rmin) and (Rpart < Rmax) and (zpart > zmin) and (zpart < zmax):
                jr = int((Rpart - Rmin)/deltaR)
                jz = int((zpart - zmin)/deltaz)
                dep_matrix[jsrc][jcomp][jr, jz] += 1

        dep_R = {}
        res_R = {}
        for jsrc in src_arr:
            dep_R[jsrc] = {}
            res_R[jsrc] = {}
            for jcomp in comp:
                res_R[jsrc][jcomp] = np.zeros(nr)
                dep_R[jsrc][jcomp] = np.sum(dep_matrix[jsrc][jcomp], axis=1) # z-sum
                res_R[jsrc][jcomp][nr-1] = np.sum(dep_R[jsrc][jcomp])
                for jr in range(nr-2, -1, -1):
                    res_R[jsrc][jcomp][jr] = res_R[jsrc][jcomp][jr+1] - dep_R[jsrc][jcomp][jr]

# NBI geometry for plots
        print('RUNID = %s' %runid)
        aug = nbi_geo.AUG_TR(0, runid=runid, raug=False)

#=======
# Plots
#=======

# Above view

        ntheta = 101
        Rmaj = self.rsurf[0,0]
        Rtor_in  = np.min(self.rsurf)
        Rtor_out = np.max(self.rsurf)
        phi_tor = np.linspace(0,2*np.pi,ntheta)
        phi_dep = np.radians(cdf_d['zeta'])
        cols = ('r', 'b', 'g', 'm', 'y')

# Component by component: v_par/v_perp

        nrows = 2
        ncols = n_comp + 1

        delta_pit = 0.2

        jfig = 12
        fig = {}
# One figure for each source

        cosp = np.cos(phi_tor)
        sinp = np.sin(phi_tor)

        xpol = (90,250)
        ypol = (-125,125)

        aug_geo = {}
        xend = np.array(4*[-60]+[40,250,250,40])

        for jsrc in src_arr:
            jnb = jsrc - 1
            aug_geo['xend']  = xend[jnb]
            aug_geo['x0']    = aug.xsrc[jnb]
            aug_geo['y0']    = aug.ysrc[jnb]
            aug_geo['z0']    = aug.xybsca[jnb]
            aug_geo['theta'] = aug.theta_los[jnb]
            aug_geo['phi']   = aug.phi_los[jnb]
            nbi_los = los.PHILOS(aug_geo)
            xlin = nbi_los.xline
            ylin = nbi_los.yline
            rlin = nbi_los.rline
            zlin = nbi_los.zline

            lbl = 'NBI #%d' %jsrc

            fig[jfig] = plt.figure(jfig, figsize=(5*ncols, 5.5*nrows))
            fig[jfig].clf()
            fig[jfig].subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.95)

            jsplot = 1

            fig[jfig].text(0.5, 0.975, lbl               , ha='center')
            fig[jfig].text(0.5, 0.95 , 'Top view'        , ha='center')
            fig[jfig].text(0.5, 0.49 , 'Poloidal section', ha='center')

# Overplot 3 species mix
# Above view
            axtop = fig[jfig].add_subplot(nrows, ncols, jsplot, aspect='equal')
            axtop.set_title('All energy components', fontsize=self.fsize)
            axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
            axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
            axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')
            axtop.plot(xlin, ylin, 'g-', linewidth=2.5)
            try:
                for key, tor_pl in self.tor_d.iteritems():
                    axtop.plot(100*tor_pl.x, 100*tor_pl.y, 'b-')
            except:
                print('')

# Poloidal section

            axpol = fig[jfig].add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')
            axpol.set_title('All energy components', fontsize=self.fsize)
            axpol.set_xlabel(self.Rlbl, fontsize=self.fsize)
            axpol.set_ylabel(self.zlbl, fontsize=self.fsize)
            axpol.set_xlim(xpol)
            axpol.set_ylim(ypol)
            axpol.plot(rlin, zlin, 'g-', linewidth=2.5) 

            try:
                for key in self.kgc.gc_x.iterkeys():
                    axpol.plot(100*self.kgc.gc_x[key], 100*self.kgc.gc_y[key], 'b-')
            except:
                print('')

# Overlay cells
            for irho in self.ind_zone:
                axpol.plot(self.rsurf[irho, :], self.zsurf[irho, :], color=self.MCgrid_col)
            for jbar, myr in enumerate(self.rbar):
                axpol.plot(myr, self.zbar[jbar], color=self.MCgrid_col)

            jcol=0
            for jcomp in comp:
                ind = np.where((mix_lbl == jcomp) & (cdf_d['ib'] == jsrc))
                axtop.plot(Rj[ind]*np.cos(phi_dep[ind]), Rj[ind]*np.sin(phi_dep[ind]), cols[jcol]+'o', label=jcomp)
                axpol.plot(Rj[ind], zj[ind], '%so' %cols[jcol], label=jcomp)
                jcol += 1
            axtop.legend()
            axpol.legend()
 
            jcol = 0
            jsplot = 2

# For each species overplot 5 pitch angle range
            for jcomp in comp:
                axtop = fig[jfig].add_subplot(nrows, ncols, jsplot      , aspect='equal')
                axpol = fig[jfig].add_subplot(nrows, ncols, jsplot+ncols, aspect='equal')

                axtop.set_title(jcomp+' energy components', fontsize=self.fsize)
                axtop.plot(Rtor_in *cosp, Rtor_in *sinp, 'r-')
                axtop.plot(Rtor_out*cosp, Rtor_out*sinp, 'r-')
                axtop.plot(Rmaj*cosp, Rmaj*sinp, 'r--')
                axtop.plot(xlin, ylin, 'g-', linewidth=2.5) 
                try:
                    for key, tor_pl in self.tor_d.iteritems():
                        axtop.plot(100*tor_pl.x, 100*tor_pl.y, 'b-')
                except:
                    print('')

                axpol.set_title('%s energy components' %jcomp, fontsize=self.fsize)
                axpol.set_xlabel(self.Rlbl, fontsize=self.fsize)
                axpol.set_ylabel(self.zlbl, fontsize=self.fsize)
                axpol.set_xlim(xpol)
                axpol.set_ylim(ypol)
                axpol.plot(rlin, zlin, 'g-', linewidth=2.5) 
                try:
                    for key in self.kgc.gc_x.iterkeys():
                        axpol.plot(100*self.kgc.gc_x[key], 100*self.kgc.gc_y[key], 'b-')
                except:
                    print('')

# Overlay cells
                for irho in self.ind_zone:
                    axpol.plot(self.rsurf[irho, :], self.zsurf[irho,:], color=self.MCgrid_col)
                for jbar, myr in enumerate(self.rbar):
                    axpol.plot(myr, self.zbar[jbar], color=self.MCgrid_col)

                p1 = 0
                jcol = 0
                while p1 <= 1-delta_pit:
                    p2 = p1 + delta_pit
                    ind = np.where((mix_lbl == jcomp) & (cdf_d['ib'] == jsrc) & \
                                   (cdf_d['xksid'] > p1) & (cdf_d['xksid'] < p2)) 
                    axtop.plot(Rj[ind]*np.cos(phi_dep[ind]), \
                               Rj[ind]*np.sin(phi_dep[ind]), '%so' %cols[jcol], \
                               label='%3.1f < p.a. < %3.1f' %(p1, p2))
                    axpol.plot(Rj[ind], zj[ind], '%so' %cols[jcol], \
                               label='%3.1f < p.a. < %3.1f' %(p1, p2))
                    p1 += delta_pit
                    jcol += 1
#                axtop.legend()
#                axpol.legend()
                jsplot += 1

            fig[jfig].canvas.mpl_connect('button_press_event', fconf.on_click)
            jfig += 1

#-------------------------        
# Deposition & attenuation
#-------------------------

        nrows = 2
        ncols = 3

        xgrid, ygrid = np.meshgrid(R_grid, z_grid)
        jfig += 1
        for jsrc in src_arr:
            jnb = jsrc - 1
            aug_geo['xend']  = xend[jnb]
            aug_geo['x0']    = aug.xsrc[jnb]
            aug_geo['y0']    = aug.ysrc[jnb]
            aug_geo['z0']    = aug.xybsca[jnb]
            aug_geo['theta'] = aug.theta_los[jnb]
            aug_geo['phi']   = aug.phi_los[jnb]
            nbi_los = los.PHILOS(aug_geo)
            fig[jfig] = plt.figure(jfig, figsize=(6*ncols, 6*nrows))
            fig[jfig].clf()
            fig[jfig].subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.95, \
                                      wspace=0.1, hspace=0.1)
# 2D deposition, poloidal section
            jsplot = 1
            for jcomp in comp:
                zgrid = dep_matrix[jsrc][jcomp].T
                ind = np.where(zgrid == 0)
                zgrid[ind] = None
                axpol=fig[jfig].add_subplot(nrows, ncols, jsplot, aspect='equal')
                axpol.set_xlabel(self.Rlbl, fontsize=self.fsize)
                axpol.set_ylabel(self.zlbl, fontsize=self.fsize)
                axpol.set_xlim(xpol)
                axpol.set_ylim(ypol)
                axpol.plot(nbi_los.rline, nbi_los.zline, 'g-', linewidth=2.5) 
                try:
                    for key in self.kgc.gc_x.iterkeys():
                        axpol.plot(100*self.kgc.gc_x[key], 100*self.kgc.gc_y[key], 'b-')
                except:
                    print('')

                ctr = axpol.contourf(xgrid, ygrid, zgrid)
                fig[jfig].colorbar(ctr, shrink=0.8, aspect=10)

# Overlay cells
                for irho in self.ind_zone:
                    axpol.plot(self.rsurf[irho, :], self.zsurf[irho, :], color=self.MCgrid_col)
                for jbar, myr in enumerate(self.rbar):
                    axpol.plot(myr, self.zbar[jbar], color=self.MCgrid_col)

                jsplot += 1

            for jcomp in comp:
                ax2 = fig[jfig].add_subplot(nrows, ncols, jsplot)
                ax2.set_xlabel(self.Rlbl, fontsize=self.fsize)
                ax2.set_ylabel('NBI attenuation', fontsize=self.fsize)
                ax2.plot(R_grid, res_R[jsrc][jcomp])
                jsplot += 1
            fig[jfig].canvas.mpl_connect('button_press_event', fconf.on_click)
            jfig += 1
            
        plt.show()


    def plot_cell(self, jcell, rhoint=0, thint=0, logscale=0):

# "None" better then "0" for clarity in contour plot
        fnc_d = {'cell':self.fbm, 'zone':self.dens_zone, 'vol':self.dens_vol}
        title_d = {'vol': r'$\rho_{tor}$ averaged, $\theta$ averaged', \
                   'zone': r'$\rho_{tor}$'+r' = %8.4f, $\theta$ averaged' %self.x2d[jcell], \
                   'cell': r'$\rho_{tor} = $ %8.4f $\theta = $%8.4f deg' \
                   %(self.x2d[jcell], np.degrees(self.th2d[jcell])) }

        jrho = np.where(self.rho_grid == self.x2d[jcell])[0][0]

        if rhoint:
            lbl = 'vol'
            zarr = fnc_d[lbl]
        else:
            if thint:
                lbl = 'zone'
                zarr = fnc_d[lbl][jrho]
            else:
                lbl = 'cell'
                zarr = fnc_d[lbl][jcell]

        zmin = 1e4
        indzero = np.where(zarr <= zmin)
        levels = np.linspace(zmin, np.max(zarr), 10)
        zarr[indzero] = np.nan
        if logscale:
            zarr = np.log10(zarr)
            levels = np.linspace(np.log10(zmin), 8.5, 10)

        print levels
        fig_rho = plt.figure(1)
        fig_rho.clf()
        ax = fig_rho.add_subplot(1, 1 ,1)
        ax.set_xlabel('Energy [keV]', fontsize=self.fsize)
        ax.set_ylabel('Pitch angle' , fontsize=self.fsize)
        ax.set_xlim((0, 500))
        ax.set_title(title_d[lbl], fontsize=self.fsize)
        ctr = ax.contourf(self.e_d, self.a_d, zarr, levels)
        fig_rho.colorbar(ctr, aspect=10)
        ax.plot([0, np.max(self.e_d)], [0, 0], 'k-')
        if lbl == 'cell':
            ax.plot([0, np.max(self.e_d)], [ self.trap_pit[jcell],  self.trap_pit[jcell]], 'g-')
            ax.plot([0, np.max(self.e_d)], [-self.trap_pit[jcell], -self.trap_pit[jcell]], 'g-')
        plt.show()


if __name__ == "__main__":

    FBM()
