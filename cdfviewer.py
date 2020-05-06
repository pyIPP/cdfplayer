#!/usr/bin/env python

import sys, os, shutil

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

import numpy as np
import ufiles, tr_read_ctr, tkhyper, ctr2rz
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
try:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk as nt2tk
except:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as nt2tk

from matplotlib.figure import Figure
import matplotlib.pylab as plt
from scipy.io import netcdf


class VIEWER:


    timeout = 1e-10
    dir_set = '%s/python/cdfp' %os.getenv('HOME')

    inNextOrPrev = False


    def __init__(self, f_cdf_in, list_file='', tok='AUGD'):


        locdir = os.path.dirname(os.path.realpath(__file__))
        self.tok = tok

        self.runid     = f_cdf_in[-12:-4]
        self.cdf_file  = f_cdf_in
        self.list_file = list_file

        self.jtab = 0

# Widget frame

        print('Setting GUI frame')
        self.viewframe = tk.Tk()
        xmax = self.viewframe.winfo_screenwidth()
        ymax = self.viewframe.winfo_screenheight()
        width  = min(1600, int(0.95*xmax))
        height = min(960, int(0.88*ymax)) 
        self.viewframe.title('NetCDF viewer')
        self.viewframe.geometry('%dx%d' %(width, height))
        self.viewframe.option_add("*Font", "Helvetica")
# Menu bar

        menubar = tk.Menu(self.viewframe)

        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Load run...", command=self.callcdf)
        filemenu.add_command(label="Display variables...", command=self.display_vars)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)

        editmenu = tk.Menu(menubar, tearoff=0)
        editmenu.add_command(label="Load setup...", command=self.callset)
        editmenu.add_command(label="Save setup...", command=self.callsave)
        editmenu.add_command(label="Edit setup...", command=self.editlist)

        savemenu = tk.Menu(menubar, tearoff=0)
        savemenu.add_command(label="Save 1D u-file(t)"     , command=self.uf1t)
        savemenu.add_command(label="Save 1D u-file(rho)"   , command=self.uf1r)
        savemenu.add_command(label="Save 2D u-file(t, rho)", command=self.uf2)
        savemenu.add_command(label="Save 2D as .avi"       , command=self.smov)
        if self.tok == 'AUGD':
            savemenu.add_command(label="Save eqdsk"                 , command=self.eqdisk)
            savemenu.add_command(label="Save shotfile (all signals)", command=self.tr2sf)
            savemenu.add_command(label="Save equil shotfile"        , command=self.cdf2tre)

        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.about)

        menubar.add_cascade(label="File"  ,  menu=filemenu)
        menubar.add_cascade(label="Setup" ,  menu=editmenu)
        menubar.add_cascade(label="Output",  menu=savemenu)
        menubar.add_cascade(label="Help"  ,  menu=helpmenu)

        self.viewframe.config(menu = menubar)

# Frames in main widet

        canv_frame = ttk.Frame(self.viewframe, height=height-60)
        canv_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        canv_frame.pack_propagate(0)

        toolframe = ttk.Frame(self.viewframe, height=45)
        toolframe.pack(side=tk.BOTTOM, fill=tk.X)
        toolframe.pack_propagate(0)

# Plots

        nbpol = ttk.Notebook(canv_frame, name='nbpol')
        nbpol.pack(side=tk.TOP, fill=tk.X)
        nbpol.bind('<Button-1>', self.on_click)

        frame_eq = ttk.Frame(nbpol)
        frame_1d = ttk.Frame(nbpol)
        frame_2d = ttk.Frame(nbpol)

        nbpol.add(frame_eq, text='Equilibirum')
        nbpol.add(frame_1d, text='Scalars    ')
        nbpol.add(frame_2d, text='Profiles   ')

        self.fig_eq = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_1d = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_2d = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_eq.subplots_adjust(left=0.01, bottom=0.1, right=0.95, top=0.95)
        self.fig_1d.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.92, \
            wspace=0.3, hspace=0.3)
        self.fig_2d.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.92, \
            wspace=0.3, hspace=0.3)
        canv_eq = FigureCanvasTkAgg(self.fig_eq, master=frame_eq)
        canv_1d = FigureCanvasTkAgg(self.fig_1d, master=frame_1d)
        canv_2d = FigureCanvasTkAgg(self.fig_2d, master=frame_2d)
        for canv in canv_eq, canv_1d, canv_2d:
            canv._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Navigation toolbars

        toolbar = nt2tk( canv_eq, frame_eq)
        toolbar.update()
        toolbar = nt2tk( canv_1d, frame_1d)
        toolbar.update()
        toolbar = nt2tk( canv_2d, frame_2d)
        toolbar.update()

# Toolbar

        self.playfig  = tk.PhotoImage(file='%s/ButtonPlay.gif'  %locdir)
        self.pausefig = tk.PhotoImage(file='%s/ButtonPause.gif' %locdir)
        rewfig  = tk.PhotoImage(file='%s/ButtonBack.gif'     %locdir)
        prevfig = tk.PhotoImage(file='%s/ButtonPrevious.gif' %locdir)
        nextfig = tk.PhotoImage(file='%s/ButtonNext.gif'     %locdir)

        self.playbt = ttk.Button(toolframe, command=self.play , image=self.playfig)
        rewbt  = ttk.Button(toolframe, command=self.rewind, image=rewfig )
        prevbt = ttk.Button(toolframe, command=self.prev  , image=prevfig)
        nextbt = ttk.Button(toolframe, command=self.next  , image=nextfig)
        for but in self.playbt, rewbt, prevbt, nextbt:
            but.pack(side=tk.LEFT)

        self.crntsc = tk.Scale(toolframe, command=self.jump, \
                               orient=tk.HORIZONTAL, length=300)
        self.crntsc.pack(side=tk.LEFT)

# Set up plot
        self.load()

        self.viewframe.mainloop()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/aug/manuals/transp/output/cdfplayer.html">cdfplayer homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "340x60")


    def on_click(self, event):
        if event.widget.identify(event.x, event.y) == 'label':
            self.jtab = event.widget.index('@%d,%d' % (event.x, event.y))
        self.update_plot()


    def read(self, f_cdf):

        if not os.path.isfile(f_cdf):
            print('NetCDF file %s not found' %f_cdf)
            return

        if hasattr(self, 'surf'):
            del self.surf
        print('Reading profiles from %s' %f_cdf)
        self.cv = netcdf.netcdf_file(f_cdf, 'r', mmap=False).variables

        print('Reading equilibrium from %s' %f_cdf)
        rho = np.linspace(0.1, 1.0, 10)
        if 'TIME3' in self.cv.keys():
            self.surf = tr_read_ctr.TR_READ_CTR(f_cdf, rho=rho, nthe=201, endpoint=True)


    def set_plots(self, pr_list=None):

        self.mark1d  = {}
        self.equline = {}
        self.line1d  = {}
        self.line2d  = {}

        if os.path.isfile(self.list_file) and pr_list is None:
            f = open(self.list_file)
            pr_list = [s.strip() for s in f.readlines()]
            f.close()
            print('Using profile list from file %s' %self.list_file)
        elif pr_list is None:
            pr_list = ['CUR', 'TE', 'NE', 'V', 'BPOL', 'NEUTT', 'NIMP', 'PINJ', 'TAUEA']
 
        if 'TIME3' in self.cv.keys():
            self.time = self.cv['TIME3'].data
            self.rho  = self.cv['X'][0, :]
            self.rhob = self.cv['XB'][0, :]
        else: 
            self.time = self.cv['TIME'].data
            self.rho  = self.cv['XRHO'].data
        self.lbl  = {}
        self.dim  = {}
        self.unit = {}
        self.prof_list  = []
        self.trace_list = []
        for sig in pr_list:
            sig = sig.upper()
            if sig in self.cv.keys():
                self.lbl[sig]  = self.cv[sig].long_name[:-2]
                self.dim[sig]  = self.cv[sig].dimensions
                self.unit[sig] = self.cv[sig].units
                self.nt = self.cv[sig].data.shape[0]
                if self.cv[sig].data.ndim > 1:
                    self.prof_list.append(sig)
                else:
                    self.trace_list.append(sig)
            else:
                print('No variable %s in CDF file, ignoring!' %sig)

        if len(self.trace_list) > 15:
            print('Too many time traces, taking only the first 15')
            self.trace_list = self.trace_list[:15]
        if len(self.prof_list) > 15:
            print('Too many profiles, taking only the first 15')
            self.prof_list = self.prof_list[:15]
        nrow = 3
        ncol = 5

# Equilibrium

        if hasattr(self, 'surf'):
            self.R = self.surf.Rsurf
            self.Z = self.surf.Zsurf
        else:
            self.R = self.cv['Rsurf'].data
            self.Z = self.cv['Zsurf'].data
        self.n_rho = self.R.shape[1]
        print('Nrho_eq = %d' %self.n_rho)
        print(self.R.shape)
# Plots

        fsize = 12

        self.fig_eq.clf()
        self.fig_1d.clf()
        self.fig_2d.clf()

        self.fig_eq.text(0.5, 0.99, self.runid, ha='center', va='top')
        self.fig_1d.text(0.5, 0.99, self.runid, ha='center', va='top')
        self.fig_2d.text(0.5, 0.99, self.runid, ha='center', va='top')
        self.txt_eq = self.fig_eq.text(0.5, 0.97, '', ha='center', va='top')
        self.txt_1d = self.fig_1d.text(0.5, 0.97, '', ha='center', va='top')
        self.txt_2d = self.fig_2d.text(0.5, 0.97, '', ha='center', va='top')

        ax_eq = self.fig_eq.add_subplot(1, 1, 1, aspect='equal')
        for jrho in range(self.n_rho):
            self.equline[jrho],  = ax_eq.plot([], [], 'g-')
        self.equline['axis'], = ax_eq.plot([], [], 'g+')
        ax_eq.set_xlabel('R [cm]', fontsize=fsize)
        ax_eq.set_ylabel('z [cm]', fontsize=fsize)

# Plot AUG structures
#        nshot = int(self.runid[0:5])
        try:
            import get_gc
            if hasattr(self, 'surf'):
                fac=100.
            else:
                fac=1.
            gc_r, gc_z = get_gc.get_gc()
            for key in gc_r.keys():
                ax_eq.plot(fac*gc_r[key], fac*gc_z[key], 'b-')

        except:
            pass

        jplot = 1
        for trace in self.trace_list:
            ylab = '%s [%s]' %(trace, self.unit[trace].strip())

            ax_1d = self.fig_1d.add_subplot(nrow, ncol, jplot)
            ax_1d.set_xlabel('t [s]', fontsize=fsize)
            ax_1d.set_ylabel(ylab, fontsize=fsize)
            self.line1d[trace], = ax_1d.plot(self.time, self.cv[trace][:], 'b-')
            ax_1d.grid('on')
            self.mark1d[trace], = ax_1d.plot(self.time[0], self.cv[trace][0], 'go')
            jplot += 1
            ax_1d.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))
            ax_1d.yaxis.major.formatter._useMathText = True

        jplot = 1
        for prof in self.prof_list:
            if 'X' in self.dim[prof]:
                xrho = self.rho
                xlab = r'$\rho_{tor}$(X)'
            elif 'XB' in self.dim[prof]:
                xrho = self.rhob
                xlab = r'$\rho_{tor}$(XB)'
            elif 'XRHO' in self.dim[prof]:
                xrho = self.rho
                xlab = r'$\rho_{tor}$(XRHO)'
            ylab = '%s [%s]' %(prof, self.unit[prof].strip())

            ax_2d = self.fig_2d.add_subplot(nrow, ncol, jplot)
            ax_2d.set_xlim([0, 1.])
            ax_2d.set_ylim([0, 1.1*np.max(self.cv[prof].data)])
            ax_2d.set_xlabel(xlab, fontsize=fsize)
            ax_2d.set_ylabel(ylab, fontsize=fsize)
            self.line2d[prof], = ax_2d.plot(xrho, self.cv[prof][0, :], 'b-')
            ax_2d.grid('on')
            ax_2d.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))
            ax_2d.yaxis.major.formatter._useMathText = True
            jplot += 1

        self.update_plot()


    def eqdisk(self):

        import eqdsk, get_sf_grid

        cdf1d = ['PCUR', 'BZ', 'BZXR', 'RAXIS', 'YAXIS']
        cdf2d = ['PLFLX', 'Q', 'GFUN', 'PMHD_IN']
        coord = ['XB', 'X']

#==================
# Psi(R, z), j(R, z)
#==================

        cdf_d = {}
        cdf_d['n_the'] = '101'
        cdf_d['time']  = self.time[self.jt]

        eq = get_sf_grid.get_grid(self.runid)
        Rgrid = eq['Ri'][:, 0]
        zgrid = eq['Zj'][:, 0]
        rz = ctr2rz.CTR2RZ(self.cdf_file, it=self.jt)

        data = {}
        for lbl in cdf1d + cdf2d + coord:
            data[lbl] = self.cv[lbl].data[self.jt]

        data['fpol'] = -0.01*data['GFUN']*data['BZXR']

        xb_grid = np.append(0, data['XB'])
        n_xb = len(xb_grid)

        pfl  = np.append(0, data['PLFLX'])
        qpsi = np.interp(xb_grid, data['XB'], data['Q'])
        fpol = np.interp(xb_grid, data['XB'], data['fpol'])
        Pres = np.interp(xb_grid, data['X'] , data['PMHD_IN'])

# Derivatives d/dPsi

        dpsi  = np.gradient(pfl)
        dfpol = np.gradient(fpol)/dpsi
        dPres = np.gradient(Pres)/dpsi

        eqd_d = {}
        eqd_d['Rgrid']   = Rgrid
        eqd_d['zgrid']   = zgrid
        eqd_d['Raxis']   = 0.01*data['RAXIS']
        eqd_d['zaxis']   = 0.01*data['YAXIS']
        eqd_d['psi_ax']  = pfl[0]
        eqd_d['psi_bdy'] = pfl[-1]
        eqd_d['Ipl']     = data['PCUR']
        eqd_d['bcentr']  = data['BZ']
        eqd_d['rcentr']  = 0.01*data['BZXR']/eqd_d['bcentr']
        eqd_d['qpsi']    = qpsi
        eqd_d['PFL']     = pfl
        eqd_d['fpol']    = fpol
        eqd_d['pres']    = Pres
        eqd_d['fprim']   = dfpol
        eqd_d['pprim']   = dPres
        eqd_d['r_bdy']   = 0.01*self.R[self.jt, -1, :]
        eqd_d['z_bdy']   = 0.01*self.Z[self.jt, -1, :]
        eqd_d['psi']     = -rz.pfm/(2.*np.pi)

        f_eqdsk = self.cdf_file.replace('.CDF', '.eqdsk')
        eqdsk.write_eqdsk(eqd_d, f_eqdsk, self.time[self.jt], f_cdf=self.cdf_file)


    def uf1t(self):

        self.ff = False
        self.rw = False

        sshot = self.runid[0:5]
        tlbl = 'Time:'.ljust(20) + 's'
        for trace in self.trace_list:
            dlbl = trace.ljust(20) + self.unit[trace].strip()
            uf_d = { 'pre':'A', 'ext':trace, 'shot':sshot, \
                     'grid': {'X': {'lbl':tlbl, 'arr':self.time} }, \
                     'data': {'lbl': dlbl, 'arr': self.cv[trace][:]} }
            ufiles.WU(uf_d)


    def uf1r(self):

        self.ff = False
        self.rw = False

        sshot = self.runid[0:5]
        xlbl = 'rho_tor'
        for prof in self.prof_list:
            if 'XB' in self.dim[prof]:
                xrho = self.rhob
            else:
                xrho = self.rho
            dlbl = prof.ljust(20) + self.unit[prof].strip()
            uf_d = { 'pre': 'A', 'ext': prof, 'shot': sshot, \
                     'scal': [['Time:'.ljust(20) + 's', self.time[self.jt]]], 
                     'grid': {'X': {'lbl': xlbl, 'arr': xrho} }, \
                     'data': {'lbl': dlbl, 'arr': self.cv[prof][self.jt, :]} }
            ufiles.WU(uf_d)


    def uf2(self):

        self.ff = False
        self.rw = False
        
        sshot = self.runid[0:5]
        tlbl = 'Time'.ljust(20) + 'Seconds'
        xlbl = 'rho_tor'
        for prof in self.prof_list:
            if 'XB' in self.dim[prof]:
                xrho = self.rhob
            else:
                xrho = self.rho
            dlbl = prof.ljust(20) + self.unit[prof].rstrip()
            uf_d = { 'pre': 'A', 'ext': prof + '2', 'shot': sshot, \
                     'grid': {'X':{'lbl': tlbl, 'arr': self.time}, \
                              'Y':{'lbl': xlbl, 'arr': xrho} }, \
                     'data': {'lbl': dlbl, 'arr': self.cv[prof].data} }
            ufiles.WU(uf_d)


    def cdf2tre(self):

        import cdf2sf
        cdf2sf.cdf2tre(self.runid)


    def tr2sf(self):

        import cdf2sf
        cdf2sf.cdf2tra(self.runid)


    def smov(self):

        self.ff = False
        self.rw = False

# Store png's
        order = int(np.log10(self.nt)) + 1 # To have nice alphabetic sorted files
        fmt = '%' + str(order) + 'd'
        for jtim in range(self.nt):
            self.jt = jtim
            self.update_plot()
            jstr = (fmt % self.jt).replace(' ', '0')
            fout = '%s-%s.png' %(self.runid, jstr)
            print('Stored %s' %fout)
            self.viewfig.savefig(fout, dpi=100)

# Concatenate png's to avi
        out_mov = '%s.avi' %self.runid
        command = ('mencoder', 
                   'mf://*.png', 
                   '-mf', 
                   'type=png:w=800:h=600:fps=25', 
                   '-ovc', 
                   'lavc', 
                   '-lavcopts', 
                   'vcodec=mpeg4', 
                   '-oac', 
                   'copy', 
                   '-o', 
                   out_mov)

        encoder = '/afs/@cell/common/soft/visualization/mencoder/svn-2011-05-17/@sys/bin/mencoder'
        os.spawnvp(os.P_WAIT, encoder, command)
        os.system('rm *.png')
        print('Stored movie %s' %out_mov)


    def update_plot(self, reset=True):

#        if hasattr(self, 'crntsc') and reset:
        if reset:
            self.crntsc.set(float(self.jt)/(self.nt - 1)*100.)

        strtim = 't =%6.3f s' %self.time[self.jt]
        if self.jtab == 0:
            for jrho in range(self.n_rho):
                self.equline[jrho].set_data(self.R[self.jt, jrho, :], self.Z[self.jt, jrho, :])
            self.txt_eq.set_text(strtim)
# Plot mag. axis:
            if hasattr(self, 'surf'):
                self.equline['axis'].set_data(self.surf.r_mag[self.jt], self.surf.z_mag[self.jt])
            self.fig_eq.canvas.draw()

        elif self.jtab == 1:
            for trace in self.trace_list:
                self.mark1d[trace].set_data(self.time[self.jt], self.cv[trace][self.jt])
            self.txt_1d.set_text(strtim)
            self.fig_1d.canvas.draw()

        elif self.jtab == 2:
            for prof in self.prof_list:
                self.line2d[prof].set_ydata(self.cv[prof][self.jt, :])
            self.txt_2d.set_text(strtim)
            self.fig_2d.canvas.draw()


    def prev(self):

        self.ff = False
        self.rw = False
        if self.jt > 0:
            self.jt -= 1
            self.inNextOrPrev = True
            self.update_plot()
            self.inNextOrPrev = False


    def next(self):

        self.ff = False
        self.rw = False
        if self.jt < self.nt - 1:
            self.jt += 1
            self.inNextOrPrev = True
            self.update_plot()
            self.inNextOrPrev = False


    def jump(self, arg):

        if self.inNextOrPrev:
            return
        self.jt = int(float(arg)/100. * (self.nt-1))
        self.update_plot(reset=False)


    def pause(self):

        self.ff = False
        self.rw = False
        self.playbt['image']   = self.playfig
        self.playbt['command'] = self.play
        self.update_plot()


    def play(self):

        self.ff = True
        self.rw = False
        self.playbt['image']   = self.pausefig
        self.playbt['command'] = self.pause
        while self.ff and self.jt < self.nt - 1:
            self.jt += 1
            self.update_plot()
            self.fig_eq.canvas.start_event_loop(self.timeout)


    def rewind(self):

        self.ff = False
        self.rw = True
        while self.rw and self.jt > 0:
            self.jt -= 1
            self.update_plot()
            self.fig_eq.canvas.start_event_loop(self.timeout)


    def callcdf(self):

        if hasattr(self, 'cdf_file'):
            dir_init = os.path.dirname(self.cdf_file)
        else:
            dir_init = '%s/tr_client/%s' %(os.getenv('HOME'), self.tok)
        self.cdf_file = tkfd.askopenfilename(initialdir=dir_init, filetypes=[("All formats", "*.CDF")])

        self.load()


    def callset(self):

        self.list_file = tkfd.askopenfilename(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        self.set_plots()


    def callsave(self):

        f_out = tkfd.asksaveasfile(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        out_txt = ''
        for prof in self.prof_list:
            out_txt += prof + '\n'
        for sig in self.trace_list:
            out_txt += sig + '\n'
        f_out.write(out_txt)
        f_out.close()


    def editlist(self, profiles=None):

        sig_list = self.prof_list + self.trace_list
        self.topl = tk.Toplevel(self.viewframe)
        self.topl.resizable(0, 0)
        ttk.Label(self.topl, text="Profiles to be viewed").grid(row=0, columnspan=3)
        self.txt = tk.Text(self.topl, width=40, height=10)
        self.txt.grid(row=1, columnspan=2)
        self.txt.insert('end', "\n".join(sig_list if profiles is None else profiles))
        scr = tk.Scrollbar(self.topl)
        scr.config(command=self.txt.yview)
        scr.grid(row=1, column=2, sticky=tk.N+tk.S)
        self.txt.config(yscrollcommand=scr.set)
        b = ttk.Button(self.topl, text="OK", command=self.dialog_ok)
        b.grid(row=2, column=0)
        b = ttk.Button(self.topl, text="Cancel", command=self.topl.destroy)
        b.grid(row=2, column=1)


    def dialog_ok(self):

        newsigs = self.txt.get("0.0", "end-1c").split()
        oldprofiles = self.prof_list
        oldtraces   = self.trace_list
        try:
            self.set_plots(pr_list=newsigs)
        except ValueError:
            print('Error!')
            self.prof_list  = oldprofiles
            self.trace_list = oldtraces
        self.topl.destroy()


    def load(self):

        self.jt = 0
        self.runid = self.cdf_file[-12:-4]
        self.read(self.cdf_file)
        self.set_plots()


    def display_vars(self):

        dispframe = tk.Toplevel()
#        dispframe.geometry('240x200')
        dispframe.title('Advanced search')
        entries = ['Var_name', 'Unit', 'Description']
        self.disp_d = {}

        nrow = 0
        for key in entries:
            lbl = ttk.Label(dispframe, text=key)
            var = ttk.Entry(dispframe, width=12)
            lbl.grid(row=nrow, column=0, sticky=tk.W+tk.E)
            var.grid(row=nrow, column=1, sticky=tk.W+tk.E)
            var.insert(0, '')
            self.disp_d[key] = var
            nrow += 1
        dim = tk.IntVar()
        dim.set(3)
        rb1 = ttk.Radiobutton(dispframe, text='Scalar' , variable=dim, val=1)
        rb2 = ttk.Radiobutton(dispframe, text='Profile', variable=dim, val=2)
        rb3 = ttk.Radiobutton(dispframe, text='All', variable=dim, val=3)
        rb1.grid(row=nrow, column=0, sticky=tk.W+tk.E)
        rb2.grid(row=nrow, column=1, sticky=tk.W+tk.E)
        rb3.grid(row=nrow, column=2, sticky=tk.W+tk.E)
        self.disp_d['dim'] = dim
        nrow += 1
        but1 = ttk.Button(dispframe, text='Print variable', command=self.print_vars)
        but2 = ttk.Button(dispframe, text='Close', command=dispframe.destroy)
        but1.grid(row=nrow, column=0)
        but2.grid(row=nrow, column=1)

        dispframe.mainloop()


    def print_vars(self):

        skey  = self.disp_d['Var_name'].get().upper().strip()
        sunit = self.disp_d['Unit'].get().upper().strip()
        sdesc = self.disp_d['Description'].get().upper().strip()
        dim = self.disp_d['dim'].get()
        for key, val in self.cv.items():
            unit  = val.units.strip()
            descr = val.long_name.strip()
            key_flag  = (skey  == '') or (skey  in key)
            unit_flag = (sunit == '') or (sunit in unit)
            desc_flag = (sdesc == '') or (sdesc in descr)
            if key_flag and unit_flag and desc_flag:
                if val.data.ndim == dim or dim == 3:
                    print(key.ljust(10), unit.ljust(16), descr.ljust(20))
        print('')


if __name__ == "__main__":


    cdf_file  = '/afs/ipp/home/g/git/tr_client/AUGD/23076/W01/23076W01.CDF'
    list_file = ''

    if len(sys.argv) > 2:
        f = sys.argv[2]
        if f in ('-h', '--help', '-?'):
            print('Usage: cdfplayer [<cdf file> [<profile list file>]]')
            sys.exit()
        else:
            import tr_path
            tr = tr_path.TR_PATH(f)
            f = tr.fcdf
        if os.path.isfile(f):
            cdf_file = f
    if len(sys.argv) > 3:
            f = sys.argv[3]
            if os.path.isfile(f):
                list_file = f
            else:
                f_path = os.getenv('HOME') + '/python/cdfp/' + f
                if os.path.isfile(f_path):
                    list_file = f_path
    VIEWER(cdf_file, list_file=list_file)
