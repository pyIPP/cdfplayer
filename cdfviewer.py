import sys, os, shutil
sys.path.append('/afs/ipp/home/g/git/python/repository')
import Tkinter as tk
import ttk
from tkFileDialog import askopenfilename, asksaveasfile
import tkMessageBox

import numpy as np
import ufiles, read_equ, mom2rz, tkhyper, rz2psi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
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

# Widget frame

        print('Setting GUI frame')
        self.viewframe = tk.Tk()
        self.viewframe.title('Profile viewer')
        self.viewframe.geometry('1600x960')

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
        savemenu.add_command(label="Save 1D u-file(t)"          , command=self.uf1t)
        savemenu.add_command(label="Save 1D u-file(rho)"        , command=self.uf1r)
        savemenu.add_command(label="Save 2D u-file(t, rho)"     , command=self.uf2)
        savemenu.add_command(label="Save 2D as .avi"            , command=self.smov)
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

        self.canv_frame = tk.Frame(self.viewframe, height=900)
        self.canv_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canv_frame.pack_propagate(0)

        toolframe = tk.Frame(self.viewframe, height=45)
        toolframe.pack(side=tk.BOTTOM, fill=tk.X)
        toolframe.pack_propagate(0)

# Plots

        self.nbpol = ttk.Notebook(self.canv_frame, name='nbpol')
        self.nbpol.pack(side=tk.TOP, fill=tk.X)
        self.nbpol.bind('<Button-1>', self.on_click)
  
        frame_eq = tk.Frame(self.nbpol)
        frame_1d = tk.Frame(self.nbpol)
        frame_2d = tk.Frame(self.nbpol)

        self.nbpol.add(frame_eq, text='Equilibirum')
        self.nbpol.add(frame_1d, text='Scalars    ')
        self.nbpol.add(frame_2d, text='Profiles   ')

        self.fig_eq = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_1d = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_2d = Figure(figsize=(14., 8.5), dpi=100)
        self.fig_eq.subplots_adjust(left=0.01, bottom=0.1, right=0.95, top=0.95)
        self.fig_1d.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.92, \
            wspace=0.3, hspace=0.3)
        self.fig_2d.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.92, \
            wspace=0.3, hspace=0.3)
        self.canv_eq = FigureCanvasTkAgg(self.fig_eq, master=frame_eq)
        self.canv_1d = FigureCanvasTkAgg(self.fig_1d, master=frame_1d)
        self.canv_2d = FigureCanvasTkAgg(self.fig_2d, master=frame_2d)
        self.canv_eq._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canv_1d._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canv_2d._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Navigation toolbars

        toolbar = NavigationToolbar2TkAgg( self.canv_eq, frame_eq)
        toolbar.update()
        toolbar = NavigationToolbar2TkAgg( self.canv_1d, frame_1d)
        toolbar.update()
        toolbar = NavigationToolbar2TkAgg( self.canv_2d, frame_2d)
        toolbar.update()

# Toolbar

        self.playfig  = tk.PhotoImage(file='%s/ButtonPlay.gif'  %locdir)
        self.pausefig = tk.PhotoImage(file='%s/ButtonPause.gif' %locdir)
        rewfig  = tk.PhotoImage(file='%s/ButtonBack.gif'     %locdir)
        prevfig = tk.PhotoImage(file='%s/ButtonPrevious.gif' %locdir)
        nextfig = tk.PhotoImage(file='%s/ButtonNext.gif'     %locdir)

        self.playbt = tk.Button(toolframe, command=self.play , image=self.playfig)
        rewbt  = tk.Button(toolframe, command=self.rewind, image=rewfig )
        prevbt = tk.Button(toolframe, command=self.prev  , image=prevfig)
        nextbt = tk.Button(toolframe, command=self.next  , image=nextfig)
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
        self.replot()


    def read(self, f_cdf):

        if not os.path.isfile(f_cdf):
            print('NetCDF file %s not found' %f_cdf)
            return

        print('Reading profiles from %s' %f_cdf)
        self.cv = netcdf.netcdf_file(f_cdf, 'r', mmap=False).variables
        for key, val in self.cv.iteritems():
            if key == 'rho':
                print key, val
                if key == 'Ts':
                    print key, val[0, :], val

        print 'Reading equilibrium from '+f_cdf
        rho = np.linspace(0.1, 1.0, 10)
        self.surf = read_equ.READ_EQU(f_cdf, rho=rho)


    def set_plots(self, pr_list=None):

        if os.path.isfile(self.list_file) and pr_list == None:
            f = open(self.list_file)
            pr_list = [s.strip() for s in f.readlines()]
            f.close()
            print('Using profile list from file %s' %self.list_file)
        elif pr_list == None:
            pr_list = ['CUR', 'TE', 'NE', 'V', 'BPOL', 'NEUTT', 'NIMP', 'PINJ', 'TAUEA']
 
        self.time = self.cv['TIME3'].data
        self.rho  = self.cv['X'][0, :]
        self.rhob = self.cv['XB'][0, :]
        self.ynp  = {}
        self.lbl  = {}
        self.dim  = {}
        self.unit = {}
        self.prof_list  = []
        self.trace_list = []
        for sig in pr_list:
            sig = sig.upper()
            if sig in self.cv.iterkeys():
                self.ynp[sig]  = self.cv[sig].data
                self.lbl[sig]  = self.cv[sig].long_name[:-2]
                self.dim[sig]  = self.cv[sig].dimensions
                self.unit[sig] = self.cv[sig].units
                self.nt = self.ynp[sig].shape[0]
                if len(self.cv[sig].data.shape) > 1:
                    self.prof_list.append(sig)
                else:
                    self.trace_list.append(sig)
            else:
                print('No variable %s in CDF file, ignoring!' %sig)

        if len(self.trace_list) > 15:
            print('Too many time traces, taking only the first 15'
            self.trace_list = self.trace_list[:15]
        if len(self.prof_list) > 15:
            print('Too many profiles, taking only the first 15'
            self.prof_list = self.prof_list[:15]
        nrow = 3
        ncol = 5

# Equilibrium

        npoints = 201
        Rin, Zin = mom2rz.mom2rz(self.surf.rc, \
            self.surf.rs, self.surf.zc, self.surf.zs, nthe=npoints)

        nt, self.n_rho, npoints = Rin.shape
        self.R = np.zeros((nt, self.n_rho, npoints+1))
        self.Z = np.zeros((nt, self.n_rho, npoints+1))
        self.R[:, :, :-1] = Rin
        self.Z[:, :, :-1] = Zin
        self.R[:, :, -1] = self.R[:, :, 0]
        self.Z[:, :, -1] = self.Z[:, :, 0]

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
            self.equline[jrho],  = ax_eq.plot(self.R[self.jt, jrho, :], self.Z[self.jt, jrho, :], 'g-')
        self.equline['axis'], = ax_eq.plot(self.surf.r_mag[self.jt], self.surf.z_mag[self.jt], 'g+')
        ax_eq.set_xlabel('R [cm]', fontsize=fsize)
        ax_eq.set_ylabel('z [cm]', fontsize=fsize)

# Plot AUG structures
        nshot = int(self.runid[0:5])
        try:
            import plot_aug
            plot_aug.vessel_pol(ax_eq, fac=100.)
        except:
            pass

        jplot = 1
        for trace in self.trace_list:
            ylab = '%s [%s]' %(trace, self.unit[trace].strip())

            ax_1d = self.fig_1d.add_subplot(nrow, ncol, jplot)
            ax_1d.set_xlabel('t [s]', fontsize=fsize)
            ax_1d.set_ylabel(ylab, fontsize=fsize)
            self.line1d[trace], = ax_1d.plot(self.time, self.ynp[trace][:], 'b-')
            ax_1d.grid('on')
            self.mark1d[trace], = ax_1d.plot(self.time[0], self.ynp[trace][0], 'go')
            jplot += 1
            ax_1d.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))
            ax_1d.yaxis.major.formatter._useMathText = True

        jplot = 1
        for prof in self.prof_list:
            if 'X' in self.dim[prof]:
                xrho = self.rho
                xlab = r'$\rho_{tor}$(X)'
            if 'XB' in self.dim[prof]:
                xrho = self.rhob
                xlab = r'$\rho_{tor}$(XB)'
            ylab = '%s [%s]' %(prof, self.unit[prof].strip())

            ax_2d = self.fig_2d.add_subplot(nrow, ncol, jplot)
            ax_2d.set_xlim([0, 1.])
            ax_2d.set_ylim([0, 1.1*np.max(self.ynp[prof])])
            ax_2d.set_xlabel(xlab, fontsize=fsize)
            ax_2d.set_ylabel(ylab, fontsize=fsize)
            self.line2d[prof], = ax_2d.plot(xrho, self.ynp[prof][0, :], 'b-')
            ax_2d.grid('on')
            ax_2d.ticklabel_format(axis='y', style='sci', scilimits=(-4,-4))
            ax_2d.yaxis.major.formatter._useMathText = True
            jplot += 1

        self.canv_eq.draw()
        self.canv_1d.draw()
        self.canv_2d.draw()
        self.jtab = 0

        self.replot()


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
        rz = rz2psi.RZ2PSI(self.cdf_file, it=self.jt)

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
        eqd = eqdsk.EQDSK()
        eqd.write_eqdsk(eqd_d, f_eqdsk, self.time[self.jt], f_cdf=self.cdf_file)


    def uf1t(self):

        self.ff = False
        self.rw = False

        sshot = self.runid[0:5]
        tlbl = 'Time:'.ljust(20) + 's'
        for trace in self.trace_list:
            dlbl = trace.ljust(20) + self.unit[trace].strip()
            uf_d = { 'pre':'A', 'ext':trace, 'shot':sshot, \
                     'grid': {'X': {'lbl':tlbl, 'arr':self.time} }, \
                     'data': {'lbl':dlbl, 'arr':self.ynp[trace][:]} }
            ufiles.WU(uf_d)


    def uf1r(self):

        self.ff = False
        self.rw = False

        sshot = self.runid[0:5]
        xlbl = 'rho_tor'
        for prof in self.prof_list:
            dlbl = prof.ljust(20) + self.unit[prof].strip()
            uf_d = { 'pre':'A', 'ext':prof, 'shot':sshot, \
                     'scal': [['Time:'.ljust(20) + 's', self.time[self.jt]]], 
                     'grid': {'X': {'lbl':xlbl, 'arr':self.rho} }, \
                     'data': {'lbl':dlbl, 'arr':self.ynp[prof][self.jt, :]} }
            ufiles.WU(uf_d)


    def uf2(self):

        self.ff = False
        self.rw = False
        
        sshot = self.runid[0:5]
        tlbl = 'Time'.ljust(20) + 'Seconds'
        xlbl = 'rho_tor'
        for prof in self.prof_list:
            dlbl = prof.ljust(20) + self.unit[prof].rstrip()
            uf_d = { 'pre':'A', 'ext':prof + '2', 'shot':sshot, \
                     'grid': {'X':{'lbl':tlbl, 'arr':self.time}, \
                              'Y':{'lbl':xlbl, 'arr':self.rho} }, \
                     'data': {'lbl':dlbl, 'arr':self.ynp[prof]} }
            ufiles.WU(uf_d)


    def cdf2tre(self):

        homdir = os.getenv('HOME')
        os.system('mkdir -p %s/shotfiles/TRE' %homdir)
        fsfh = '%s/tr_client/AUGD/TRE00000.sfh' %homdir
        source = '/afs/ipp/home/t/transp/pub/TRE00000.sfh.temp'
        if not os.path.isfile(self.cdf_file):
            print('%s not found' %self.cdf_file)
            return
        try:
            shutil.copy2(source, fsfh)
        except IOError, e:
            print('Unable to copy file %s\n%s' %(source, e))

        import cdf2sf
        cdf2sf.CDF2TRE(self.cdf_file, sfh_tre=fsfh, vplot=True)


    def tr2sf(self):

        import cdf2sf
        cdf2sf.CDF2TRA(self.runid)


    def smov(self):

        self.ff = False
        self.rw = False

# Store png's
        order = int(np.log10(self.nt)) + 1 # To have nice alphabetic sorted files
        fmt = '%' + str(order) + 'd'
        for jtim in range(self.nt):
            self.jt = jtim
            self.replot()
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


    def replot(self, reset=True):

#        if hasattr(self, 'crntsc') and reset:
        if reset:
            self.crntsc.set(float(self.jt)/(self.nt - 1)*100.)

        strtim = 't =%6.3f s' %self.time[self.jt]
        if self.jtab == 0:
            for jrho in range(self.n_rho):
                self.equline [jrho].set_data(self.R[self.jt, jrho, :], self.Z[self.jt, jrho, :])
            self.txt_eq.set_text(strtim)
# Plot mag. axis:
            self.equline['axis'].set_data(self.surf.r_mag[self.jt], self.surf.z_mag[self.jt])
            self.canv_eq.draw()
            
        elif self.jtab == 1:
            for trace in self.trace_list:
                self.mark1d[trace].set_data(self.time[self.jt], self.ynp[trace][self.jt])
            self.txt_1d.set_text(strtim)
            self.canv_1d.draw()

        elif self.jtab == 2:
            for prof in self.prof_list:
                self.line2d[prof].set_ydata(self.ynp[prof][self.jt, :])
            self.txt_2d.set_text(strtim)
            self.canv_2d.draw()


    def prev(self):

        self.ff = False
        self.rw = False
        if self.jt > 0:
            self.jt -= 1
            self.inNextOrPrev = True
            self.replot()
            self.inNextOrPrev = False


    def next(self):

        self.ff = False
        self.rw = False
        if self.jt < self.nt - 1:
            self.jt += 1
            self.inNextOrPrev = True
            self.replot()
            self.inNextOrPrev = False


    def jump(self, arg):

        if self.inNextOrPrev:
            return
        self.jt = int(float(arg)/100. * (self.nt-1))
        self.replot(reset=False)


    def pause(self):

        self.ff = False
        self.rw = False
        self.playbt['image']   = self.playfig
        self.playbt['command'] = self.play
        self.replot()


    def play(self):

        self.ff = True
        self.rw = False
        self.playbt['image']   = self.pausefig
        self.playbt['command'] = self.pause
        while self.ff and self.jt < self.nt - 1:
            self.jt += 1
            self.replot()
            self.canv_eq.start_event_loop(self.timeout)


    def rewind(self):

        self.ff = False
        self.rw = True
        while self.rw and self.jt > 0:
            self.jt -= 1
            self.replot()
            self.canv_eq.start_event_loop(self.timeout)


    def callcdf(self):

        dir_in = '%s/tr_client/%s' %(os.getenv('HOME'), self.tok)
        self.cdf_file = askopenfilename(initialdir=dir_in, filetypes=[("All formats", "*.CDF")])
        self.load()


    def callset(self):

        self.list_file = askopenfilename(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        self.jt = 0
        self.equline  = {}
        self.mark1d = {}
        self.line1d = {}
        self.line2d = {}
        self.set_plots()


    def callsave(self):

        f_out = asksaveasfile(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        out_txt = ''
        for prof in self.prof_list:
            out_txt += prof + '\n'
        for sig in self.trace_list:
            out_txt += sig + '\n'
        f_out.write(out_txt)
        f_out.close()


    def editlist(self, profiles=None):


        class ProfilesDialog(object):

            def __init__(self, parent, viewer, profiles=None):
                self.viewer = viewer
                sig_list = self.viewer.prof_list + self.viewer.trace_list
                top = self.top = tk.Toplevel(parent)
                top.resizable(0, 0)
                tk.Label(top, text="Profiles to be viewed").grid(row=0, columnspan=3)
                self.t = tk.Text(top, width=40, height=10)
                self.t.grid(row=1, columnspan=2)
                self.t.insert('end', "\n".join(sig_list if profiles==None else profiles))
                scr = tk.Scrollbar(top)
                scr.config(command=self.t.yview)
                scr.grid(row=1, column=2, sticky=tk.N+tk.S)
                self.t.config(yscrollcommand=scr.set)
                b = tk.Button(top, text="OK", command=self.ok)
                b.grid(row=2, column=0)
                b = tk.Button(top, text="Cancel", command=self.top.destroy)
                b.grid(row=2, column=1)

            def ok(self):
                newsigs = self.t.get("0.0", "end-1c").split()
                oldprofiles = self.viewer.prof_list
                oldtraces   = self.viewer.trace_list
                try:
                    self.viewer.jt = 0
                    self.viewer.equline = {}
                    self.viewer.line1d = {}
                    self.viewer.line2d = {}
                    self.viewer.set_plots(pr_list=newsigs)
                    self.top.destroy()
                except Exception, e:
                    tkMessageBox.showerror('Error', 'This went wrong: %s'%repr(e))
                    self.viewer.prof_list  = oldprofiles
                    self.viewer.trace_list = oldtraces

        dialog = ProfilesDialog(self.viewframe, self, profiles=profiles)


    def load(self):

        self.jt = 0
        self.runid = self.cdf_file[-12:-4]
        self.equline  = {}
        self.mark1d = {}
        self.line1d = {}
        self.line2d = {}
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
            lbl = tk.Label(dispframe, text=key, fg="#000000")
            var = tk.Entry(dispframe, width=12)
            lbl.grid(row=nrow, column=0, sticky=tk.W+tk.E)
            var.grid(row=nrow, column=1, sticky=tk.W+tk.E)
            var.insert(0, '')
            self.disp_d[key] = var
            nrow += 1
        dim = tk.IntVar()
        dim.set(3)
        rb1 = tk.Radiobutton(dispframe, text='Scalar' , variable=dim, val=1)
        rb2 = tk.Radiobutton(dispframe, text='Profile', variable=dim, val=2)
        rb3 = tk.Radiobutton(dispframe, text='All', variable=dim, val=3)
        rb1.grid(row=nrow, column=0, sticky=tk.W+tk.E)
        rb2.grid(row=nrow, column=1, sticky=tk.W+tk.E)
        rb3.grid(row=nrow, column=2, sticky=tk.W+tk.E)
        self.disp_d['dim'] = dim
        nrow += 1
        but1 = tk.Button(dispframe, text='Print variable', command=self.print_vars)
        but2 = tk.Button(dispframe, text='Close', command=dispframe.destroy)
        but1.grid(row=nrow, column=0)
        but2.grid(row=nrow, column=1)

        dispframe.mainloop()


    def print_vars(self):

        skey  = self.disp_d['Var_name'].get().upper().strip()
        sunit = self.disp_d['Unit'].get().upper().strip()
        sdesc = self.disp_d['Description'].get().upper().strip()
        dim = self.disp_d['dim'].get()
        for key, val in self.cv.iteritems():
            unit  = val.units.strip()
            descr = val.long_name.strip()
            key_flag  = (skey  == '') or (skey  in key)
            unit_flag = (sunit == '') or (sunit in unit)
            desc_flag = (sdesc == '') or (sdesc in descr)
            if key_flag and unit_flag and desc_flag:
                if len(val.data.shape) == dim or dim == 3:
                    print key.ljust(10), unit.ljust(16), descr.ljust(20)
        print('')


if __name__ == "__main__":


    cdf_file  = '/afs/ipp/home/g/git/tr_client/AUGD/23076/W01/23076W01.CDF'
    list_file = ''

    if len(sys.argv) > 2:
        f = sys.argv[2]
        print f
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
