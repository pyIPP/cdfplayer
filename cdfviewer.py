import sys, os, shutil
import Tkinter as tk
from tkFileDialog import askopenfilename, asksaveasfile
import tkMessageBox

import numpy as np
import ufiles, read_equ, mom2rz, tkhyper, rz2psi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.pylab as plt
from netCDF4 import Dataset


nxcol = '#50c050'
qtcol = '#f00000'
svcol = '#a0a0d0'
encol = '#dddddd'


class VIEWER:


    timeout = 1e-10
    dir_set = '%s/python/cdfp' %os.getenv('HOME')

    inNextOrPrev = False


    def __init__(self, f_cdf_in, list_file='', tok='AUGD'):

        self.tok = tok

        self.runid     = f_cdf_in[-12:-4]
        self.cdf_file  = f_cdf_in
        self.list_file = list_file

# Widget frame

        print('Setting GUI frame')
        self.viewframe = tk.Tk()
        self.viewframe.title('Profile viewer')

# Menu bar

        menubar = tk.Menu(self.viewframe)

        filemenu = tk.Menu(menubar, tearoff=0)
        filemenu.add_command(label="Play"       , command=self.play)
        filemenu.add_command(label="Load run...", command=self.callcdf)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=sys.exit)

        editmenu = tk.Menu(menubar, tearoff=0)
        editmenu.add_command(label="Load setup...", command=self.callback)
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

        menubar.add_cascade(label="File"  , menu=filemenu)
        menubar.add_cascade(label="Setup" , menu=editmenu)
        menubar.add_cascade(label="Output", menu=savemenu)
        menubar.add_cascade(label="Help"  , menu=helpmenu)

        self.viewframe.config(menu = menubar)

# Figure frame

        self.outframe = tk.Frame(self.viewframe)
        self.viewfig = Figure()
        self.canvas = FigureCanvasTkAgg(self.viewfig, master=self.outframe)

# Navigation toolbar
        toolbar = NavigationToolbar2TkAgg( self.canvas, self.viewframe)
        toolbar.update()

# Toolbar

        toolframe = tk.Frame(self.viewframe, height=45)
        toolframe.pack_propagate(0)
        toolframe.pack(side=tk.BOTTOM, fill=tk.X)

        self.playfig  = tk.PhotoImage(file='ButtonPlay.gif')
        self.pausefig = tk.PhotoImage(file='ButtonPause.gif')
        rewfig  = tk.PhotoImage(file='ButtonRewind.gif')
        begfig  = tk.PhotoImage(file='ButtonFirst.gif')
        endfig  = tk.PhotoImage(file='ButtonLast.gif')
        prevfig = tk.PhotoImage(file='ButtonPrevious.gif')
        nextfig = tk.PhotoImage(file='ButtonNext.gif')

        self.playbt = tk.Button(toolframe, command=self.play , image=self.playfig)
        rewbt  = tk.Button(toolframe, command=self.rewind, image=rewfig )
        begbt  = tk.Button(toolframe, command=self.begin , image=begfig )
        endbt  = tk.Button(toolframe, command=self.end   , image=endfig )
        prevbt = tk.Button(toolframe, command=self.prev  , image=prevfig)
        nextbt = tk.Button(toolframe, command=self.next  , image=nextfig)
        for but in self.playbt, rewbt, begbt, endbt, prevbt, nextbt:
            but.pack(side=tk.LEFT)

        self.crntsc = tk.Scale(toolframe, command=self.jump, \
                               orient=tk.HORIZONTAL, length=300)
        self.crntsc.pack(side=tk.LEFT)

        tk.Button(toolframe, text='Equ', bg=nxcol, command=self.view_equ).pack(side=tk.RIGHT)
        tk.Button(toolframe, text='2D', bg=nxcol, command=self.view2d).pack(side=tk.RIGHT)
        tk.Button(toolframe, text='1D', bg=nxcol, command=self.view1d).pack(side=tk.RIGHT)

        self.switch = 0

# Set up plot
        self.load()
#        try:
#            self.load()
#        except:
#            print('Select a CDF file via File->Load run...')

        self.viewframe.mainloop()


    def about(self):

        mytext = 'Documentation at <a href="http://www.aug.ipp.mpg.de/aug/manuals/transp/output/cdfplayer.html">cdfplayer homepage</a>'
        h = tkhyper.HyperlinkMessageBox("Help", mytext, "340x60")


    def read(self, f_cdf):

        if not os.path.isfile(f_cdf):
            print('CDF file %s not found' %f_cdf)
            return

        print('Reading profiles from %s' %f_cdf)
        self.cv = Dataset(f_cdf, 'r', format='NETCDF4').variables
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
            pr_list = ['CUR', 'TE', 'NE', 'V', 'BPOL', 'NEUTT']
 
        self.xnp  = self.cv['X'][:]
        self.xbnp = self.cv['XB'][:]
        self.time = self.cv['TIME3'][:]
        self.rho  = self.xnp[0, :]
        self.rhob = self.xbnp[0, :]
        self.ynp  = {}
        self.lbl  = {}
        self.dim  = {}
        self.unit = {}
        self.prof_list  = []
        self.trace_list = []
        for sig in pr_list:
            sig = sig.upper()
            if sig in self.cv.iterkeys():
                self.ynp[sig]  = self.cv[sig][:]
                self.lbl[sig]  = self.cv[sig].long_name[:-2]
                self.dim[sig]  = self.cv[sig].dimensions
                self.unit[sig] = self.cv[sig].units
                self.nt = self.ynp[sig].shape[0]
                if len(self.cv[sig].shape) > 1:
                    self.prof_list.append(sig)
                else:
                    self.trace_list.append(sig)
            else:
                print('%s not in CDF variables, dropping!' %sig)

        nprof = max(len(self.prof_list), len(self.trace_list))
        self.nrow = 3
        cols = float(nprof)/float(self.nrow)
        self.ncol = int(cols) + 1
        if cols == int(cols):
            self.ncol += -1

        self.R, self.Z = mom2rz.mom2rz(self.surf.rc, \
            self.surf.rs, self.surf.zc, self.surf.zs, nthe=201)

        self.n_rho = self.R.shape[1]


    def view1d(self):

        self.switch = 1
        self.set_can()


    def view2d(self):

        self.switch = 0
        self.set_can()


    def view_equ(self):

        self.switch = 2
        self.set_can()


    def set_can(self):

        self.viewfig.clf()
        self.viewfig.set_canvas(self.canvas)
        self.txt = self.viewfig.text(0.5, 0.98, '', ha='center', va='top')
        fsize = 12

        if self.switch == 1:
            jplot = 1
            for trace in self.trace_list:
                unit = ' [%s]' %self.unit[trace].strip()
                if unit.strip() == '[]':
                    unit = ''
                ylab = trace + unit

                self.axes1d = self.viewfig.add_subplot(self.nrow, self.ncol, jplot)
                self.axes1d.set_xlabel('t [s]', fontsize=fsize)
                self.axes1d.set_ylabel(ylab, fontsize=fsize)
                self.line1d[trace], = self.axes1d.plot(self.time, self.ynp[trace][:], 'b-')
                self.axes1d.grid('on')
                self.mark1d[trace], = self.axes1d.plot(self.time[0], self.ynp[trace][0], 'go')
                jplot += 1

        elif self.switch == 0:
            jplot = 1
            for prof in self.prof_list:
                if 'X' in self.dim[prof]:
                    xrho = self.rho
                    xlab = r'$\rho_{tor}$(X)'
                if 'XB' in self.dim[prof]:
                    xrho = self.rhob
                    xlab = r'$\rho_{tor}$(XB)'
                unit = ' [%s]' %self.unit[prof].strip()
                if unit.strip() == '[]':
                    unit = ''
                ylab = prof + unit

                self.axes2d = self.viewfig.add_subplot(self.nrow, self.ncol, jplot)
                self.axes2d.set_xlim([0, 1.])
                self.axes2d.set_ylim([0, 1.1*np.max(self.ynp[prof])])
                self.axes2d.set_xlabel(xlab, fontsize=fsize)
                self.axes2d.set_ylabel(ylab, fontsize=fsize)
                self.line2d[prof], = self.axes2d.plot(xrho, self.ynp[prof][0, :], 'b-')
                self.axes2d.grid('on')
                jplot += 1


        else:

            self.equ = self.viewfig.add_subplot(1, 1, 0, aspect='equal')
            for jrho in range(self.n_rho):
                self.equline[jrho],  = self.equ.plot(self.R[self.jt, jrho, :], self.Z[self.jt, jrho, :], 'g-')
                self.equline2[jrho], = self.equ.plot(self.R[self.jt, jrho, [-1, 0]], self.Z[self.jt, jrho, [-1, 0]], 'g-')
            self.equline['axis'], = self.equ.plot(self.surf.r_mag[self.jt], self.surf.z_mag[self.jt], 'g+')

# Plot AUG structures
            nshot = int(self.runid[0:5])
            if self.tok == 'AUGD':
                import kk_20140416 as kk
                kgc = kk.kkGCd0(nshot)
                for key in kgc.gc_x.iterkeys():
                    self.equ.plot(100*kgc.gc_x[key], 100*kgc.gc_y[key], 'b-')

        self.replot()


    def eqdisk(self):

        import eqdsk, read_eq_sf

        cdf = Dataset(self.cdf_file, 'r', format='NETCDF4')
        cdf1d = ['PCUR', 'BZ', 'BZXR', 'RAXIS', 'YAXIS']
        cdf2d = ['PLFLX', 'Q', 'GFUN', 'PMHD_IN']
        coord = ['XB', 'X']

#==================
# Psi(R, z), j(R, z)
#==================

        cdf_d = {}
        cdf_d['n_the'] = '101'
        cdf_d['time']  = self.time[self.jt]

        eq = read_eq_sf.get_grid(self.runid)
        Rgrid = eq['Ri'][:, 0]
        zgrid = eq['Zj'][:, 0]
        rz = rz2psi.RZ2PSI(self.cdf_file, it=self.jt)

        data = {}
        for lbl in cdf1d + cdf2d + coord:
            data[lbl] = cv[lbl][self.jt]

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


    def replot(self):

        if hasattr(self, 'crntsc'):
            self.crntsc.set(float(self.jt)/(self.nt - 1)*100.)
        
        if self.switch == 0:
            #profile plot:
            strtim = '%s    Time:%9.3f s' %(self.runid, self.time[self.jt])
            for prof in self.prof_list:
                self.line2d[prof].set_ydata(self.ynp[prof][self.jt, :])
        elif self.switch == 1:
            strtim = '%s    Time:%9.3f s' %(self.runid, self.time[self.jt])
            for trace in self.trace_list:
                self.mark1d[trace].set_data(self.time[self.jt], self.ynp[trace][self.jt])

        else:
            #Equilibrium plot:
            strtim = '%s    Time:%9.3f s' %(self.runid, self.surf.time[self.jt])
            for jrho in range(self.n_rho):
                #plot flux surfaces:
                self.equline [jrho].set_data(self.R[self.jt, jrho, :], self.Z[self.jt, jrho, :])
                self.equline2[jrho].set_data(self.R[self.jt, jrho, [-1, 0]], self.Z[self.jt, jrho, [-1, 0]])
            
            #plot mag. axis:
            self.equline['axis'].set_data(self.surf.r_mag[self.jt], self.surf.z_mag[self.jt])

        self.txt.set_text(strtim)
        self.canvas.draw()


    def begin(self):

        self.ff = False
        self.rw = False
        self.jt = 0
        self.replot()


    def end(self):

        self.ff = False
        self.rw = False
        self.jt = self.nt-1
        self.replot()


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
        if self.jt < self.nt-1:
            self.jt += 1
            self.inNextOrPrev = True
            self.replot()
            self.inNextOrPrev = False


    def jump(self, arg):

        if self.inNextOrPrev:
            return
        self.jt = int(float(arg)/100. * (self.nt-1))
        self.replot()


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
        while self.ff and self.jt < self.nt-1:
            self.jt += 1
            self.replot()
            self.canvas.start_event_loop(self.timeout)


    def rewind(self):

        self.ff = False
        self.rw = True
        while self.rw and self.jt > 0:
            self.jt -= 1
            self.replot()
            self.canvas.start_event_loop(self.timeout)


    def callcdf(self):

        dir_in = os.getenv('HOME') + '/tr_client/AUGD'
        self.cdf_file = askopenfilename(initialdir=dir_in, filetypes=[("All formats", "*.CDF")])
        self.load()


    def callback(self):

        self.list_file = askopenfilename(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        self.jt = 0
        self.equline  = {}
        self.equline2 = {}
        self.mark1d = {}
        self.line1d = {}
        self.line2d = {}
        self.set_plots()
        self.resize_frame()
        self.set_can()


    def callsave(self):

        f_out = asksaveasfile(initialdir=self.dir_set, filetypes=[("All formats", "*.txt")])
        out_txt = ''
        for prof in self.prof_list:
            out_txt += prof + '\n'
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
                    self.viewer.equline2 = {}
                    self.viewer.line1d = {}
                    self.viewer.line2d = {}
                    self.viewer.set_plots(pr_list=newsigs)
                    self.viewer.resize_frame()
                    self.viewer.set_can()
                    self.top.destroy()
                except Exception, e:
                    tkMessageBox.showerror('Error', 'This went wrong: %s'%repr(e))
                    self.viewer.prof_list  = oldprofiles
                    self.viewer.trace_list = oldtraces

        dialog = ProfilesDialog(self.viewframe, self, profiles=profiles)


    def resize_frame(self):

        geo = str(max(self.ncol*350, 700))+'x950'
        print('Frame geometry %s' %geo)
        self.viewframe.geometry(geo)
        self.outframe.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.viewfig.set_size_inches((4.5*self.ncol, 10.8))
        self.viewfig.subplots_adjust(left=0.08, bottom=0.05, right=0.95, top=0.95)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)


    def load(self):

        self.jt = 0
        self.runid = self.cdf_file[-12:-4]
        self.equline  = {}
        self.equline2 = {}
        self.mark1d = {}
        self.line1d = {}
        self.line2d = {}
        self.read(self.cdf_file)
        self.set_plots()
        self.resize_frame()
        self.set_can()


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
            sshot = f[:5]
            tail  = f[5:]
            f = '%s/tr_client/AUGD/%s/%s/%s.CDF' %(os.getenv('HOME'), sshot, tail, f)
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
