import Tkinter as tk

bgcol = '#aabbaa'
fgcol = '#000000'
lbcol = '#dddddd'
encol = '#fcfcfc'
rbcol = '#cceecc'
ckcol = '#cceecc'
mncol = '#cceecc'
nxcol = '#50c050'
excol = '#90a090'
qtcol = '#f00000'
tbcol = '#ececec'

lbwid = 10
enwid = 7
btwid = 8
rbwid = 8
ckwid = 9
mnwid = 8
icon_dir = '/afs/ipp/home/g/git/python/qticons/'
xpad = 2
ypad = 5

padx = 2
pady = 5


class TKSTY:

    nrowf = 1


    def myLabel(self, frame, txt, nrow=0, ncol=0, pos=tk.W, wid=lbwid, \
                bgc=lbcol, fgc=fgcol, xpad=xpad, ypad=ypad, \
                justify=tk.CENTER, reli=tk.GROOVE):

        lbl = tk.Label(frame, width=wid, text=txt, bg=bgc, fg=fgc, bd=2, \
            relief=reli, justify=justify)
        lbl.grid(row=nrow, column=ncol, sticky=pos, padx=xpad, pady=ypad)
        return lbl


    def myButton(self, frame, txt, comm, bgc=nxcol, nrow=0, ncol=0, \
                 pos=tk.W, xpad=xpad, ypad=ypad, img=None, justify=tk.CENTER):

        but = tk.Button(frame, text=txt, bg=bgc, command=comm, image=img, \
            justify=justify)
        but.grid(row=nrow, column=ncol, sticky=pos, padx=xpad, pady=ypad)
        return but


    def myEntry(self, frame, txt, wid=enwid, bgc=encol, nrow=0, ncol=0, \
                pos=tk.W, xpad=xpad, ypad=ypad, justify=tk.CENTER):

        ent = tk.Entry(frame, width=wid, bg=bgc, justify=justify)
        ent.grid(row=nrow, column=ncol, sticky=pos, padx=xpad, pady=ypad)
        ent.insert(0, txt)
        return ent


    def myRb(self, frame, var, val, wid=rbwid, bgc=rbcol, txt='', nrow=0, \
             ncol=0, pos=tk.W, xpad=xpad, ypad=ypad, justify=tk.CENTER):

        rb = tk.Radiobutton(frame, width=wid, text=txt, variable=var, \
            value=val, bg=bgc, bd=2, relief=tk.GROOVE, justify=justify)
        rb.grid(row=nrow, column=ncol, sticky=tk.W, padx=xpad, pady=ypad)
        return rb


    def myCb(self, frame, var, onval=True, offval=False, wid=ckwid, bgc=ckcol, \
             txt='', nrow=0, ncol=0, pos=tk.W, xpad=xpad, ypad=ypad, justify=tk.CENTER):

        cb = tk.Checkbutton(frame, width=wid, text=txt, variable=var, onvalue=onval, \
            offvalue=offval, bg=ckcol, bd=2, relief=tk.GROOVE, justify=justify)
        cb.grid(row=nrow, column=ncol, sticky=tk.W, padx=5, pady=2)
        return cb


    def myMb(self, frame, txt, wid=mnwid, bgc=mncol, nrow=0, ncol=0, \
             pos=tk.E+tk.W, xpad=xpad, ypad=ypad):

        mb = tk.Menubutton(frame, width=wid, text=txt, bg=bgc, bd=2, relief=tk.GROOVE)
        mb.grid(row=nrow, column=ncol, sticky=pos, padx=xpad, pady=ypad)
        return mb


    def myToolbar(self, frame, bgc=tbcol):

        toolbar = tk.Frame(frame, bg=bgc, bd=2, relief=tk.GROOVE)
        toolbar.grid(row=0, column=0, sticky=tk.W+tk.E)
        return toolbar


    def mySubframe(self, frame, nrow=-1, ncol=0, bgc=bgcol):

        if nrow == -1:
            nrow = self.nrowf
        subframe = tk.Frame(frame, bg=bgc)
        subframe.grid(row=nrow, column=ncol, sticky=tk.W+tk.E)

        self.nrowf += 1
        return subframe


    def myCanvas(self, frame, nrow=-1, ncol=0, bgc=bgcol, height=1):

        if nrow == -1:
            nrow = self.nrowf
        canvas=tk.Canvas(frame, bg=bgc, height=height)
        canvas.grid(row=nrow, column=ncol, sticky=tk.W+tk.E)

        self.nrowf += 1
        return canvas
