import os, re
import webbrowser
try:
    import Tkinter as tk
except:
    import tkinter as tk

os.environ['BROWSER'] = '/usr/bin/firefox'

 
class HyperlinkMessageBox():


    hyperlinkPattern = re.compile(r'<a href="(?P<address>.*?)">(?P<title>.*?)'

                                  '</a>')
    def __init__(self, title=None, message=None, geom="300x80", **options):

        master = tk.Tk()
        master.geometry(geom)
        self.text = tk.Text(master, wrap=tk.WORD, bg=master.cget('bg'),
                         height=master.cget('height'))
        self._formatHyperLink(message)
        self.text.config(state=tk.DISABLED)
        self.text.pack(side=tk.TOP, fill=tk.X)
        tk.Button(master, text="Ok",
               command=lambda *a, **k: master.destroy()).pack()
        master.mainloop()

 
    def _formatHyperLink(self, message):

        text = self.text
        start = 0
        for index, match in enumerate(self.hyperlinkPattern.finditer(message)):
            groups = match.groupdict()
            text.insert("end", message[start: match.start()])
            text.insert("end", groups['title'])
            text.tag_add(str(index),
                         "end-%dc" % (len(groups['title']) + 1),
                         "end",)
            text.tag_config(str(index),
                            foreground="blue",
                            underline=1)
            text.tag_bind(str(index),
                          "<Enter>",
                          lambda *a, **k: text.config(cursor="arrow"))
            text.tag_bind(str(index),
                          "<Leave>",
                          lambda *a, **k: text.config(cursor="arrow"))
            text.tag_bind(str(index),
                          "<Button-1>",
                          self._callbackFactory(groups['address']))
            start = match.end()
        else:
            text.insert("end", message[start:])
 

    def _callbackFactory(self, url):

        return lambda *args, **kwargs: webbrowser.open(url)


if __name__ == "__main__":
 
    HyperlinkMessageBox(title="My App", message='Some message <a href="http://www.google.com">Google</a>.')
