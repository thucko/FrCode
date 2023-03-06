"""
Author: Tim Hucko
Version 0
"""


from PyPDF2 import PdfFileMerger
from tkinter.filedialog import askopenfilenames
from tkinter.filedialog import asksaveasfile
from tkinter import simpledialog
from tkinter import *


def save_file():
    f = asksaveasfile(mode='w', defaultextension=".pdf")
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f

Tk().withdraw()

files = askopenfilenames(title='Choose Files')
merger = PdfFileMerger()
for pdf in files:
    merger.append(pdf)
root = Tk()
root.withdraw()
path = simpledialog.askstring(title='Save File', prompt='Enter File Name:')
path = path+'.pdf'

merger.write(path)
merger.close()

