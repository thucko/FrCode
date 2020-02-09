import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import askdirectory
from tkinter import *




def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd



file =file_dialog()

d1 = pd.read_csv(file, dtype='a')