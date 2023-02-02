import pickle as pl
from tkinter import filedialog
import matplotlib.pyplot as plt
import tkinter as tk

import matplotlib

matplotlib.interactive(True)

file = filedialog.askopenfile(mode='r',
                              filetypes=[('Pickle File', '*.pickle')])

figx = pl.load(open(file.name, 'rb'))
figx.show()
tk.mainloop()
print("Yuh")