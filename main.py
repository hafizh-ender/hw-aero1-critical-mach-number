from tkinter import *
from tkinter import ttk

import matplotlib.pyplot as plt
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure

from numpy import sqrt
from math import floor

import lib

# NOTE: To skip GUI implementation and directly observe the calculation, go to the buttonCommand function

####################################################################################################
#                                                                                                  #
#                                       GUI Implementation                                         #
#                                                                                                  #
####################################################################################################
root = Tk()
root.title("Critical Mach Number Calculator")

content = ttk.Frame(root, padding="3 3 12 12")
content.grid(column=0, row=0, sticky=(N, S, E, W), columnspan=2, rowspan=10)

gamma = StringVar()
ttk.Label(content, text=r"Specific Heat Ratio").grid(column=0, row=1, sticky=(W, E))

gamma_entry = ttk.Entry(content, width=7, textvariable=gamma)
gamma_entry.grid(column=1, row=1, sticky=(W, E))

Cp0 = StringVar()
ttk.Label(content, text="Incompressible Pressure Coefficient").grid(column=0, row=2, sticky=(W, E))

Cp0_entry = ttk.Entry(content, width=7, textvariable=Cp0)
Cp0_entry.grid(column=1, row=2, sticky=(W, E))

ttk.Label(content, text="\n").grid(column=0, row=4, columnspan=2)
ttk.Label(content, text="Critical Mach Number").grid(column=0, row=5, columnspan=2)

ttk.Label(content, text="Prandtl-Glauert Rule").grid(column=0, row=6, sticky=(W, E))
ttk.Label(content, text="Laitone Rule").grid(column=0, row=7, sticky=(W, E))
ttk.Label(content, text="Karman-Tsien Rule").grid(column=0, row=8, sticky=(W, E))

labelPrandtlGlauert_text = StringVar()
labelLaitone_text = StringVar()
labelKarmanTsien_text = StringVar()

labelPrandtlGlauert = ttk.Entry(content, textvariable=labelPrandtlGlauert_text, state="readonly")
labelPrandtlGlauert.grid(column=1, row=6, sticky=(W, E))

labelLaitone = ttk.Entry(content, textvariable=labelLaitone_text, state="readonly")
labelLaitone.grid(column=1, row=7, sticky=(W, E))

labelKarmanTsien = ttk.Entry(content, textvariable=labelKarmanTsien_text, state="readonly")
labelKarmanTsien.grid(column=1, row=8, sticky=(W, E))

ttk.Label(content, text="\n").grid(column=0, row=9, columnspan=2)
ttk.Label(content, text="Made by Hafizh Renanto Akhmad").grid(column=0, row=12, columnspan=2)

def buttonCommand(nIter=1000):
    ####################################################################################################
    #                                       Calculation                                                #
    ####################################################################################################
    Cp0 = float(Cp0_entry.get())
    gamma = float(gamma_entry.get())
    
    McrPrandtlGlauert = lib.bisectionMethodPrandtlGlauert(lib.minimum_float, lib.MmaxPrandtlGlauert, nIter, Cp0, gamma)
    McrLaitone = lib.bisectionMethodLaitone(lib.minimum_float, lib.MmaxLaitone(Cp0, gamma), nIter, Cp0, gamma)
    McrKarmanTsien = lib.bisectionMethodKarmanTsien(lib.minimum_float, lib.MmaxKarmanTsien(Cp0), nIter, Cp0, gamma)

    labelPrandtlGlauert_text.set(f"{McrPrandtlGlauert}")
    labelLaitone_text.set(f"{McrLaitone}")
    labelKarmanTsien_text.set(f"{McrKarmanTsien}")

    ####################################################################################################
    #                                       Plot                                                       #
    ####################################################################################################
    fig, ax = plt.subplots(1, 1)
    ax.set_xlabel(r"$M_\infty$")
    ax.set_ylabel(r"$C_p$")
    ax.grid()

    detail = 1000

    minCp = 0
    maxCp = lib.criticalPressureCoefficients(McrLaitone, gamma) * 1.5

    minM = 0
    maxM = 1

    MPrandtlGlauerts = [i/detail for i in range(minM, floor(lib.MmaxPrandtlGlauert*detail))]
    MLaitones = [i/detail for i in range(minM, floor(lib.MmaxLaitone(Cp0, gamma)*detail))]
    MKarmanTsiens = [i/detail for i in range(minM, floor(lib.MmaxKarmanTsien(Cp0)*detail))]
    MCrits = [i/detail for i in range(maxM*detail, minM, -1)]

    CpPrandtlGlauert = [lib.prandtlGlauert(M, Cp0) for M in MPrandtlGlauerts]
    CpLaitone = [lib.laitone(M, Cp0, gamma) for M in MLaitones]
    CpKarmanTsien = [lib.karmanTsien(M, Cp0) for M in MKarmanTsiens]
    CpCrs = [lib.criticalPressureCoefficients(M, gamma) for M in MCrits]

    ax.plot(MPrandtlGlauerts, CpPrandtlGlauert, label="Prandtl-Glauert Rule", color="tab:blue")
    ax.plot(MLaitones, CpLaitone, label="Laitone Rule", color="tab:orange")
    ax.plot(MKarmanTsiens, CpKarmanTsien, label="Karman-Tsien Rule", color="tab:green")
    ax.plot(MCrits, CpCrs, label="Critical Pressure Coefficient", color="black")

    resultMcrs = [McrPrandtlGlauert, McrLaitone, McrKarmanTsien]
    resultCpcrs = [lib.criticalPressureCoefficients(M, gamma) for M in resultMcrs]

    ax.scatter(resultMcrs[0], resultCpcrs[0], color="tab:blue")
    ax.scatter(resultMcrs[1], resultCpcrs[1], color="tab:orange")
    ax.scatter(resultMcrs[2], resultCpcrs[2], color="tab:green")

    for i in range(0, len(resultMcrs)):
        ax.annotate(f"({round(resultMcrs[i], 5)}, {round(resultCpcrs[i], 5)})", (resultMcrs[i], resultCpcrs[i]))

    ax.set_xlim(minM, maxM)
    ax.set_ylim(maxCp, minCp)
    ax.invert_yaxis()
    ax.legend()

    fig.tight_layout

    canvas = FigureCanvasTkAgg(fig, master=content)
    canvas.draw()
    canvas.get_tk_widget().grid(column=0, row=11, columnspan=2, sticky=(N, S, E, W))

    ttk.Label(content, text="Plot").grid(column=0, row=10, columnspan=2)

calculateButton = ttk.Button(content, text="Calculate", command=buttonCommand)
calculateButton.grid(column=0, row=3, columnspan=2)

root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

content.columnconfigure(1, weight=1)

root.mainloop()
root.quit()