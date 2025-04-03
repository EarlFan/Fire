#!/usr/bin/python3.12

import matplotlib.font_manager
import numpy as np
# import cantera as ct
import pandas as pd
import math

import matplotlib.pyplot as plt
import matplotlib

T1 = 747.672441027809
p1 = 35.594e3
rho1 = 0.18075
u1 = 487.34

x_min = 0.0
x_max = 0.1

headerList2 = ["x", "rhoMix", "p", "T", "rho_XE", "rho_H2", "rho_H", "rho_OH",
               "rho_O", "rho_H2O", "rho_HO2", "rho_H2O2", "rho_O2", "rho_N2", "pMax",
               "u", "v", "w", "gradRho_norm", "tau_c", "tau_r", "lvl", "Xi"]

MW_list = [39.948, 2.01588, 1.00794, 17.00734,
           15.9994, 18.01528, 33.00677, 34.0147, 31.9988]

startFrame = 0

endFrame = 45

frameGap = 5

frameDt = 0.1

df_list = [pd.DataFrame() for x in range(int(endFrame/frameGap))]
df_list_raw = [pd.DataFrame() for x in range(int(endFrame/frameGap))]

parentDir = ""

for i in range(startFrame, int(endFrame/frameGap), 1):

    fileName = parentDir+"cuts/cut1D_CPU0_TIME"+str(i*frameGap)+".out"

    df_list_raw[i] = pd.read_csv(
        fileName, skip_blank_lines=True, sep=r'\s+', names=headerList2)

    df_list_raw[i]['Y_OH'] = df_list_raw[i]['rho_OH']/df_list_raw[i]['rhoMix']

    df_list_raw[i]['Y_H2O'] = df_list_raw[i]['rho_H2O']/df_list_raw[i]['rhoMix']

    df_list[i] = df_list_raw[i].sort_values(by=["x"])


pparam = dict(xlabel='x', ylabel='T')

m_markersize = 3
m_markevery = 0.1

# temperature history
if(1):
    fig, ax = plt.subplots()

    for i in range(startFrame, int(endFrame/frameGap), 1):
        labelList = f"t = {i*frameGap*frameDt:.1f}ms"
        ax.plot(df_list[i]["x"], df_list[i]["T"], label=labelList,
                markersize=m_markersize, marker="s", markevery=0.08)

    ax.set_xlim([x_min, x_max])

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylabel('T(K)')
    ax.set_xlabel('x(m)')

    fig.suptitle(r'Temperature distribution')
    fig.set_size_inches(6, 4)
    fig.tight_layout()
    fig.savefig("T.jpg", dpi=300)

# Y_OH history
if(1):
    fig, ax = plt.subplots()

    for i in range(startFrame, int(endFrame/frameGap), 1):
        labelList = f"t = {i*frameGap*frameDt:.1f}ms"
        ax.plot(df_list[i]["x"], df_list[i]["Y_OH"], label=labelList,
                markersize=m_markersize, marker="s", markevery=0.08)

    ax.set_xlim([x_min, x_max])

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylabel('Y_OH')
    ax.set_xlabel('x(m)')

    fig.suptitle(r'OH mass fraction distribution')
    fig.set_size_inches(6, 4)
    fig.tight_layout()
    fig.savefig("Y_OH.jpg", dpi=300)

# Y_H2O history
if(1):
    fig, ax = plt.subplots()

    for i in range(startFrame, int(endFrame/frameGap), 1):
        labelList = f"t = {i*frameGap*frameDt:.1f}ms"
        ax.plot(df_list[i]["x"], df_list[i]["Y_H2O"], label=labelList,
                markersize=m_markersize, marker="s", markevery=0.08)

    ax.set_xlim([x_min, x_max])

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylabel('Y_H2O')
    ax.set_xlabel('x(m)')

    fig.suptitle(r'H2O mass fraction distribution')
    fig.set_size_inches(6, 4)
    fig.tight_layout()
    fig.savefig("Y_H2O.jpg", dpi=300)

# u history
if(1):
    fig, ax = plt.subplots()

    for i in range(startFrame, int(endFrame/frameGap), 1):
        labelList = f"t = {i*frameGap*frameDt:.1f}ms"
        ax.plot(df_list[i]["x"], df_list[i]["u"], label=labelList,
                markersize=m_markersize, marker="s", markevery=0.08)

    ax.set_xlim([x_min, x_max])

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 1.0, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax.set_ylabel('u(m/s)')
    ax.set_xlabel('x(m)')

    fig.suptitle(r'u distribution')
    fig.set_size_inches(6, 4)
    fig.tight_layout()
    fig.savefig("u.jpg", dpi=300)