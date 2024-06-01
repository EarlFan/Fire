#!/usr/bin/python3.9

import meshio
import numpy as np
# import cantera as ct
import pandas as pd
import math

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.font_manager

T1 = 747.672441027809
p1 = 35.594e3
rho1 = 0.18075
u1 = 487.34

headerList2 = ["x", "rhoMix", "p", "T", "rho_AR", "rho_H2", "rho_H", "rho_OH",
               "rho_O", "rho_H2O", "rho_HO2", "rho_H2O2", "rho_O2", "pMax",
               "u", "v", "w", "gradRho_norm", "tau_c", "tau_r", "lvl", "Xi"]

MW_list = [39.948,2.01588,1.00794,17.00734,15.9994,18.01528,33.00677,34.0147,31.9988]

totalFrames = 24

frameGap = 1

frameDt = 10

df_list = [pd.DataFrame() for x in range(totalFrames)]
df_list_raw = [pd.DataFrame() for x in range(totalFrames)]

parentDir = ""

for i in range(0,totalFrames,):

    fileName = parentDir+"cuts/cut1D_CPU0_TIME"+str(i*frameGap)+".out"

    df_list_raw[i] = pd.read_csv(fileName,skip_blank_lines=1, delim_whitespace=1,names=headerList2)

    df_list[i] = df_list_raw[i].sort_values(by=["x"])


pparam = dict(xlabel='x', ylabel='T')

m_markersize = 3
m_markevery = 0.1

# 180us

header_Samuel= ["x","y"]

rho_180_raw = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/180us_rho.dat",skip_blank_lines=1, delim_whitespace=1)
u_180 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/180us_u.dat",skip_blank_lines=1, delim_whitespace=1)
p_180 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/180us_p.dat",skip_blank_lines=1, delim_whitespace=1)
T_180 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/180us_T.dat",skip_blank_lines=1, delim_whitespace=1)

#230us

rho_230 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/230us_rho.dat",skip_blank_lines=1, delim_whitespace=1)
u_230 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/230us_u.dat",skip_blank_lines=1, delim_whitespace=1)
p_230 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/230us_p.dat",skip_blank_lines=1, delim_whitespace=1)
T_230 = pd.read_csv("../Detonation_1D_post-process/Paolucci_dat/230us_T.dat",skip_blank_lines=1, delim_whitespace=1)

rho_180 = rho_180_raw.sort_values(by=["x"])

fig, ax = plt.subplots(2,2)
# fig.subplots_adjust(right=0.75)
frameIndex = 18

m_markersize = 3
m_every = 3

ax[0,0].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["rhoMix"]/rho1,label="Current",color="k",marker="None")
ax[0,1].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["u"]/u1,color="k",marker="None")
ax[1,0].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["p"]/p1,color="k",marker="None")
ax[1,1].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["T"]/T1,color="k",marker="None")

ax[0,0].plot(rho_180["x"], rho_180["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[0,1].plot(u_180["x"], u_180["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[1,0].plot(p_180["x"], p_180["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[1,1].plot(T_180["x"], T_180["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')


ax[0,0].legend(fontsize=6,loc='upper left')
# ax.legend(bbox_to_anchor=(0.9, 0.5, 0.5, 0.5))

# rho
# ax[0,0].autoscale(tight=True)
ax[0,0].set_xlim([0,1])
# ax[0,0].set_ylim([0,5])
ax[0,0].set_yticks([0,1,2,3,4,5])
ax[0,0].set_ylabel(r'$\rho/\rho_1$')
# ax[0,0].set_xlabel(r'$x/L$')

# velocity
ax[0,1].set_xlim([0,1])
ax[0,1].set_ylim([-1.5,1.5])
ax[0,1].set_yticks([-1.5,-1,-0.5,0,0.5,1,1.5])
ax[0,1].set_ylabel(r'$u/u_1$')
# ax[0,1].set_xlabel(r'$x/L$')

# Pressure
ax[1,0].set_xlim([0,1])
ax[1,0].set_ylim([0,15])
ax[1,0].set_yticks([0,5,10,15])
ax[1,0].set_ylabel(r'$p/p_1$')
ax[1,0].set_xlabel(r'$x/L$')

# Temperature
ax[1,1].set_xlim([0,1])
ax[1,1].set_ylim([0,4])
ax[1,1].set_yticks([0,1,2,3,4])
ax[1,1].set_ylabel(r'$T/T_1$')
ax[1,1].set_xlabel(r'$x/L$')



fig.set_size_inches(4.75, 4)
fig.tight_layout()
fig.savefig("180us-5L.jpg", dpi=600)



fig, ax = plt.subplots(2,2)
# fig.subplots_adjust(right=0.75)
frameIndex = 23

ax[0,0].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["rhoMix"]/rho1,label="Current study",color="k",marker="None")
ax[0,1].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["u"]/u1,color="k",marker="None")
ax[1,0].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["p"]/p1,color="k",marker="None")
ax[1,1].plot(df_list[frameIndex]["x"]/0.12, df_list[frameIndex]["T"]/T1,color="k",marker="None")

ax[0,0].plot(rho_230["x"], rho_230["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[0,1].plot(u_230["x"], u_230["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[1,0].plot(p_230["x"], p_230["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')
ax[1,1].plot(T_230["x"], T_230["y"],label="Paolucci et al.",linestyle="None",color="b",markevery=m_every,marker="s",markersize=m_markersize,markerfacecolor='none')


ax[0,0].legend(fontsize=6)
# ax.legend(bbox_to_anchor=(0.9, 0.5, 0.5, 0.5))

# rho
ax[0,0].set_xlim([0,1])
ax[0,0].set_ylim([0,4])
ax[0,0].set_ylabel(r'$\rho/\rho_1$')
# ax[0,0].set_xlabel(r'$x/L$')

# velocity
ax[0,1].set_xlim([0,1])
ax[0,1].set_ylim([-1.5,1.5])
ax[0,1].set_yticks([-1.5,-1,-0.5,0,0.5,1,1.5])
ax[0,1].set_ylabel(r'$u/u_1$')
# ax[0,1].set_xlabel(r'$x/L$')

# Pressure
ax[1,0].set_xlim([0,1])
ax[1,0].set_ylim([0,12])
ax[1,0].set_yticks([0,4,8,12])
# ax[1,0].set_yticks([0,5,10,15])
ax[1,0].set_ylabel(r'$p/p_1$')
ax[1,0].set_xlabel(r'$x/L$')

# Temperature
ax[1,1].set_xlim([0,1])
ax[1,1].set_ylim([0,4])
ax[1,1].set_yticks([0,1,2,3,4])
ax[1,1].set_ylabel(r'$T/T_1$')
ax[1,1].set_xlabel(r'$x/L$')



fig.set_size_inches(4.75, 4)
fig.tight_layout()
fig.savefig("230us-5L.jpg", dpi=600)

fig, ax = plt.subplots(2,1)

m_markersize=1

ax[0].plot(df_list[18]["x"]/0.12, df_list[18]["lvl"],label=r'$180 \mu s$',color="b",linestyle = "None",marker="o",markersize=m_markersize,markerfacecolor='b')
ax[1].plot(df_list[23]["x"]/0.12, df_list[23]["lvl"],label=r'$230 \mu s$',color="b",linestyle = "None",marker="o",markersize=m_markersize,markerfacecolor='b')


ax[0].set_xlim([0,1])
ax[0].set_xticks([])
ax[0].set_ylim([0,6])
ax[0].set_yticks([0,2,4,6])
ax[0].set_ylabel(r'Level')
ax[0].minorticks_off()
# ax[0].set_xlabel(r'$x/L$')

# velocity
# ax[0,1].autoscale(tight=True)
ax[1].set_xlim([0,1])
ax[1].set_ylim([0,6])
ax[1].set_yticks([0,2,4,6])
ax[1].set_ylabel(r'Level')
ax[1].set_xlabel(r'$x/L$')
ax[1].minorticks_off()

# fig.suptitle(r'Unsteady reaction fronts')
# figureName = "TIME-"+str(i)+".jpg"
fig.set_size_inches(4, 3)
fig.tight_layout()
fig.savefig("level-5L.jpg", dpi=600)
# fig.savefig("level-8L.tiff", dpi=600)

