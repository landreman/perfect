from p_1d_plot import perfect_1d_plot
from p_2d_plot import perfect_2d_plot

import sys
from scipy.constants import pi


dirlist=["../2","."]

dirlist=[x[:-1] if x[-1]=='/' else x for x in dirlist ]
#vlines2=[94.927395957025573, 97.463697978512787]
vlines=[0.94927395957025573, 0.97463697978512787]
vlines2=[0.94927395957025573, 0.97463697978512787]


xlims = [0.77,1.04]

colors=['blue',"red"]
linestyles= ['solid',"solid"]



attribs=["particle_source_over_m2nPed","heat_source_over_m2nPed","momentum_source_over_m52nPed"]
ylabels=[r"$\hat{S}_p/n_{0}$",r"$\hat{S}_h/n_{0}$",r"$\hat{S}_m/n_{0}$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="source",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)

attribs=["charge_source"]
ylabels=[r"$\sum Z\hat{S}_p$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="charge_source",vlines=vlines,hlines=[0],generic_labels=False,xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)

attribs=["particle_flux_over_nPed","conductive_heat_flux_over_nPed","m_momentum_flux_over_nPed"]
ylabels=[r"$\hat{\Gamma}/n_{0}$",r"$\hat{q}/n_{0}$",r"$\hat{m}\hat{\Pi}/n_{0}$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="flux",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)

attribs=["jHat","conductive_heat_flux_sum","momentum_flux_sum"]
ylabels=[r"$\sum_a Z_a \hat{\Gamma}_a$",r"$\sum_a \hat{q}_a$",r"$\sum_a \hat{m}_a\hat{\Pi}_a$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="total_flux",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)


#TODO: FSAFlow is currently not an output
attribs=["FSAkPar","FSAFlow"]
ylabels=[r"$\langle k_\| \rangle$",r"$\langle V_\| \rangle$"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="FSAkParFlow",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)

attribs=["toroidal_flow_outboard","toroidal_flow_inboard"]
ylabels=[r"$\hat{V}_t$ Outboard",r"$\hat{V}_t$ Inboard"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="VtInOut",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)

attribs=["poloidal_flow_outboard","poloidal_flow_inboard"]
ylabels=[r"$\hat{V}_p$ Outboard",r"$\hat{V}_p$ Inboard"]
perfect_1d_plot(dirlist,attribs,same_plot=True,ylabels=ylabels,outputname="VpInOut",vlines=vlines,hlines=[0],generic_labels=False,share_scale=["D","He"],xlims=xlims,xattr="actual_psiN",colors=colors,linestyles=linestyles,lg=False)
