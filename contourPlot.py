"""
==============================================================================
AE588 Project 2D contour plots
==============================================================================
@File    :   contourPlot.py
@Date    :   2021/04/20
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os

# ==============================================================================
# External Python modules
# ==============================================================================
import matplotlib.pyplot as plt
import niceplots
import numpy as np
import pickle
from pyoptsparse import History
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ==============================================================================
# Extension modules
# ==============================================================================

niceplots.setRCParams()
niceColours = niceplots.get_niceColors()
colourList = list(niceColours.values())

with open("FineContour2/ContourResults.pkl", "rb") as file:
    contourData = pickle.load(file)

with open("LineSweeps/LSkinSweepResults.pkl", "rb") as file:
    LSkinSweepData = pickle.load(file)

with open("LineSweeps/USkinSweepResults.pkl", "rb") as file:
    USkinSweepData = pickle.load(file)

fig, ax = plt.subplots(figsize=(11.75 * 1.5, 8.25 * 1.5))
# create new axes on the right and on the top of the current axes
divider = make_axes_locatable(ax)
# below height and pad are in inches
ax_sweepx = divider.append_axes("top", 1.5 * 3, pad=0.1, sharex=ax)
ax_sweepy = divider.append_axes("right", 1.5 * 3, pad=0.4, sharey=ax)
ax_sweepx.xaxis.set_tick_params(labelbottom=False)
ax_sweepy.yaxis.set_tick_params(labelleft=False)

# --- Contour plot of Mass with constrain boundaries ---


ax.set_xlabel("$t_L$ (mm)")
ax.set_ylabel("$t_U$ (mm)", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax, outward=True)
ax.contour(
    1e3 * contourData["tLSkin"],
    1e3 * contourData["tUSkin"],
    contourData["2.5g_TotalMass"],
    levels=30,
    cmap=niceplots.parula_map,
)

ax.contour(
    1e3 * contourData["tLSkin"],
    1e3 * contourData["tUSkin"],
    contourData["2.5g_USkinKSFailure"],
    levels=[1.0],
    colors=niceColours["Orange"],
)
ax.contourf(
    1e3 * contourData["tLSkin"],
    1e3 * contourData["tUSkin"],
    contourData["2.5g_USkinKSFailure"],
    levels=[1.0, np.inf],
    colors=niceColours["Orange"],
    alpha=0.2,
)

ax.contour(
    1e3 * contourData["tLSkin"],
    1e3 * contourData["tUSkin"],
    contourData["2.5g_LSkinKSFailure"],
    levels=[1.0],
    colors=niceColours["Purple"],
)
ax.contourf(
    1e3 * contourData["tLSkin"],
    1e3 * contourData["tUSkin"],
    contourData["2.5g_LSkinKSFailure"],
    levels=[1.0, np.inf],
    colors=niceColours["Purple"],
    alpha=0.2,
)

for data in [LSkinSweepData, USkinSweepData]:
    ax.plot(1e3 * data["tLSkin"], 1e3 * data["tUSkin"], color=niceColours["Grey"])

# --- Contour plot of Stress constraints along upper skin thickness sweep ---


ax_sweepy.set_xlabel("$KS(\sigma)$")
ax_sweepy.set_ylabel("$t_U$\n(mm)", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax_sweepy, spines=["bottom"], outward=True)
ax_sweepy.plot(
    USkinSweepData["2.5g_LSkinKSFailure"],
    1e3 * USkinSweepData["tUSkin"],
    label="$\sigma_L$ Constraint",
    color=niceColours["Purple"],
)
ax_sweepy.plot(
    USkinSweepData["2.5g_USkinKSFailure"],
    1e3 * USkinSweepData["tUSkin"],
    label="$\sigma_U$ Constraint",
    color=niceColours["Orange"],
)
ax_sweepy.axvline(1.0, linestyle="--", color=niceColours["Grey"], alpha=0.6)
ax_sweepy.legend(labelcolor="linecolor", loc="center", bbox_to_anchor=(0.5, 1.35), fontsize=32)

# --- Contour plot of Stress constraints along lower skin thickness sweep ---


ax_sweepx.set_xlabel("$t_L$ (mm)")
ax_sweepx.set_ylabel("$KS(\sigma)$", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax_sweepx, spines=["left"], outward=True)
ax_sweepx.plot(
    1e3 * LSkinSweepData["tLSkin"],
    LSkinSweepData["2.5g_LSkinKSFailure"],
    label="$\sigma_L$ Constraint",
    color=niceColours["Purple"],
)
ax_sweepx.plot(
    1e3 * LSkinSweepData["tLSkin"],
    LSkinSweepData["2.5g_USkinKSFailure"],
    label="$\sigma_U$ Constraint",
    color=niceColours["Orange"],
)
ax_sweepx.axhline(1.0, linestyle="--", color=niceColours["Grey"], alpha=0.6)

# --- Noise Plot ---

fig.savefig("Figures/2dOptProblem.pdf")
plt.show()
# ==============================================================================
# 2d Optimisation results
# ==============================================================================
pathfig, pathaxes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True, figsize=(16, 9))
pathaxes[1, 0].set_xlabel("$t_L$ (mm)")
pathaxes[1, 0].set_ylabel("$t_U$\n(mm)", rotation="horizontal", ha="right")
pathaxes[1, 0].set_xticks([25, 41.75, 70])
pathaxes[1, 0].set_yticks([15, 26.8, 50])
pathaxes = pathaxes.flatten()


histFig, histAx = plt.subplots()
histAx.set_xlabel("Iterations")
histAx.set_ylabel("$M$ (kg)", rotation="horizontal", ha="center")
niceplots.adjust_spines(histAx, outward=True)

altPathFig, altPathAx = plt.subplots()
altPathAx.set_xlabel(r"$\sigma_{max}$")
altPathAx.set_ylabel("$M$ (kg)", rotation="horizontal", ha="center")
niceplots.adjust_spines(altPathAx, outward=True)

optimisers = ["SNOPT", "IPOPT", "NLPQLP", "ParOptSL1", "ParOptFilter", "ParOptMMA"]
optNames = ["SNOPT", "IPOPT", "NLPQLP", "ParOpt S$\mathcal{L}_1$", "ParOpt Filter", "ParOpt MMA"]
xTicks = [0]
for i in range(len(optimisers)):
    niceplots.adjust_spines(pathaxes[i], outward=True)
    pathaxes[i].contour(
        1e3 * contourData["tLSkin"],
        1e3 * contourData["tUSkin"],
        contourData["2.5g_TotalMass"],
        levels=30,
        cmap=niceplots.parula_map,
    )

    pathaxes[i].contour(
        1e3 * contourData["tLSkin"],
        1e3 * contourData["tUSkin"],
        contourData["2.5g_USkinKSFailure"],
        levels=[1.0],
        colors=niceColours["Orange"],
    )
    pathaxes[i].contourf(
        1e3 * contourData["tLSkin"],
        1e3 * contourData["tUSkin"],
        contourData["2.5g_USkinKSFailure"],
        levels=[1.0, np.inf],
        colors=niceColours["Orange"],
        alpha=0.4,
    )

    pathaxes[i].contour(
        1e3 * contourData["tLSkin"],
        1e3 * contourData["tUSkin"],
        contourData["2.5g_LSkinKSFailure"],
        levels=[1.0],
        colors=niceColours["Purple"],
    )
    pathaxes[i].contourf(
        1e3 * contourData["tLSkin"],
        1e3 * contourData["tUSkin"],
        contourData["2.5g_LSkinKSFailure"],
        levels=[1.0, np.inf],
        colors=niceColours["Purple"],
        alpha=0.4,
    )
    histFile = os.path.join("2dOpts", optimisers[i], f"{optimisers[i]}StructOpt.hst")
    hist = History(histFile)
    values = hist.getValues()
    if optimisers[i] == "IPOPT":
        # For some reason, every IPOPT iteration seems to be stored in the history twice so we need to correct this
        for f, vals in values.items():
            values[f] = vals[::2]

    pathaxes[i].plot(
        1e3 * values["tLSkin"],
        1e3 * values["tUSkin"],
        marker="o",
        markeredgecolor="w",
        markersize=6,
        label=optNames[i],
        color=colourList[i],
    )
    pathaxes[i].set_title(optNames[i], color=colourList[i])
    histAx.plot(
        values["2.5g_TotalMass"],
        marker="o",
        markeredgecolor="w",
        markersize=8,
        label=optNames[i],
        clip_on=False,
        zorder=20 - len(values["2.5g_TotalMass"]),
    )
    if len(values["2.5g_TotalMass"]) - 1 not in xTicks:
        xTicks.append(len(values["2.5g_TotalMass"]) - 1)

    altPathAx.plot(
        values["2.5g_TotalMaxFailure"],
        values["2.5g_TotalMass"],
        # marker="o",
        # markeredgecolor="w",
        label=optNames[i],
        clip_on=False,
    )

# pathax.legend(labelcolor="linecolor", frameon=True, framealpha=0.9)
histAx.set_xticks(xTicks)
histAx.set_yticks([values["2.5g_TotalMass"][0], values["2.5g_TotalMass"][-1]])
histAx.legend(labelcolor="linecolor")

pathfig.savefig("Figures/2dOptPaths.pdf")
histFig.savefig("Figures/2dOptHist.pdf")
plt.show()
