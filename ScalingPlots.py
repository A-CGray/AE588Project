"""
==============================================================================
AE588 Project Plotting scaling test results
==============================================================================
@File    :   ScalingPlots.py
@Date    :   2021/04/23
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
import numpy as np
import matplotlib.pyplot as plt
import niceplots
from pyoptsparse import History

# ==============================================================================
# Extension modules
# ==============================================================================

niceplots.setRCParams()
colourList = list(niceplots.get_niceColors().values())

# ==============================================================================
# Get history data
# ==============================================================================
numDVs = [5, 45, 90, 162, 234, 450, 810]
optimisers = ["SNOPT", "IPOPT", "NLPQLP", "ParOptSL1", "ParOptFilter"]
optNames = ["SNOPT", "IPOPT", "NLPQLP", "ParOpt S$\mathcal{L}_1$", "ParOpt Filter"]
data = {}

for i in range(len(optimisers)):
    data[optimisers[i]] = {}
    if "paropt" in optimisers[i].lower():
        n = 5
    else:
        n = 7
    for j in range(n):
        hist = History(f"scalingTest/{optimisers[i]}/{j}/{optimisers[i]}StructOpt.hst")
        varsToGet = hist.getObjNames() + hist.getConNames() + ["struct"]
        data[optimisers[i]][j] = hist.getValues(varsToGet, major=False)
        if optimisers[i] == "IPOPT":
            # For some reason, every IPOPT iteration seems to be stored in the history twice so we need to correct this
            for f, vals in data[optimisers[i]][j].items():
                data[optimisers[i]][j][f] = vals[::2]

# ==============================================================================
# Make some plots
# ==============================================================================

# ==============================================================================
# Plot 1, number of func evals vs number of DV's
# ==============================================================================
scalingFig, scalingAx = plt.subplots()
scalingAx.set_xlabel("$N_x$")
scalingAx.set_ylabel("Function Evaluations")
scalingAx.set_yscale("log")
scalingAx.set_xscale("log")
scalingAx.set_xticks(numDVs)
scalingAx.set_xticklabels(numDVs)
niceplots.adjust_spines(scalingAx, outward=True)


for i in range(len(optimisers)):
    fEvals = []
    numOpts = len(data[optimisers[i]])
    for j in range(numOpts):
        fEvals.append(len(data[optimisers[i]][j]["2.5g_TotalMass"]))
    scalingAx.plot(numDVs[:numOpts], fEvals, "-o", markeredgecolor="w", markersize=10, clip_on=False, label=optNames[i])

# scalingAx.plot(numDVs, numDVs, "--", color="gray", clip_on=False, label="$N_x$ Scaling")
scalingAx.plot(numDVs, 10 * np.sqrt(np.array(numDVs)), "--", color="gray", clip_on=False, label="$\sqrt{N_x}$ Scaling")

scalingAx.legend(labelcolor="linecolor")
scalingFig.savefig("FuncEvalScaling.pdf")

# ==============================================================================
# Plot 2 & 3, Mass and failure histories from level 5 optimisation, 1d and 2d views
# ==============================================================================
histFig, histAx = plt.subplots(nrows=2, sharex=True)
histAx[1].set_xlabel("Function Evaluations")
histAx[0].set_ylabel("$M$ (kg)", rotation="horizontal", ha="right")
# histAx[1].set_xscale("log")
histAx[1].set_ylabel(r"$\sigma_{max}$", rotation="horizontal", ha="right")
# histAx[1].set_yscale("log")

niceplots.adjust_spines(histAx[0], spines=["left"], outward=True)
niceplots.adjust_spines(histAx[1], outward=True)

pathFig, pathAx = plt.subplots()
pathAx.set_xlabel("$M$ (kg)")
pathAx.set_ylabel(r"$\sigma_{max}$", rotation="horizontal", ha="right")
pathAx.set_yscale("log")

niceplots.adjust_spines(pathAx, outward=True)


for i in range(len(optimisers)):
    dat = data[optimisers[i]][4]
    histAx[0].plot(dat["2.5g_TotalMass"], clip_on=False, color=colourList[i], zorder=300 - len(dat["2.5g_TotalMass"]))

    maxFailure = np.zeros(len(dat["2.5g_TotalMass"]))
    for f in dat:
        if "KS" in f:
            KSFunc = dat[f].flatten()
            maxFailure = np.where(maxFailure > KSFunc, maxFailure, KSFunc)
    histAx[1].plot(
        maxFailure, label=optNames[i], clip_on=False, color=colourList[i], zorder=300 - len(dat["2.5g_TotalMass"])
    )
    pathAx.plot(
        dat["2.5g_TotalMass"],
        np.where(maxFailure < 1.000000001, 1.000000001, maxFailure - 1.0),
        clip_on=False,
        color=colourList[i],
    )


histAx[1].legend(labelcolor="linecolor")
pathAx.legend(labelcolor="linecolor")
histFig.savefig("ScalingHistories.pdf")
plt.show()