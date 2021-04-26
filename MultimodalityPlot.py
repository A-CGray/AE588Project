"""
==============================================================================
AE588 Project Plotting Multimodality test results
==============================================================================
@File    :   MultimodalityPlot.py
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

adjConData = {}
noAdjConData = {}
for i in range(5):
    hist = History(f"MultiModalTest/AdjCon/{i}/SNOPTStructOpt.hst")
    adjConData[i] = hist.getValues()
    if "xAvg" not in adjConData:
        adjConData["xAvg"] = adjConData[i]["struct"][-1] * 1 / 5
    else:
        adjConData["xAvg"] += adjConData[i]["struct"][-1] * 1 / 5

    hist = History(f"MultiModalTest/NoAdjCon/{i}/SNOPTStructOpt.hst")
    noAdjConData[i] = hist.getValues()
    if "xAvg" not in noAdjConData:
        noAdjConData["xAvg"] = noAdjConData[i]["struct"][-1] * 1 / 5
    else:
        noAdjConData["xAvg"] += noAdjConData[i]["struct"][-1] * 1 / 5

# ==============================================================================
# DV Scatter Plot
# ==============================================================================
scatterFig, scatterAx = plt.subplots(nrows=2, figsize=(9, 12))
scatterAx[0].set_xlabel(r"$t_{av}$ (mm)")
scatterAx[1].set_xlabel(r"$t_{av}$ (mm)")
scatterAx[0].set_ylabel(r"$t$ (mm)", rotation="horizontal", ha="right")
scatterAx[1].set_ylabel(r"$t$ (mm)", rotation="horizontal", ha="right")
scatterAx[0].set_yscale("log")
scatterAx[0].set_xscale("log")
niceplots.adjust_spines(scatterAx[0], outward=True)
scatterAx[1].set_yscale("log")
scatterAx[1].set_xscale("log")
niceplots.adjust_spines(scatterAx[1], outward=True)
scatterAx[0].set_title("With Adjacency Constraints")
scatterAx[1].set_title("Without Adjacency Constraints")

scatterAx[0].set_xticks([1e3 * np.min(adjConData["xAvg"]), 1e3 * np.max(adjConData["xAvg"])])
scatterAx[0].set_yticks([1e3 * np.min(adjConData["xAvg"]), 1e3 * np.max(adjConData["xAvg"])])

scatterAx[1].set_xticks([1e3 * np.min(noAdjConData["xAvg"]), 1e3 * np.max(noAdjConData["xAvg"])])
scatterAx[1].set_yticks([1e3 * np.min(noAdjConData["xAvg"]), 1e3 * np.max(noAdjConData["xAvg"])])

scatterAx[0].set_xticklabels([f"{(1e3 * np.min(adjConData['xAvg'])):.01f}", f"{1e3 * np.max(adjConData['xAvg']):.01f}"])
scatterAx[0].set_yticklabels([f"{1e3 * np.min(adjConData['xAvg']):.01f}", f"{1e3 * np.max(adjConData['xAvg']):.01f}"])

scatterAx[1].set_xticklabels(
    [f"{1e3 * np.min(noAdjConData['xAvg']):.01f}", f"{1e3 * np.max(noAdjConData['xAvg']):.01f}"]
)
scatterAx[1].set_yticklabels(
    [f"{1e3 * np.min(noAdjConData['xAvg']):.01f}", f"{1e3 * np.max(noAdjConData['xAvg']):.01f}"]
)

for delta in [20, 10, 5, 1]:
    scatterAx[0].fill_between(
        1e3 * np.sort(adjConData["xAvg"]),
        1e3 * np.sort((1.0 - delta / 100.0) * adjConData["xAvg"]),
        1e3 * np.sort((1.0 + delta / 100.0) * adjConData["xAvg"]),
        alpha=0.3,
        label=f"$\pm {delta}$%",
        clip_on=False,
    )
    scatterAx[1].fill_between(
        1e3 * np.sort(noAdjConData["xAvg"]),
        1e3 * np.sort((1.0 - delta / 100.0) * noAdjConData["xAvg"]),
        1e3 * np.sort((1.0 + delta / 100.0) * noAdjConData["xAvg"]),
        alpha=0.3,
        label=f"$\pm {delta}$%",
        clip_on=False,
    )

for i in range(5):
    scatterAx[0].plot(
        1e3 * adjConData["xAvg"],
        1e3 * adjConData[i]["struct"][-1],
        "+",
        markersize=10,
        markeredgewidth=0.75,
        alpha=0.7,
        clip_on=False,
    )
    scatterAx[1].plot(
        1e3 * noAdjConData["xAvg"],
        1e3 * noAdjConData[i]["struct"][-1],
        "+",
        markersize=10,
        markeredgewidth=0.75,
        alpha=0.7,
        clip_on=False,
    )
scatterAx[0].legend()
scatterFig.savefig("Figures/MultimodalityCorrelation.pdf")

# ==============================================================================
# Cummulative Distribution plots
# ==============================================================================
distFig, distAx = plt.subplots(nrows=2, sharex=True, figsize=(9, 12))
distAx[1].set_xlabel(r"$|\frac{t - t_{av}}{t_{av}}|$")
distAx[0].set_ylabel("Prob", rotation="horizontal", ha="right")
distAx[1].set_ylabel("Prob", rotation="horizontal", ha="right")
# distAx[0].set_yscale("log")
distAx[0].set_xscale("log")
niceplots.adjust_spines(distAx[0], spines=["left"], outward=True)
# distAx[1].set_yscale("log")
distAx[1].set_xscale("log")
niceplots.adjust_spines(distAx[1], outward=True)
distAx[0].set_title("With Adjacency Constraints")
distAx[1].set_title("Without Adjacency Constraints")

bins = 10 ** np.linspace(-10, -2, 100)
for i in range(5):
    relDiff = np.abs((adjConData[i]["struct"][-1] - adjConData["xAvg"]) / adjConData["xAvg"])
    distAx[0].hist(relDiff, bins, density=True, histtype="step", cumulative=True)

    relDiff = np.abs((noAdjConData[i]["struct"][-1] - noAdjConData["xAvg"]) / noAdjConData["xAvg"])
    distAx[1].hist(relDiff, bins, density=True, histtype="step", cumulative=True)
distFig.savefig("Figures/MultimodalityDistributions.pdf")

# ==============================================================================
# Mass and Max Stress path
# ==============================================================================
pathFig, pathAx = plt.subplots(nrows=2, figsize=(9, 12))

pathAx[1].set_xlabel(r"$\sigma_{max}$")
pathAx[0].set_ylabel("$M$ (kg)", rotation="horizontal", ha="right")
pathAx[1].set_ylabel("$M$ (kg)", rotation="horizontal", ha="right")
# pathAx[0].set_yscale("log")
# pathAx[0].set_xscale("log")
niceplots.adjust_spines(pathAx[0], spines=["left"], outward=True)
# pathAx[1].set_yscale("log")
# pathAx[1].set_xscale("log")
niceplots.adjust_spines(pathAx[1], outward=True)
pathAx[0].set_title("With Adjacency Constraints")
pathAx[1].set_title("Without Adjacency Constraints")

for i in range(5):
    pathAx[0].plot(
        adjConData[i]["2.5g_TotalMaxFailure"],
        adjConData[i]["2.5g_TotalMass"],
        # marker="o",
        # markeredgecolor="w",
        clip_on=False,
    )
    pathAx[0].plot(
        adjConData[i]["2.5g_TotalMaxFailure"][0],
        adjConData[i]["2.5g_TotalMass"][0],
        "o",
        markeredgecolor="w",
        color="k",
        clip_on=False,
    )
    pathAx[1].plot(
        noAdjConData[i]["2.5g_TotalMaxFailure"],
        noAdjConData[i]["2.5g_TotalMass"],
        # marker="o",
        # markeredgecolor="w",
        clip_on=False,
    )
    pathAx[1].plot(
        noAdjConData[i]["2.5g_TotalMaxFailure"][0],
        noAdjConData[i]["2.5g_TotalMass"][0],
        "o",
        markeredgecolor="w",
        color="k",
        clip_on=False,
    )

pathAx[0].plot(
    adjConData[i]["2.5g_TotalMaxFailure"][-1],
    adjConData[i]["2.5g_TotalMass"][-1],
    "x",
    # markeredgecolor="w",
    markersize=10,
    color="k",
    clip_on=False,
)
pathAx[1].plot(
    noAdjConData[i]["2.5g_TotalMaxFailure"][-1],
    noAdjConData[i]["2.5g_TotalMass"][-1],
    "x",
    # markeredgecolor="w",
    markersize=10,
    color="k",
    clip_on=False,
)
pathFig.savefig("Figures/MultimodalityPaths.pdf")
plt.show()