"""
==============================================================================
AE588 Project Main Runscript
==============================================================================
@File    :   struct_opt.py
@Date    :   2021/04/20
@Author  :   Alasdair Christison Gray
@Description :
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import argparse
import ast
import os
import time
from pprint import pprint as pp

# ==============================================================================
# External Python modules
# ==============================================================================
from mpi4py import MPI
from baseclasses import StructProblem
from pyoptsparse import Optimization, OPT, History
from SETUP import setup_tacs
import matplotlib.pyplot as plt
import numpy as np
import pickle

# ==============================================================================
# Extension modules
# ==============================================================================
from OptimiserOptions import getOptOptions


parser = argparse.ArgumentParser()
parser.add_argument(
    "--task", type=str, choices=["analysis", "opt", "contour", "2dOpt", "line", "noiseTest"], default="opt"
)
parser.add_argument("--output", type=str, default="debug")
parser.add_argument("--optOptions", type=ast.literal_eval, default={}, help="additional optimizer options to be added")
parser.add_argument("--noPlot", action="store_true", help="disable plotting")
parser.add_argument("--saveOptIters", type=str, choices=["major", "all", "none"], default="none")
parser.add_argument("--adjCon", type=float, default=-1)
parser.add_argument("--gLoad", action="store_true", help="apply gravity loads")
parser.add_argument("--randomX0", action="store_true", help="Use random initial design")
parser.add_argument("--nContour", type=int, default=101)
parser.add_argument(
    "--optimiser",
    type=str,
    choices=["SNOPT", "IPOPT", "NLPQLP", "ParOptSL1", "ParOptFilter", "ParOptMMA"],
    default="SNOPT",
)
parser.add_argument("--nSkin", type=int, default=360)
parser.add_argument("--nSpar", type=int, default=36)
parser.add_argument("--nRib", type=int, default=18)
parser.add_argument("--optIter", type=int, default=400)
parser.add_argument("--hessianUpdate", type=int, default=1000)
parser.add_argument("--feasibility", type=float, default=1e-6)
parser.add_argument("--optimality", type=float, default=1e-5)

args = parser.parse_args()

comm = MPI.COMM_WORLD

if comm.rank == 0:
    os.system(f"mkdir -p {args.output}")

# ==============================================================================
#       Initialize TACS
# ==============================================================================

structOptions = {
    "transferSize": 1e-1,
    "gravityVector": [0, -9.81, 0],
    "projectVector": [0, 1, 0],
    "outputDir": args.output,
}
FEASolver, dvDictMap = setup_tacs.setup(
    meshFile="geometry/wingbox.bdf",
    options=structOptions,
    nSkinGroups=args.nSkin,
    nSparGroups=args.nSkin,
    nRibGroups=args.nRib,
)

# Write out components to visualize design variables (contour of dv1)
FEASolver.writeDVVisualization(os.path.join(args.output, "dv_visual.f5"))

# ==============================================================================
#       Set up structural problem
# ==============================================================================

dispFuncs = FEASolver.functionList
objConFuncs = ["TotalMass"] + [f for f in dispFuncs if "KSFailure" in f]
sp = StructProblem("2.5g", loadFile="SETUP/forces.txt", loadFactor=2.5, evalFuncs=dispFuncs)

# Add inertial (gravity) loads
if args.gLoad:
    FEASolver.addInertialLoad(sp)

# Add some adjacency constraints between skin and spar thicknesses
if args.adjCon > 0.0:
    FEASolver.addAdjacencyConstraints(delta=args.adjCon, include=["SKIN", "SPAR"])


# ==============================================================================
#       Set up optimization
# ==============================================================================


def obj(x):
    funcs = {}
    FEASolver.setDesignVars(x)
    FEASolver(sp)
    FEASolver.evalFunctions(sp, funcs)
    if comm.rank == 0:
        print("\n\n\n--------")
        print("RESULTS:")
        print("--------")
        for funcType in ["Mass", "KSFailure", "MaxFailure", "Compliance"]:
            print(f"\n{funcType} Functions:")
            print("-----------------------------")
            for f in funcs:
                if funcType in f:
                    print(f"{f:30s}: {funcs[f]:6f}")
        print("-------------------------------------------\n\n\n")

    if args.saveOptIters.lower() == "all":
        FEASolver.writeSolution(outputDir=args.output)
    elif FEASolver.curSP.tacsData.callCounter == 0:
        FEASolver.writeSolution(baseName="%s_init" % sp.name)
    return funcs


def sens(x, funcs):
    funcsSens = {}
    FEASolver.evalFunctionsSens(sp, funcsSens, evalFuncs=objConFuncs)

    if args.saveOptIters.lower() == "major":
        FEASolver.writeSolution(outputDir=args.output)
    return funcsSens


if "opt" in args.task.lower():

    if args.task == "2dOpt":

        def obj2d(x):
            if comm.rank == 0:
                print(x)
            xCopy = {"struct": np.array([x["tUSkin"], x["tLSkin"], 3e-2, 3e-2, 3e-3])}
            return obj(xCopy)

        def sens2d(x, funcs):
            xCopy = {"struct": np.array([x["tUSkin"], x["tLSkin"], 3e-2, 3e-2, 3e-3])}
            funcsSens = sens(xCopy, funcs)
            funcsSens2d = {}
            for f in funcsSens:
                funcsSens2d[f] = {"tUSkin": funcsSens[f]["struct"][0], "tLSkin": funcsSens[f]["struct"][1]}
            return funcsSens2d

        objFun = obj2d
        sensFun = sens2d
    else:
        objFun = obj
        sensFun = sens

    # Set up the optimization problem
    optProb = Optimization("Mass minimization", objFun)
    optProb.addObj("2.5g_TotalMass", scale=1e-3)
    if "2d" not in args.task:
        FEASolver.addVariablesPyOpt(optProb)
        FEASolver.addConstraintsPyOpt(optProb)
    else:
        optProb.addVar("tUSkin", lower=2e-3, upper=1e-1, scale=1e2)
        optProb.addVar("tLSkin", lower=2e-3, upper=1e-1, scale=1e2)
        optProb.setDVs({"tUSkin": 45e-3, "tLSkin": 65e-3})

    # Paropt uses an L1 norm for optimality so we need to do a rough scaling to try and make the L1 optimality criterion roughly equivalent to the L2 criteria
    optimalityScale = np.sqrt(FEASolver.getNumDesignVars()) if "paropt" in args.optimiser.lower() else 1.0
    optOptions = getOptOptions(
        alg=args.optimiser,
        outputDir=args.output,
        majIterLimit=args.optIter,
        hessianUpdate=args.hessianUpdate,
        feasibility=args.feasibility,
        optimality=args.optimality * optimalityScale,
    )
    optOptions.update(args.optOptions)

    if "paropt" in args.optimiser.lower():
        alg = "paropt"
    else:
        alg = args.optimiser.lower()
    opt = OPT(alg, options=optOptions)

    for f in FEASolver.functionList:
        if "KSFailure" in f and "total" not in f.lower() and ("skin" in f.lower() or "2d" not in args.task.lower()):
            optProb.addCon(f"{sp.name}_{f}", upper=1.0)

    if comm.rank == 0:
        print(optProb)
    optProb.printSparsity()

    if args.randomX0:
        minThick = 2e-3
        maxThick = 60e-3
        randomX = minThick + (maxThick - minThick) * np.random.rand(len(FEASolver.scaleList))
        optProb.setDVs({"struct": randomX})

    sol = opt(optProb, sens=sensFun, storeHistory=os.path.join(args.output, f"{args.optimiser}StructOpt.hst"))

    # Write the final solution
    FEASolver.writeOutputFile(os.path.join(args.output, sp.name + "_final.f5"))

    # Save the structural DVs in a dictionary so we can more easily use it
    x = np.zeros(FEASolver.getNumDesignVars())
    FEASolver.structure.getDesignVars(x)

    if comm.rank == 0:
        dvDict = {}
        for cd in dvDictMap:
            dvDict[cd] = {}
            for dv in dvDictMap[cd]:
                idv = dvDictMap[cd][dv]
                dvDict[cd][dv] = x[idv]
        varFile = os.path.join(args.output, "struct_vars.pkl")
        output = open(varFile, "wb")
        pickle.dump(dvDict, output)
        pp(dvDict)
        output.close()

    if comm.rank == 0 and not args.noPlot and False:
        # load the optimisation history
        optHist = History(os.path.join(args.output, f"{args.optimiser}StructOpt.hst"))
        histValues = optHist.getValues()

        # Plot Feasibility, optimality, wingbox mass and failure constraints over optimisation history
        fig, axes = plt.subplots(nrows=3, sharex=True, constrained_layout=True, figsize=(14, 10))
        try:
            axes[0].plot("nMajor", "optimality", data=histValues, label="Optimality")
            axes[0].plot("nMajor", "feasibility", data=histValues, label="Feasibility")
            axes[0].set_yscale("log")
            axes[0].axhline(1e-6, linestyle="--", color="gray")
            axes[0].annotate("Convergence Criteria", xy=(10, 1e-6), ha="left", va="top", fontsize=24, color="gray")
            axes[0].legend(fontsize=20, labelcolor="linecolor", loc="upper right", frameon=False)
            axes[0].autoscale(enable=True, tight=True)
        except:
            pass

        for f in histValues.keys():
            if "mass" in f.lower():
                legName = f.replace("2.5g_", "").replace("Mass", "")
                axes[1].plot("nMajor", f, data=histValues, label=legName)
            if "ksfailure" in f.lower():
                legName = f.replace("2.5g_", "").replace("KSFailure", "")
                axes[2].plot("nMajor", f, data=histValues, label=legName)

        axes[1].legend(fontsize=20, labelcolor="linecolor", loc="upper right", frameon=False, ncol=3)
        axes[2].legend(fontsize=20, labelcolor="linecolor", loc="upper right", frameon=False, ncol=3)
        axes[1].autoscale(enable=True, tight=True)
        axes[2].autoscale(enable=True, tight=True)

        axes[1].set_ylabel("Mass\n(kg)", rotation="horizontal", ha="right", fontsize=24)
        axes[2].set_ylabel("Failure\nConstraint", rotation="horizontal", ha="right", fontsize=24)
        axes[2].set_xlabel("Major Iterations", fontsize=24)
        axes[2].axhline(1, linestyle="--", color="gray")

        for ax in axes:
            ax.tick_params(axis="both", labelsize=20)

        plt.show()

if args.task == "contour":

    uSkinThickness = np.linspace(1.5e-2, 5e-2, args.nContour)
    lSkinThickness = np.linspace(2.5e-2, 7e-2, args.nContour)
    tUSkin, tLSkin = np.meshgrid(uSkinThickness, lSkinThickness)
    contourResults = {"tUSkin": tUSkin, "tLSkin": tLSkin}
    for f in FEASolver.functionList:
        contourResults[f"{sp.name}_{f}"] = np.zeros_like(tUSkin)

    for i in range(args.nContour):
        for j in range(args.nContour):
            if comm.rank == 0:
                print(f"Upper Skin Thickness = {1e3*tUSkin[i, j]:.02f} mm")
                print(f"Lower Skin Thickness = {1e3*tLSkin[i, j]:.02f} mm")
            x = {"struct": np.array([tUSkin[i, j], tLSkin[i, j], 3e-2, 3e-2, 3e-3])}
            startTime = time.time()
            funcs = obj(x)
            evalTime = time.time() - startTime
            if comm.rank == 0:
                print(f"\nAnalysis time = {evalTime}")
            for key in funcs:
                contourResults[key][i, j] = funcs[key]
    with open(os.path.join(args.output, f"ContourResults.pkl"), "wb") as file:
        pickle.dump(contourResults, file, protocol=-1)

if args.task == "line":

    tUSkin = np.linspace(1.5e-2, 5e-2, args.nContour)
    tLSkin = np.linspace(2.5e-2, 7e-2, args.nContour)

    contourResults = {"tUSkin": tUSkin, "tLSkin": 50e-3 * np.ones_like(tUSkin)}
    for f in FEASolver.functionList:
        contourResults[f"{sp.name}_{f}"] = np.zeros_like(tUSkin)

    for i in range(args.nContour):
        if comm.rank == 0:
            print(f"Upper Skin Thickness = {1e3*tUSkin[i]:.02f} mm")
            print("Lower Skin Thickness = 50 mm")
        x = {"struct": np.array([tUSkin[i], 50e-3, 3e-2, 3e-2, 3e-3])}
        funcs = obj(x)
        for key in funcs:
            contourResults[key][i] = funcs[key]
    with open(os.path.join(args.output, f"USkinSweepResults.pkl"), "wb") as file:
        pickle.dump(contourResults, file, protocol=-1)

    contourResults = {"tUSkin": 35e-3 * np.ones_like(tLSkin), "tLSkin": tLSkin}
    for f in FEASolver.functionList:
        contourResults[f"{sp.name}_{f}"] = np.zeros_like(tLSkin)
    for i in range(args.nContour):
        if comm.rank == 0:
            print(f"Upper Skin Thickness = 35 mm")
            print(f"Lower Skin Thickness = {1e3*tLSkin[i]:.02f} mm")
        x = {"struct": np.array([35e-3, tLSkin[i], 3e-2, 3e-2, 3e-3])}
        funcs = obj(x)
        for key in funcs:
            contourResults[key][i] = funcs[key]
    with open(os.path.join(args.output, f"LSkinSweepResults.pkl"), "wb") as file:
        pickle.dump(contourResults, file, protocol=-1)

# if args.task == "noiseTest":