import pytacs
from tacs_orig import functions, constitutive
import numpy as np


def setup(meshFile, options, nSparGroups=1, nSkinGroups=1, nRibGroups=1):

    FEASolver = pytacs.pyTACS(meshFile, options=options)

    # rst end setup_tacs-init

    # ==============================================================================
    #       Set up design variable groups
    # ==============================================================================
    dvGroups = {
        "USkin": {"compName": "U_SKIN", "numDV": nSkinGroups},
        "LSkin": {"compName": "L_SKIN", "numDV": nSkinGroups},
        "LESpar": {"compName": "SPAR.00", "numDV": nSparGroups},
        "TESpar": {"compName": "SPAR.09", "numDV": nSparGroups},
        "Ribs": {"compName": "RIB", "numDV": nRibGroups},
    }

    for group, info in dvGroups.items():
        FEASolver.addDVGroup(userDVName=group, include=info["compName"], nGroup=info["numDV"])

    # ==============================================================================
    #       Set-up constitutive properties for each DVGroup
    # ==============================================================================
    def conCallBack(dvNum, compDescripts, userDescript, specialDVs, **kargs):

        # Define Aluminium material properties and shell thickness limits
        rho_2024 = 2780  # Density, kg/m^3
        E_2024 = 73.1e9  # Elastic Modulus, Pa
        ys_2024 = 324e6  # Yield Strength, Pa
        nu = 0.33  # Poisson's ratio
        kcorr = 5.0 / 6.0  # Shear correction factor for isotropic shells
        tMin = 0.002
        tMax = 0.10

        # Set shell thickness depending on component type
        if "RIB" in userDescript.upper():
            t = 0.004
        else:
            t = 0.02

        # Create constitutive object
        con = constitutive.isoFSDTStiffness(rho_2024, E_2024, nu, kcorr, ys_2024, t, dvNum, tMin, tMax)
        scale = [100.0]
        return con, scale

    # rst end conCallBack

    # TACS will now create the constitutive objects for the different DV groups
    FEASolver.createTACSAssembler(conCallBack)

    dvDictMap = {}
    for i in range(len(FEASolver.DVGroups)):
        for compNum in FEASolver.DVGroups[i]:
            dvDictMap[FEASolver.compDescript[compNum]] = {
                "skin_thickness": i,
            }

    # rst end createAssembler

    # ==============================================================================
    #       Add functions
    # ==============================================================================

    funcGroups = {
        "USkin": "U_SKIN",
        "LSkin": "L_SKIN",
        "LESpar": "SPAR.00",
        "TESpar": "SPAR.09",
        "Rib": "RIBS",
    }
    safetyFactor = 1.5
    KSWeight = 100.0

    # Add mass, compliance and failure functions for the whole wing and for each component group
    FEASolver.addFunction("TotalMass", functions.StructuralMass)
    FEASolver.addFunction("TotalCompliance", functions.Compliance)
    FEASolver.addFunction(
        "TotalKSFailure",
        functions.AverageKSFailure,
        KSWeight=KSWeight,
        safetyFactor=safetyFactor,
    )
    FEASolver.addFunction(
        "TotalMaxFailure",
        functions.AverageMaxFailure,
        safetyFactor=safetyFactor,
    )
    for group, info in dvGroups.items():
        FEASolver.addFunction(group + "Mass", functions.StructuralMass, include=info["compName"])

        FEASolver.addFunction(
            group + "MaxFailure", functions.AverageMaxFailure, include=info["compName"], safetyFactor=safetyFactor
        )

        FEASolver.addFunction(
            group + "KSFailure",
            functions.AverageKSFailure,
            include=info["compName"],
            KSWeight=KSWeight,
            safetyFactor=safetyFactor,
        )

        FEASolver.addFunction(group + "Compliance", functions.Compliance, include=info["compName"])

    return FEASolver, dvDictMap
