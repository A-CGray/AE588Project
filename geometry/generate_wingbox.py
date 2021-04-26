# rst begin import
# ==============================================================================
#       Import modules
# ==============================================================================
import numpy
from pygeo import pyGeo, geo_utils
from pylayout import pyLayout

# rst end import
# ==============================================================================
#       Specify wingbox properties
# ==============================================================================
chords = [5, 1.25]  # root and tip chords
sweep = [0, 7.5]  # root and tip sweep
span = 14.001  # wingspan
sob = 1.5  # span location of side-of-body
ncols = 37  # number of columns (aligned with ribs)
nrows = 10  # number of rows (aligned with spars)
nbreak = 5  # column index of side-of-body kink
# rst start refinement
# Number of quad elements in each component
colSpace = 4 * numpy.ones(ncols - 1, "intc")  # elements between columns
rowSpace = 3 * numpy.ones(nrows + 1, "intc")  # elements between rows
verticalSpace = 6  # elements between skins
ribStiffenerSpace = 2  # elements along rib stiffener depth
stringerSpace = 2  # elements along stringer depth
# rst end refinement
# ==============================================================================
#       Set up blanking arrays
# ==============================================================================
# rst start ribBlank
# Blanking for ribs
ribBlank = numpy.ones((ncols, nrows - 1), "intc")
ribBlank[0, :] = 0  # Blank root rib
ribBlank[1::2, :] = 0
# rst end ribBlank
# Blanking for spars
sparBlank = numpy.zeros((nrows, ncols - 1), "intc")
sparBlank[0, :] = 1  # Keep First
sparBlank[-1, :] = 1  # Keep Last

# Blanking for upper skin stringers:
topStringerBlank = numpy.zeros((nrows, ncols - 1), "intc")
# topStringerBlank[:, :] = 1  # NO Blanking

# Blanking for lower skin stringers:
botStringerBlank = numpy.zeros((nrows, ncols - 1), "intc")
# botStringerBlank[:, :] = 1  # NO Blanking

# Blanking for rib stiffeners:
ribStiffenerBlank = numpy.zeros((ncols, nrows), "intc")  # No rib stiffeners
teEdgeList = []
# rst end blanking
# ==============================================================================
#       Set up array of grid coordinates for ribs, spars, and stringers
# ==============================================================================
leList = [
    [sweep[0] + 0.25 * chords[0], 0, 0.01],
    [sweep[0] + 0.25 * chords[0], 0, sob],
    [sweep[1] + 0.15 * chords[1], 0, span],
]

teList = [
    [sweep[0] + 0.80 * chords[0], 0, 0.01],
    [sweep[0] + 0.80 * chords[0], 0, sob],
    [sweep[1] + 0.75 * chords[1], 0, span],
]

# Initialize grid coordinate matrix
X = numpy.zeros((ncols, nrows, 3))

# Fill in LE and TE coordinates from root to side-of-body
X[0:nbreak, 0] = geo_utils.linearEdge(leList[0], leList[1], nbreak)
X[0:nbreak, -1] = geo_utils.linearEdge(teList[0], teList[1], nbreak)

# Fill in LE and TE coordinates from side-of-body to tip
X[nbreak - 1 : ncols, 0] = geo_utils.linearEdge(leList[1], leList[2], ncols - nbreak + 1)
X[nbreak - 1 : ncols, -1] = geo_utils.linearEdge(teList[1], teList[2], ncols - nbreak + 1)

# Finally fill in chord-wise with linear edges
for i in range(ncols):
    X[i, :] = geo_utils.linearEdge(X[i, 0], X[i, -1], nrows)

# Boundary conditions
ribBC = {
    0: {"all": "345"},
    4: {"edge": "12"},
}
# rst end setup
# ==============================================================================
#       Generate wingbox
# ==============================================================================
# Get surface definition to use for projections
surfFile = "wing.igs"
geo = pyGeo("iges", fileName=surfFile)

# Initialize pyLayout
layout = pyLayout.Layout(
    geo,
    teList=[],
    nribs=ncols,
    nspars=nrows,
    elementOrder=2,
    X=X,
    ribBlank=ribBlank,
    sparBlank=sparBlank,
    topStringerBlank=topStringerBlank,
    botStringerBlank=botStringerBlank,
    ribStiffenerBlank=ribStiffenerBlank,
    minStringer_height=0.025,
    maxStringer_height=0.025,
    spanSpace=colSpace,
    ribSpace=rowSpace,
    vSpace=verticalSpace,
    stringerSpace=stringerSpace,
    ribStiffenerSpace=ribStiffenerSpace,
    flipRibStiffner=False,
    flipUp=False,
    ribBC=ribBC,
)
# rst begin writefiles
# Write bdf file
layout.finalize("wingbox.bdf")

# Write a tecplot file
layout.writeTecplot("wingbox.dat")
# rst end file
