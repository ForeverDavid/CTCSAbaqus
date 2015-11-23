import numpy
from abaqus import mdb
import abaqusConstants as aq 
from hAssembly import *
from hCoordinates import *
from hJob import *
from hMaterial import *
from hMesh import *
from hModel import *
from hPart import *
from hProperty import *
from hStep import *

#### execfile('hPart.py')

materials = getMaterialList() # Load in material Data
matrix = "ESBR" # Choose the first experiment from TCNanoFillers
fillers = ["Alumina", "ZincOxide"] # Two fillers associated with first experiment

modelObject, modelName = createModel(2) # Create model database "Model-1"
side, radius, portions, dP, dM, cP, cM = defExperiment(modelObject, matrix, fillers[0]) # Define material attributes for specified matrix, fillers

seed = numpy.random.randint(1000) # Random seed for coordinate generation
radius, number = invPHR(portions[2], dP, radius, dM, side) # Returns specific radius size and number of inclusions for closest PHR value.
xVals, yVals, warningPoints, number = getPoints2D(seed, side, radius, number) # returns coordinates for inclusion locations. 
part = createMatrix(modelObject, side) # Create the matrix
edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our part
matrixSet, fillerSet, allSet = createCircleInclusion(modelObject, part, radius, number, xVals, yVals, side) # Draw inclusions in the matrix
createSection(modelObject, part, matrix, matrixSet) # Create section for matrix material
createSection(modelObject, part, fillers[0], fillerSet) # Create section for filler material
assemblyTop, assemblyBottom = makeAssembly2D(modelObject, part) # Create assembly and return references to assembly sets
temp1, temp2 = 328.15, 298.15 # Assign heat temperature to be used in experiment
heatStep2D(modelObject, assemblyBottom, assemblyTop, temp1, temp2) # apply heat BC
#elements, nodes = makeMesh2D() # Draw mesh and return number of nodes and elements

######
#warningString = submitJob() # Submit job and take note of any warnings
#avgHF, TC = getThermalProperties2D() # Extract relevant information about thermal properties
