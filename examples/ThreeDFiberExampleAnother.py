import numpy
import math
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
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

materials = getMatListCNT()
matrix = 'ESBR'
filler = 'CNT'

modelObject, modelName = createModel(2)
side, radius, length, portions, dP, dM, cP, cM = defExperimentFiber(modelObject, matrix, filler)

part = createMatrix(modelObject, side, False)
edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our matrix
part2 = createFiber(modelObject, radius, length) # Create Fiber part
edges2, vertices2, face2 = assignGeomSequence(part2) # Create references to important sets in fiber 
matrixSet, fiberSet = create3DInitialSetsFibers(part, part2, side)

createSection(modelObject, part, matrix, matrixSet) # Create section for matrix material
createSection(modelObject, part2, filler, fiberSet) # Create section for filler material

warningPoints = ""

number = calcNumberFibers(radius, length, side, portions[0])
#cylinders, warningMsg, newnum = get3DCylinders(10, side, radius, length, number)
cylinders, warningMsg, newnum = get3DCylindersBrute(10, side, radius, length, number)

modelRootAssembly, fullMatrixPart = create3DMatrixFiberInclusions(modelObject, part, part2, cylinders) # Create assembly and return references to assembly sets

