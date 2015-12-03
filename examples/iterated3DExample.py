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


trialsPer = 1 # Number of times to run each experiment
interfaceSize = 0.0 # No interface area


materials = getMaterialList() # Load in material Data
matrix = "ESBR" # Choose the first experiment from TCNanoFillers
fillers = ["Alumina", "ZincOxide"] # Two fillers associated with first experiment

fileName = matrix+fillers[0] + fillers[1] # Name file after matrix, filler materials
f = open(fileName, 'w')
f.write('Matrix\tFiller\tPortion\tRadius\tNumber\tSide\tIntSize\tDelta\tCalcPor\tIntCond\tSeed\tNodes\tElements\tDevFac\tMeshSeed\tq\tdT\tk\tNoElmWarn\tWarn\n')
f.close()

for i in range(len(materials[matrix]['fillers'])): # "For each filler material"
	for j in range(len(materials[matrix]['fillers'][fillers[i]]['phr'])): # "For each PHR/volume portion specified"
		for k in range(trialsPer): # "Run trialsPer times for each material/filler/PHR combination"
			modelObject, modelName = createModel(2) # Create model database "Model-2"
			side, radius, portions, dP, dM, cP, cM = defExperiment(modelObject, matrix, fillers[i]) # Define material attributes for specified matrix, fillers
			seed = numpy.random.randint(1000000) # Random seed for coordinate generation
			radius, number = invPHR(portions[j], dP, radius, dM, side, False) # Returns specific radius size and number of inclusions for closest PHR value.
			
			defDelta = materials[matrix]['fillers'][fillers[i]]['delta'] # radius times this fraction is the minimum distance between particles
			xVals, yVals, zVals, warningPoints, number = getPoints3D(seed, side, radius, number, 0, defDelta) # returns coordinates for inclusion locations. 
			calcPHR = calculatePHR3D(number, dP, radius, dM, side)
			part = createMatrix(modelObject, side, False) # Create the matrix
			edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our matrix
			part2 = createSphereParticle(modelObject, radius, side) # Create Particle part
			edges2, vertices2, face2 = assignGeomSequence(part2) # Create references to important sets in particle 
			matrixSet, particleSet = create3DInitialSets(part, part2, side)
			createSection(modelObject, part, matrix, matrixSet) # Create section for matrix material
			createSection(modelObject, part2, fillers[i], particleSet) # Create section for filler material
			
			modelRootAssembly, fullMatrixPart = create3DMatrixInclusions(modelObject, part, part2, number, xVals, yVals, zVals) # Create assembly and return references to assembly sets
			assemblyTop, assemblyBottom, assemblyAll = define3DAssemblySets(modelRootAssembly, side)
			temp1, temp2 = 328.15, 298.15 # Assign heat temperature to be used in experiment
			heatStep3D(modelObject, assemblyBottom, assemblyTop, temp1, temp2) # apply heat BC
			
			meshSeed = materials[matrix]['fillers'][fillers[i]]['meshSeed'] # recommended mesh
			df = materials[matrix]['fillers'][fillers[i]]['df'] # recommended deviation factor
			elements, nodes, df, meshSeed = makeMesh3D(modelObject, modelRootAssembly) # Draw mesh and return number of nodes and elements
			
			makeElementSet(fullMatrixPart, modelRootAssembly)
			warningString, noElementsWarning = submitJob(modelName, fileName) # Submit job and take note of any warnings
			avgHF, TC = getThermalProperties3D(side, temp1, temp2, fileName) # Extract relevant information about thermal properties
			af = open(fileName, 'a')
			af.write(dataString(matrix, fillers[i], portions[j], radius, number, side, interfaceSize, defDelta, calcPHR, "nil", seed, nodes, elements, df, meshSeed, avgHF, temp1, temp2, TC, warningString, warningPoints, noElementsWarning)) # Write the data to file
			af.close()
			del mdb.jobs[fileName] # Optional to clear out the existing jobs and models each iteration.
			del mdb.models[modelName]
		
	

f.close()