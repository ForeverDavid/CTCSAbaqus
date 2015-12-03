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

## NOTE: That the following way of changing interfacePortions and conductivities
## isn't as robust as generating random values based on materialList maximum portions
## and defined conductivities. See inpExample for best practices.

### Increase distance delta in coordinate generation by interface length! 
interfacePortions = [0.15, 0.2, 0.25] # 0.4 could be too large depending on phr/volume material combination
interfaceConductivities = [500000, 2000000, 3500000, 5500000, 7500000, 9500000, 11000000]

trialsPer = 1 # Number of times to run each experiment

materials = getMaterialList() # Load in material Data
matrix = "ESBR" # Choose the first experiment from TCNanoFillers
fillers = ["Alumina", "ZincOxide"] # Two fillers associated with first experiment

fileName = matrix+fillers[0] + fillers[1] +"Interface" # Name file after matrix, filler materials
f = open(fileName, 'w')
f.write('Matrix\tFiller\tPortion\tRadius\tNumber\tSide\tIntSize\tDelta\tCalcPor\tIntCond\tSeed\tNodes\tElements\tDevFac\tMeshSeed\tq\tdT\tk\tNoElmWarn\tWarn\n')
f.close()

for i in range(len(materials[matrix]['fillers'])): # "For each filler material"
	for j in range(len(materials[matrix]['fillers'][fillers[i]]['phr'])): # "For each PHR/volume portion specified"
		for l in range(len(interfacePortions)):
			for k in range(trialsPer): # "Run trialsPer times for each material/filler/PHR combination"
				for n in range(len(interfaceConductivities)):
					modelObject, modelName = createModel(2) # Create model database "Model-2"
					side, radius, portions, dP, dM, cP, cM = defExperiment(modelObject, matrix, fillers[i]) # Define material attributes for specified matrix, fillers
					defineMaterial(modelObject, "Interface", dM, interfaceConductivities[n]) # Define interface conductivity... Note that this will be generated randomly
					
					seed = numpy.random.randint(1000000) # Random seed for coordinate generation
					radius, number = invPHR(portions[j], dP, radius, dM, side, False) # Returns specific radius size and number of inclusions for closest PHR value.
					xVals, yVals, zVals, warningPoints, number = getPoints3D(seed, side, radius, number, interfacePortions[l], 0.15) # returns coordinates for inclusion locations. 
					calcPHR = calculatePHR3D(number, dP, radius, dM, side)
					
					part = createMatrix(modelObject, side, False) # Create the matrix
					edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our matrix
					part2 = createSphereParticle(modelObject, radius, side) # Create Particle part
					edges2, vertices2, face2 = assignGeomSequence(part2) # Create references to important sets in particle 
					part3 = createSphereParticle(modelObject, (radius + radius*interfacePortions[l]), side, "Interface") # Create interface part
					edges3, vertices3, face3 = assignGeomSequence(part3) # Create references to important sets in interface
					matrixSet, particleSet, interfaceSet = create3DInitialSets(part, part2, side, part3) # create sets for particle, matrix, and interface
					
					createSection(modelObject, part, matrix, matrixSet) # Create section for matrix material
					createSection(modelObject, part2, fillers[i], particleSet) # Create section for filler material
					createSection(modelObject, part3, "Interface", interfaceSet) # Create section for interface
					
					modelRootAssembly, fullMatrixPart = create3DMatrixInclusions(modelObject, part, part2, number, xVals, yVals, zVals, part3) # Create assembly and return references to assembly sets
					assemblyTop, assemblyBottom, assemblyAll = define3DAssemblySets(modelRootAssembly, side)
					temp1, temp2 = 328.15, 298.15 # Assign heat temperature to be used in experiment
					heatStep3D(modelObject, assemblyBottom, assemblyTop, temp1, temp2) # apply heat BC
					
					meshSeed = materials[matrix]['fillers'][fillers[i]]['meshSeed'] # recommended mesh
					df = materials[matrix]['fillers'][fillers[i]]['df'] # recommended deviation factor
					elements, nodes, df, meshSeed = makeMesh3D(modelObject, modelRootAssembly) # Draw mesh and return number of nodes and elements
					makeElementSet(fullMatrixPart, modelRootAssembly)
					warningString, noElementsWarning = submitJob(modelName, fileName)  # Submit job and take note of any warnings
					avgHF, TC = getThermalProperties3D(side, temp1, temp2, fileName) # Extract relevant information about thermal properties
					
					af = open(fileName, 'a')
					af.write(dataString(matrix, fillers[i], portions[j], radius, number, side, interfacePortions[l], 0.15, calcPHR, interfaceConductivities[n], seed, nodes, elements, df, meshSeed, avgHF, temp1, temp2, TC, warningString, warningPoints, noElementsWarning)) # Write the data to file
					af.close()
					
					del mdb.jobs[fileName]
					del mdb.models[modelName]
		
	

f.close()