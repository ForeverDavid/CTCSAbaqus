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
materials = getMaterialList() # Load in material Data
matrix = "ESBR" # Choose the first experiment from TCNanoFillers
fillers = ["ZincOxide", "Alumina"] # Two fillers associated with first experiment
totalIterations = 0


for i in range(len(materials[matrix]['fillers'])): # "For each filler material"
	for j in range(len(materials[matrix]['fillers'][fillers[i]]['phr'])): # "For each PHR/volume portion specified"
		for k in range(trialsPer): # "Run trialsPer times for each material/filler/PHR combination"
			combinations = len(materials[matrix]['fillers']) * len(materials[matrix]['fillers'][fillers[i]]['phr']) * trialsPer
			totalIterations = totalIterations + 1
			print("number " + str(totalIterations) + " of " + str(combinations))
			
			modelObject, modelName = createModel(2) # Create model database "Model-2"
			side, radius, portions, dP, dM, cP, cM = defExperiment(modelObject, matrix, fillers[i]) # Define material attributes for specified matrix, fillers
			
			seed = numpy.random.randint(1000000) # Random seed for coordinate generation
			radius, number = invPHR(portions[j], dP, radius, dM, side, False) # Returns specific radius size and number of inclusions for closest PHR value.
			
			delta = materials[matrix]['fillers'][fillers[i]]['delta']
			minInt = materials[matrix]['fillers'][fillers[i]]['delta']
			intPortionLimit = getInterfacePortionLimit(side, radius, number)
			interfacePortion = numpy.random.sample(1) * (intPortionLimit-minInt) + minInt # random minInt to limit inclusive
			interfacePortion = round(interfacePortion[0], 3)
			
			interfaceConductivity= numpy.random.sample(1) * (cP-cM) + cM # Between cM and cP
			interfaceConductivity = int(interfaceConductivity[0])
			
			xVals, yVals, zVals, warningPoints, number = getPoints3D(seed, side, radius, number, interfacePortion, delta) # returns coordinates for inclusion locations. 
			calcPHR = round(calculatePHR3D(number, dP, radius, dM, side))
			
			
			
			defineMaterial(modelObject, "Interface", dM, interfaceConductivity) # Define interface conductivity... Note that this will be generated randomly
			
			part = createMatrix(modelObject, side, False) # Create the matrix
			edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our matrix
			part2 = createSphereParticle(modelObject, radius, side) # Create Particle part
			edges2, vertices2, face2 = assignGeomSequence(part2) # Create references to important sets in particle 
			part3 = createSphereParticle(modelObject, (radius + radius*interfacePortion), side, "Interface") # Create interface part
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
			
			elements, nodes, df, meshSeed = makeMesh3D(modelObject, modelRootAssembly, meshSeed, df)  # Draw mesh and return number of nodes and elements
			makeElementSet(fullMatrixPart, modelRootAssembly)
			print(str(calcPHR) + " " str(portions[j]))
			print(str(seed) + " " +str(number) + " " + str(radius) + " " + str(df) + " " + str(meshSeed) + " " + str(elements) + " " +str(nodes) + " " + warningPoints)
			odbfileName = modelName
			warningString, noElementsWarning = submitJob(modelName, odbfileName)  # Submit job and take note of any warnings
			avgHF, TC = getThermalProperties3D(side, temp1, temp2, odbfileName) # Extract relevant information about thermal properties
			print(dataString(matrix, fillers[i], portions[j], radius, number, side, interfacePortion, delta, calcPHR, interfaceConductivity, seed, nodes, elements, df, meshSeed, avgHF, temp1, temp2, TC, warningString, warningPoints, noElementsWarning)) # Write the data to file
			#del mdb.jobs["Job-1"]
			#del mdb.models[modelName]
		
	

f.close()