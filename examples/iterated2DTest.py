"""
	This is an example of how to run multiple iterations of a given experiment 
	and write the results to file. The idea is to load in the relevant material
	data, choose which experiment to simulate, then generate the matrix with 
	inclusions and subject it to heat boundary conditions. Since the inclusions
	are randomly generated then run the same simulation multiple times then the
	average should yield an accurate answer.
"""

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

trialsPer = 3 # Number of times to run each experiment

materials = getMaterialList() # Load in material Data
matrix = "ESBR" # Choose the first experiment from TCNanoFillers
fillers = ["Alumina", "ZincOxide"] # Two fillers associated with first experiment

fileName = matrix+fillers[0] + fillers[1] # Name file after matrix, filler materials
f = open(fileName, 'w')
f.write('Matrix\tFiller\tPortion\tRadius\tNumber\tSide\tSeed\tNodes\tq\tdT\tk\tWarn\n')

for i in range(len(materials[matrix]['fillers'])): # "For each filler material"
	modelObject, modelName = createModel(2) # Create model database "Model-1"
	fileName = modelName
	side, radius, portions, dP, dM, cP, cM = defExperiment(modelObject, matrix, fillers[i]) # Define material attributes for specified matrix, fillers
	
	for j in range(len(portions)): # "For each PHR/volume portion specified"
		for k in range(trialsPer): # "Run trialsPer times for each material/filler/PHR combination"
			seed = numpy.random.randint(10000) # Random seed for coordinate generation
			radius, number = invPHR(portions[j], dP, radius, dM, side) # Returns specific radius size and number of inclusions for closest PHR value.
			xVals, yVals, warningPoints, number = getPoints2D(seed, side, radius, number) # returns coordinates for inclusion locations. 
			part = createMatrix(modelObject, side) # Create the matrix
			edges1, vertices1, face1 = assignGeomSequence(part) # Create references to important sets in our part
			matrixSet, fillerSet, allSet = createCircleInclusion(modelObject, part, radius, number, xVals, yVals, side) # Draw inclusions in the matrix
			createSection(modelObject, part, matrix, matrixSet) # Create section for matrix material
			createSection(modelObject, part, fillers[i], fillerSet) # Create section for filler material
			assemblyTop, assemblyBottom = makeAssembly2D(modelObject, part) # Create assembly and return references to assembly sets
			temp1, temp2 = 328.15, 298.15 # Assign heat temperature to be used in experiment
			heatStep2D(modelObject, assemblyBottom, assemblyTop, temp1, temp2) # apply heat BC
			elements, nodes = makeMesh2D(modelObject, part, allSet) # Draw mesh and return number of nodes and elements
			warningString = submitJob(modelName,fileName) # Submit job and take note of any warnings
			avgHF, TC = getThermalProperties2D(side, temp1, temp2, fileName) # Extract relevant information about thermal properties
			f.write(dataString2D(matrix, fillers[i], portions[j], radius, number, side, seed, nodes, elements, avgHF, temp1, temp2, TC, warningString[0], warningPoints)) # Write the data to file
	
f.close()