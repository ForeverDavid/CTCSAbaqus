"""
# Coordinates 
	* Add comments!
"""
from math import *
import numpy

"""
	 * Given a phr along with density, approximate radius size, and size of ESBR side
	 * returns a tuple containing two values: actual radius size and number of fillers
"""
def invPHR(phr, densityFiller, radiusFiller, densityMatrix, sideMatrix, twoD=True):
	import numpy
	
	radiusStep = 0.03*radiusFiller # Alter radius to accommodate varying phr 
	
	if twoD:
		numberFillers = numpy.arange(1.,36.) # 1 to 30 particles
		rsSquared = (numpy.arange(radiusFiller-10*radiusStep,  
			radiusFiller+10*radiusStep , radiusStep)) ** 2
		## Grid of number circles with variable radius
		nRGrid = numpy.outer(numberFillers,rsSquared) * pi 
		## Grid of differences
		diffGrid = sideMatrix ** 2 - nRGrid
	else:
		defaultN = 75
		endRange = calculatePHR3D(defaultN, densityFiller, radiusFiller+
			10*radiusStep, densityMatrix, sideMatrix)
		while phr < endRange:
			defaultN = int(defaultN * 0.1 + defaultN)
			endRange = calculatePHR3D(defaultN, densityFiller, radiusFiller+
				10*radiusStep, densityMatrix, sideMatrix)
		
		numberFillers = numpy.arange(1.,defaultN+1) # 1 to defaultN+1 particles
		rsCubed = (numpy.arange(radiusFiller-10*radiusStep,  
			radiusFiller+10*radiusStep , radiusStep)) ** 3
		## Grid of number circles with variable radius
		nRGrid = numpy.outer(numberFillers,rsCubed) * pi * 4.0/3.0
		## Grid of differences
		diffGrid = sideMatrix ** 3 - nRGrid 
	
	ratioGrid = nRGrid/diffGrid
	phrGrid = numpy.multiply((100 * densityFiller / densityMatrix), ratioGrid)
	# Matrix of distances from phr
	distn = numpy.power(numpy.power((phr-phrGrid), 2), 0.5)
	# Returns the combination yielding closest approximation to phr
	vals = numpy.nonzero(distn == distn.min()) 
	if twoD:
		r = numpy.sqrt(rsSquared[vals[1][0]]) 
	else:
		r = numpy.power((rsCubed[vals[1][0]]), 1/3.0)
	
	n = numberFillers[vals[0][0]] # number of inclusions
	return r, int(n)

# Need an input PHR and then an actual PHR same with volume.
# Add a uniform version! one that isn't so random for when i fix all 
# the inputs and just search over interface values
def getPoints2D(seed, side, radius, number):
	import random
	import numpy
	random.seed(seed)
	
	rng = side-2.2*radius
	randXs = (1.1* radius) + rng * numpy.random.rand(100000,1)
	randYs = (1.1* radius) + rng * numpy.random.rand(100000,1)
	
	delta = 2.1 * radius;
	xCoords = [randXs[0][0]]
	yCoords = [randYs[0][0]]
	numberCoords = 1
	
	for i in range(1, 100000):
		x = randXs[i]
		y = randYs[i]
		distances = numpy.sqrt(numpy.power((x-xCoords), 2) + numpy.power((y-yCoords), 2))
		mindist = numpy.min(distances)
		if numberCoords == number:
			break
		if (mindist > delta):
			xCoords.append(x[0])
			yCoords.append(y[0])
			numberCoords += 1
		
	
	warningMsg = ''
	if numberCoords != number:
		warningMsg = '*'
		number = numberCoords # Needed for updated return value.
	
	return xCoords, yCoords, warningMsg, number

# Be more descriptive. This is confusing
def getPoints3D(seed, side, radius, number, interfacePortion=0.0, deltaCoefficient=0.15):
	import random
	import numpy
	random.seed(seed)
	
	# Delta is distances between particles
	#delta = 2 * radius + (radius * deltaCoefficient) + (interfacePortion * radius) # Default delta will be 2 radius + one tenth of radius. s
	# Think this works.
	#delta = 2 * ((radius + (radius * interfacePortion)) + (deltaCoefficient*(radius + radius*interfacePortion)))
	
	# x, y inside matrix without touching sides
	# This may be an issue.
	#rngr = side - (delta + radius * deltaCoefficient + radius * interfacePortion)
	#rngr = side - delta
	#rngr = (delta / 2.0) + (delta/2.0) * 0.1
	#randXs = (radius + (radius * deltaCoefficient) + (interfacePortion * radius)) + rngr * numpy.random.rand(100000,1)
	#randYs = (radius + (radius * deltaCoefficient) + (interfacePortion * radius)) + rngr * numpy.random.rand(100000,1)
	#randZs = (radius + (radius * deltaCoefficient) + (interfacePortion * radius)) + rngr * numpy.random.rand(100000,1)
	
	#randXs = rngr + (side-(rngr)) * numpy.random.rand(100000,1)
	#randYs = rngr + (side-(rngr)) * numpy.random.rand(100000,1)
	#randZs = rngr + (side-(rngr)) * numpy.random.rand(100000,1)
	
	######
	r = radius
	i = radius * interfacePortion + radius
	dC = i * deltaCoefficient + i
	d = r + (i - r) + (dC - i)
	#d = radius + (interfacePortion * radius + radius)
	randXs = d + (side-(2 * d)) * numpy.random.rand(100000,1)
	randYs = d + (side-(2 * d)) * numpy.random.rand(100000,1)
	randZs = d + (side-(2 * d)) * numpy.random.rand(100000,1)
	
	xCoords = [randXs[0][0]]
	yCoords = [randYs[0][0]]
	zCoords = [randZs[0][0]]
	numberCoords = 1
	
	for i in range(1, 100000):
		x = randXs[i]
		y = randYs[i]
		z = randZs[i]
		distances = numpy.sqrt(numpy.power((x-xCoords), 2) + numpy.power((y-yCoords), 2) + numpy.power((z-zCoords), 2))
		mindist = numpy.min(distances)
		if numberCoords == number:
			break
		if (mindist > 2 * d):
			xCoords.append(x[0])
			yCoords.append(y[0])
			zCoords.append(z[0])
			numberCoords += 1
		
	
	warningMsg = ''
	if numberCoords != number:
		warningMsg = '?'
		number = numberCoords # Needed for updated return value.
	
	return xCoords, yCoords, zCoords, warningMsg, int(round(number))

# Alternate version of other invPHR function
# NOTE: consider having one that returns closes approximation for each number size
# so we can make sure that number of particles doesn't disproportionately influence 
# TC. 
# NOTE: This is a problem because sometimes we should have more of an error
# but be closer to the original radius then dramatically altering to fit. 
def invPHRAlternate3D(phr, densityFiller, radiusFiller, densityMatrix, sideMatrix):
	import numpy
	radiusStep = 0.025*radiusFiller # Alter radius to accommodate varying phr 
	stepSize = 15
	nS = [1,8,27,64] 
	endRange = calculatePHR3D(nS[3], densityFiller, 
		radiusFiller+stepSize*radiusStep, densityMatrix, sideMatrix)
	lowRange = calculatePHR3D(nS[0], densityFiller, 
		radiusFiller-stepSize*radiusStep, densityMatrix, sideMatrix)
	
	## NOTE SOMETHING WRONG HHERE THATT TAKES LONG TIME
	while phr > endRange and phr < lowRange:
		if phr > endRange:
			stepSize = stepSize + 3
			endRange = calculatePHR3D(nS[3], densityFiller, radiusFiller+
				stepSize*radiusStep, densityMatrix, sideMatrix)
		else:
			stepSize = stepSize + 3
			lowRange = calculatePHR3D(nS[0], densityFiller, radiusFiller+
				stepSize*radiusStep, densityMatrix, sideMatrix)
			
	
	rsCubed = (numpy.arange(radiusFiller-stepSize*radiusStep,  
		radiusFiller+stepSize*radiusStep , radiusStep)) ** 3
	## Grid of number circles with variable radius
	nRGrid = numpy.outer(nS,rsCubed) * pi * 4.0/3.0
	## Grid of differences
	diffGrid = sideMatrix ** 3 - nRGrid 
	
	ratioGrid = nRGrid/diffGrid
	phrGrid = numpy.multiply((100 * densityFiller / densityMatrix), ratioGrid)
	# Matrix of distances from phr
	distn = numpy.power(numpy.power((phr-phrGrid), 2), 0.5)
	# Returns the combination yielding closest approximation to phr
	vals = numpy.nonzero(distn == distn.min()) 
	r = numpy.power((rsCubed[vals[1][0]]), 1/3.0)
	n = nS[vals[0][0]] # number of inclusions
	return r, int(round(n))

def getInterfacePortionLimit(side, radius, number, delta=0.15):
	numerator = side - (2.0 * (number ** (1/3.0))) * (radius * delta + radius)
	denominator = 2.0 * (number ** (1/3.0)) * (radius * delta + radius)
	intPortionLimit = numerator / float(denominator)
	return intPortionLimit

def getPoints3dDeterministic(side, radius, number):
	import numpy
	xs = []
	ys = []
	zs = []
	
	impNum = (2 *(number ** (1/3.0)))
	impNum = int(round(impNum))
	nz = [y+1 for y in range(impNum) if (y+1) % 2 == 1]
	sizeInc = side / float(impNum)
	
	for i in nz:
		for j in nz:
			for k in nz:
				xs.append(sizeInc*i)
				ys.append(sizeInc*j)
				zs.append(sizeInc*k)
	
	return xs, ys, zs

## TODO:
# LIMIT THE STEP INITIALLY AND GROW ONLY IF NEEDED
# Need to ensure inputs don't exceed 1 because this gives greater than 1 volumePortions
## This is the deterministic point generation method
def invVolumeAlternate3D(volPortion, radiusFiller, sideMatrix):
	import numpy
	radiusStep = 0.025*radiusFiller # Alter radius to accommodate varying volPortion 
	stepSize = 15
	nS = [1,8,27,64] 
	endRange = calculateVolume(nS[3], radiusFiller+stepSize*radiusStep, sideMatrix)
	lowRange = calculateVolume(nS[0], radiusFiller-stepSize*radiusStep, sideMatrix)
	
	## NOTE SOMETHING WRONG HHERE THATT TAKES LONG TIME
	while volPortion > endRange and volPortion < lowRange:
		if volPortion > endRange:
			stepSize = stepSize + 3
			endRange = calculateVolume(nS[3], radiusFiller+stepSize*radiusStep, sideMatrix)
		else:
			stepSize = stepSize + 3
			lowRange = calculateVolume(nS[0], radiusFiller+stepSize*radiusStep, sideMatrix)
		
	
	rsCubed = (numpy.arange(radiusFiller-stepSize*radiusStep,  
		radiusFiller+stepSize*radiusStep , radiusStep)) ** 3
	## Grid of number circles with variable radius
	nRGrid = numpy.outer(nS,rsCubed) * pi * 4.0/3.0
	## Grid of differences
	matrixGrid = sideMatrix ** 3
	ratioGrid = nRGrid/matrixGrid
	# Matrix of distances from volPortion
	distn = numpy.power(numpy.power((volPortion-ratioGrid), 2), 0.5)
	# Returns the combination yielding closest approximation to volPortion
	vals = numpy.nonzero(distn == distn.min()) 
	r = numpy.power((rsCubed[vals[1][0]]), 1/3.0)
	n = nS[vals[0][0]] # number
	return r, n

# AM i doing this all wrong? Perhaps I could write a method to
# take side, take number of particles, delta...

""" Need to fix!
"""
"""
	Returns number of fillers
	* ToDo : Add smaller particles instead of rounding
	* Specific to Spheres, need functionality for squares, ellipses, and TMOICF
"""
def invVolumeSphere(radiusFiller, sideMatrix, volumePortion):
	import numpy
	numFillers = round((3.0 / 4.0 * volumePortion * sideMatrix ** 3) / (pi * radiusFiller ** 3))
	return numFillers

def calculatePHR3D(n, densityFiller, radiusFiller, densityMatrix, sideMatrix):
	numerator = n * pi * radiusFiller * radiusFiller * radiusFiller * (4/3.0)
	denominator = sideMatrix * sideMatrix * sideMatrix - numerator
	phr = (numerator / denominator) * (densityFiller / densityMatrix) * 100
	return phr

def calculateVolume(n, radiusFiller, sideMatrix):
	return ((n * pi * (4.0/3.0) * radiusFiller ** 3) / (sideMatrix ** 3))

# Need a better way for dealing with fractional inclusions similar to invPHR function
def invAreaCircle(radiusFiller, sideMatrix, areaRatio):
	import numpy
	numberFillers = round((areaRatio * sideMatrix ** 2) / (pi * radiusFiller ** 2))
	return numberFillers

def calcNumberFibers(radiusFiber, fiberLength, sideMatrix, volumePortion):
	import numpy
	numFibers = round((volumePortion * sideMatrix ** 3) / (fiberLength * radiusFiber ** 2 * pi))
	return int(numFibers)

def calcFiberPortion(n, radiusFiber, fiberLength, sideMatrix):
	import numpy
	return ((n * pi * radiusFiber ** 2 * fiberLength) / (sideMatrix**3))

# Add methods for fibers
def invAreaRectangle(lengthFiller, heightFiller, sideMatrix, areaRatio):
	import numpy
	numberFillers = round((areaRatio * sideMatrix ** 2) / (lengthFiller * heightFiller))
	return numberFillers

# Random orientation
def getPointsRectangle2D(seed, side, length, width, number):
	import random
	import numpy
	random.seed(seed)
	
	
	
	rng = int(side-2.2*radius)
	randXs = (1.1* radius) + rng * numpy.random.rand(100000,1)
	randYs = (1.1* radius) + rng * numpy.random.rand(100000,1)
	
	delta = 2.1 * radius;
	xCoords = [randXs[0]]
	yCoords = [randYs[0]]
	numberCoords = 1
	
	for i in range(1, 100000):
		x = randXs[i]
		y = randYs[i]
		distances = numpy.sqrt(numpy.power((x-xCoords), 2) + numpy.power((y-yCoords), 2))
		mindist = numpy.min(distances)
		if numberCoords == number:
			break
		if (mindist > delta):
			xCoords.append(x)
			yCoords.append(y)
			numberCoords += 1
		
	
	warningMsg = ''
	if numberCoords != number:
		warningMsg = '*'
	
	return xCoords, yCoords, warningMsg



def get3DCylinders(seed, side, radius, height, number):
	import random
	import numpy
	random.seed(seed)
	
	cylCoords = [Cylinder(getC(radius, height, side), getW(), radius, height)]
	numberCoords = 1
	
	for i in range(1, 100000):
		nextCylinder = Cylinder(getC(radius, height, side), getW(), radius, height)
		flag = True
		for j in range(0, len(cylCoords)):
			if not SeperatedCylinders(cylCoords[j], nextCylinder):
				flag = False
				break
		
		if flag:
			cylCoords.append(nextCylinder)
			numberCoords += 1
		if numberCoords == number:
			break
	
	
	warningMsg = ''
	if numberCoords != number:
		warningMsg = '?'
		number = numberCoords # Needed for updated return value.
	
	return cylCoords, warningMsg, int(round(number))


# c is center point of cylnder
# unit length axis direction w which is computed randomly
# radius r
# height h
# end disks are centered c +- h/2 W ...
# u, v
# parametrized = X(theta, t) = C + s cos theta U + s sin theta V + tW , 0 <= theta < 2pi , 0 <= s <= r , |t| < h/2
class Cylinder:
	def __init__(self, centerPoint, unitLengthAxisDirectionVec, radius, height):
		self.c = centerPoint
		self.w = unitLengthAxisDirectionVec
		self.r = radius
		self.h = height
	
	def getEnds(self):
		return [self.c + self.w*(self.h/2.0), self.c - self.w*(self.h/2.0)]

def getW():
	import numpy
	z = numpy.random.uniform(-1, 1, 1)[0]
	theta = numpy.random.uniform(0, 2*pi, 1)[0]
	return numpy.array([sqrt(1-z*z)*cos(theta), sqrt(1-z*z)*sin(theta), z])

def cylindricaldegree(cylinders):
	for i in range(0, len(cylinders)):
		a = cylinders[i].getEnds()[0]
		b = cylinders[i].getEnds()[1]
		degreeAB = math.degrees(math.acos(numpy.dot(a,b)/(numpy.linalg.norm(a)*numpy.linalg.norm(b))))
		print(degreeAB)

def cylindricaldegreeother(cylinders):
	for i in range(0, len(cylinders)):
		a = numpy.array(cylinders[i][0][0])
		b = numpy.array(cylinders[i][0][1])
		degreeAB = math.degrees(math.acos(numpy.dot(a,b)/(numpy.linalg.norm(a)*numpy.linalg.norm(b))))
		print(degreeAB)

def getC(r0, h0, side):
	import numpy
	x0 = numpy.random.uniform(r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
	y0 = numpy.random.uniform(r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
	z0 = numpy.random.uniform(r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
	return numpy.array([x0, y0, z0])

def SeperatedCylinders(cylinder1, cylinder2):
	import numpy
	w1 = cylinder1.w
	w2 = cylinder2.w
	delta = numpy.subtract(cylinder2.c, cylinder1.c)
	w1Xw2 = numpy.cross(w1, w2)
	lenw1Xw2 = numpy.linalg.norm(w1Xw2)
	rSum = cylinder1.r + cylinder2.r
	h1Div2 = cylinder1.h / 2.0
	h2Div2 = cylinder2.h / 2.0
	r1 = cylinder1.r
	r2 = cylinder2.r
	
	
	if lenw1Xw2 > 0:
		if r2*lenw1Xw2 + h1Div2 + h2Div2 * numpy.linalg.norm(numpy.dot(w1, w2)) - numpy.linalg.norm(numpy.dot(w1, delta)) < 0:
			return True
		if r1*lenw1Xw2 + h1Div2 * numpy.linalg.norm(numpy.dot(w1, w2)) + h2Div2  - numpy.linalg.norm(numpy.dot(w2, delta)) < 0:
			return True
		if rSum*lenw1Xw2 - numpy.linalg.norm(numpy.dot(w1Xw2, delta)) < 0:
			return True
		if SeperatedByCylinderPerpendiculars(cylinder1, cylinder2):
			return True
		if SeperatedByCylinderPerpendiculars(cylinder2, cylinder1):
			return True
		#if SeperatedByOtherDirections(cylinder1, cylinder2, delta):
		#	return True
	else:
		if h1Div2 + h2Div2 - numpy.linalg.norm(w1, delta) < 0:
			return True
		if rSum - numpy.linalg.norm(delta - numpy.dot(w1, delta)*w1) < 0:
			return True
	
	return False

def F(t, r0, r1, h1, b1, c1, a2, b2):
	import numpy
	omt = 1 - t
	tsqr = t * t
	c1sqr = c1 * c1
	omtsqr = omt * omt
	h1b1Div2 = h1 * b1 / 2.0
	term0 = r0 * sqrt(omtsqr + tsqr)
	term1 = r1 * sqrt(omtsqr + c1sqr * tsqr)
	term2 = h1b1Div2 * t
	term3 = numpy.linalg.norm(omt * a2 + t * b2)
	return term0 + term1 + term2 - term3

def FDer(t, r0, r1, h1, b1, c1, a2, b2):
	import numpy
	omt = 1 - t
	tsqr = t * t
	c1sqr = c1 * c1
	omtsqr = omt * omt
	h1b1Div2 = h1 * b1 / 2.0
	term0 = r0 * (2 * t - 1) / sqrt(omtsqr + tsqr)
	term1 = r1 *((1 + c1sqr) * t - 1) / sqrt(omtsqr + c1sqr * tsqr)
	term2 = h1b1Div2
	term3 = (b2 - a2) * numpy.sign(omt * a2 + t * b2)
	return term0 + term1 + term2 - term3

def SeperatedByCylinderPerpendiculars(cylinder1, cylinder2):
	import numpy
	w1 = cylinder1.w
	w2 = cylinder2.w
	delta = numpy.subtract(cylinder2.c, cylinder1.c)
	c1 = numpy.dot(w1, w2)
	b1 = sqrt(1-c1*c1)
	v0 = (w2 - c1*w1)/b1
	u0 = numpy.cross(v0, w1)
	a2 = numpy.dot(delta, u0)
	b2 = numpy.dot(delta, v0)
	r1 = cylinder1.r
	r2 = cylinder2.r
	h1 = cylinder1.h
	h2 = cylinder2.h
	
	if F(0, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return True
	if F(1, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return True
	if FDer(0, r1, r2, h2, b1, c1, a2, b2) >= 0:
		return False
	if FDer(1, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return False
	
	t0 = 0
	t1 = 1
	maxIterations = 25 # I picked this arbitrarily
	
	for i in range (0, maxIterations):
		tmid = 0.5 * (t0 + t1)
		if F(tmid, r1, r2, h2, b1, c1, a2, b2) <= 0:
			return True
		fdmid = FDer(tmid, r1, r2, h2, b1, c1, a2, b2)
		if (fdmid > 0):
			t1 = tmid
		elif fdmid < 0:
			t0 = tmid
		else:
			break
	
	a2 = -1 * a2
	if F(0, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return True
	if F(1, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return True
	if FDer(0, r1, r2, h2, b1, c1, a2, b2) >= 0:
		return False
	if FDer(1, r1, r2, h2, b1, c1, a2, b2) <= 0:
		return False
	
	t0 = 0
	t1 = 1
	for i in range (0, maxIterations):
		tmid = 0.5 * (t0 + t1)
		if F(tmid, r1, r2, h2, b1, c1, a2, b2) <= 0:
			return True
		fdmid = FDer(tmid, r1, r2, h2, b1, c1, a2, b2)
		if (fdmid > 0):
			t1 = tmid
		elif fdmid < 0:
			t0 = tmid
		else:
			break
	
	return False

def G(s, t, r0, h0, r1, h1, a0, b0, c0, a1, b1, c1, delta):
	import numpy
	lenDelta = numpy.norm(delta)
	h0Div2 = h0 / 2.0
	h1Div2 = h1 / 2.0
	omsmt = 1 - s - t
	ssqr = s * s
	tsqr = t * t
	omsmtsqr = omsmt * omsmt
	temp = ssqr + tsqr + omsmtsqr
	L0 = a0 * s + b0 * t + c0 * omsmt
	L1 = a1 * s + b1 * t + c1 * omsmt
	Q0 = temp - L0 * L0
	Q1 = temp - L1 * L1
	return r0 * sqrt(Q0) + r1 * sqrt(Q1) + h0Div2 * numpy.norm(L0) + h1Div2 * numpy.norm(L1) - omsmt * lenDelta

def GDer(s, t, r0, h0, r1, h1, a0, b0, c0, a1, b1, c1, delta):
	import numpy
	lenDelta = numpy.norm(delta)
	h0Div2 = h0 / 2.0
	h1Div2 = h1 / 2.0
	omsmt = 1 - s - t
	ssqr = s * s
	tsqr = t * t
	omsmtsqr = omsmt * omsmt
	temp = ssqr + tsqr + omsmtsqr
	L0 = a0 * s + b0 * t + c0 * omsmt
	L1 = a1 * s + b1 * t + c1 * omsmt
	Q0 = temp - L0 * L0
	Q1 = temp - L1 * L1
	diffS = s - omsmt
	diffT = t - omsmt
	diffa0c0 = a0 - c0
	diffa1c1 = a1 - c1
	diffb0c0 = b0 - c0
	diffb1c1 = b1 - c1
	halfQ0s = diffS - diffa0c0 * L0
	halfQ1s = diffS - diffa1c1 * L1
	halfQ0t = diffT - diffb0c0 * L0
	halfQ1t = diffT - diffb1c1 * L1
	factor0 = r0 / sqrt(Q0)
	factor1 = r1 / sqrt(Q1)
	signL0 = numpy.sign(L0)
	signL1 = numpy.sign(L1)
	
	gradient = numpy.array([0,0])
	gradient[0] += halfQ0s * factor0
	gradient[0] += halfQ1s * factor1
	gradient[0] += h0Div2 * diffa0c0 * signL0
	gradient[0] += h1Div2 * diffa1c1 * signL1
	gradient[0] += lenDelta
	gradient[1] += halfQ0t * factor0
	gradient[1] += halfQ1t * factor1
	gradient[1] += h0Div2 * diffb0c0*signL0
	gradient[1] += h1Div2 * diffb1c1*signL1
	gradient[1] += lenDelta
	
	return gradient

# points parametrized by X(theta, t) = C + s cos theta U + s sin theta V + t W, 0 < theta < 2pi, 0 <= s <= r, |t| < h/2
# so imagine sphere just inside within a margin, so like  h/10.0 < |t| < h/2 - h/10 and r/5.0 <= s <= 4r/5.0 
# this is excessive but should work fairly well...
#import numpy
#	w1 = cylinder1.w
#	w2 = cylinder2.w
#	delta = numpy.subtract(cylinder2.c, cylinder1.c)
#	w1Xw2 = numpy.cross(w1, w2)
#	lenw1Xw2 = numpy.linalg.norm(w1Xw2)
#	rSum = cylinder1.r + cylinder2.r
#	h1Div2 = cylinder1.h / 2.0
#	h2Div2 = cylinder2.h / 2.0
#	r1 = cylinder1.r
#	r2 = cylinder2.r

def getRangeInteriorCylinderLine(cylinder1, stepSize):
	# define the range of x,y,z
	import numpy
	cyl1Ends = cylinder1.getEnds()
	if cyl1Ends[0][0] < cyl1Ends[1][0]:
		w0 = numpy.abs(cylinder1.w[0])
	else:
		w0 = numpy.abs(cylinder1.w[0]) * -1

	if cyl1Ends[0][1] < cyl1Ends[1][1]:
		w1 = numpy.abs(cylinder1.w[1])
	else:
		w1 = numpy.abs(cylinder1.w[1]) * -1

	if cyl1Ends[0][2] < cyl1Ends[1][2]:
		w2 = numpy.abs(cylinder1.w[2])
	else:
		w2 = numpy.abs(cylinder1.w[2]) * -1

	x_range = numpy.linspace(cyl1Ends[0][0] + cylinder1.r * w0, cyl1Ends[1][0] - cylinder1.r * w0, stepSize)
	y_range = numpy.linspace(cyl1Ends[0][1] + cylinder1.r * w1, cyl1Ends[1][1] - cylinder1.r * w1, stepSize)
	z_range = numpy.linspace(cyl1Ends[0][2] + cylinder1.r * w2, cyl1Ends[1][2] - cylinder1.r * w2, stepSize)
	return x_range, y_range, z_range

# could we use function pointers?

# Higher step size gives higher likelihood of being correct but increases the 
# run time considerably. n * 3 * n, 3n^2 + ... + ... at worst
def BruteNonIntersectingCylinders(cylinder1, cylinder2, stepSize): 
	x_range1, y_range1, z_range1 = getRangeInteriorCylinderLine(cylinder1, stepSize)
	x_range2, y_range2, z_range2 = getRangeInteriorCylinderLine(cylinder2, stepSize*3)
	for i in range(0, stepSize):
		firstSphereCyl1 = numpy.array([x_range1[i], y_range1[i], z_range1[i]])
		for j in range(0, stepSize*3):
			sphereCyl2 = numpy.array([x_range2[j], y_range2[j], z_range2[j]])
			dist = numpy.linalg.norm(firstSphereCyl1-sphereCyl2)
			if dist <= cylinder1.r:
				return False
			
		
	return True

def get3DCylindersBrute(seed, side, radius, height, number):
	import random
	import numpy
	random.seed(seed)
	
	cylCoords = [Cylinder(getC(radius, height, side), getW(), radius, height)]
	numberCoords = 1
	
	for i in range(1, 100000):
		nextCylinder = Cylinder(getC(radius, height, side), getW(), radius, height)
		for j in range(0, len(cylCoords)):
			if not BruteNonIntersectingCylinders(cylCoords[j], nextCylinder, 7):
				break
		
		cylCoords.append(nextCylinder)
		numberCoords += 1
		
		if numberCoords == number:
			break
	
	
	warningMsg = ''
	if numberCoords != number:
		warningMsg = '?'
		number = numberCoords # Needed for updated return value.
	
	return cylCoords, warningMsg, int(round(number))

#def SeparatedByOtherDirections(cylinder1, cylinder2, delta):
	#if G(...) <= 0:
	#	return True
	
	#return False



class RandomCylinder(Cylinder):
	def __init__(self, r0, h0, side):
		Cylinder.__init__(self, RandomCylinder.getC(r0, h0, side), RandomCylinder.getW(), r0, h0)
	
	def getC(r0, h0, side):
		import numpy
		x0 = numpy.random.uniform(0 + r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
		y0 = numpy.random.uniform(0 + r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
		z0 = numpy.random.uniform(0 + r0 + h0/2.0, side-(r0+h0/2.0), 1)[0]
		return Point(x0, y0, z0)
	
	def getW():
		import numpy
		z = numpy.random.uniform(-1, 1, 1)[0]
		theta = numpy.random.uniform(0, 2*pi, 1)[0]
		return Vector(sqrt(1-z*z)*cos(theta), sqrt(1-z*z)*sin(theta), z)
	



"""

# Add a method that gives warning for possible combinations that aren't feasible.

# All unused! 

def invAreaSquare(fillerSide, siliconSide, volumePortion):

	Returns number of fillers

	import numpy
	volumePortion = volumePortion / 100.0
	numFillers = round((volumePortion * siliconSide * siliconSide) / (volumePortion * fillerSide * fillerSide + fillerSide * fillerSide))
	return numFillers
 


def invAreaEllipse(fillerLength, fillerHeight, siliconSide, volumePortion)
	import numpy
	volumePortion = volumePortion / 100.0
	numFillers = ((volumePortion * siliconSide * siliconSide) / (volumePortion * .5 * fillerLength * .5 * fillerHeight * pi + .5 * fillerLength * .5 * fillerHeight * pi))
	return numFillers

def invVolumeCube(fillerSide, siliconSide, volumePortion):
	
	import numpy
	volumePortion = volumePortion / 100.0
	numFillers = round((volumePortion * (siliconSide ** 3)) / (volumePortion * (fillerSide ** 3) + (fillerSide ** 3)))
	return numFillers
	
def invVolumeBlock(fillerLength, fillerHeight, fillerWidth, siliconSide, volumePortion):
	
	import numpy
	volumePortion = volumePortion / 100.0
	numFillers = round((volumePortion * (siliconSide ** 3)) / (volumePortion * fillerLength * fillerHeight * fillerWidth + fillerLength * fillerHeight * fillerWidth))
	return numFillers

def invVolumeElipsoid(fillerRLength, fillerRHeight, fillerRWidth, siliconSide, volumePortion):
	
	import numpy
	volumePortion = volumePortion / 100.0
	numFillers = round((volumePortion * (siliconSide ** 3)) / (volumePortion * (4/3.0) * pi * fillerRHeight * fillerRLength * fillerRWidth + fillerRLength * fillerRHeight * fillerRWidth * (4/3.0) * pi))
	return numFillers
"""