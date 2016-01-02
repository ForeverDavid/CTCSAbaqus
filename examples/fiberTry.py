# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(0.0, 225.0))
#* Rectangle cannot be created.
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(225.0, 225.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='matrix', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['matrix'].BaseSolidExtrude(depth=225.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    0.0, 0.0), point1=(0.0, 3.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='fiber', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['fiber'].BaseSolidExtrude(depth=20.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='matrix-1', 
    part=mdb.models['Model-1'].parts['matrix'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-1', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-2', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-3', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-4', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-5', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-6', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(167.992, 
    123.899, 62.3099), axisPoint=(0.0, 0.0, 0.0), instanceList=('fiber-1', ))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(132.045, 
    59.1809, 165.288), axisPoint=(0.0, 0.0, 0.0), instanceList=('fiber-2', ))
mdb.models['Model-1'].rootAssembly.deleteFeatures(('fiber-1', 'fiber-2', 
    'fiber-3', 'fiber-4', 'fiber-5', 'fiber-6'))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-1', part=
    mdb.models['Model-1'].parts['fiber'])
del mdb.models['Model-1'].rootAssembly.features['fiber-1']
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-1', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-1', ), 
    vector=(172.221, 117.0, 80.6))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(-4.229, 
    6.899, -18.2901), axisPoint=(172.221, 117.0, 80.6), instanceList=(
    'fiber-1', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-2', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-2', ), 
    vector=(119.213, 69.8651, 154.28))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(12.832, 
    -10.6842, 11.008), axisPoint=(119.213, 69.8651, 154.28), instanceList=(
    'fiber-2', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-3', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-3', ), 
    vector=(66.5381, 61.9916, 85.9758))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(2.8894, 
    15.3987, -12.4311), axisPoint=(66.5381, 61.9916, 85.9758), instanceList=(
    'fiber-3', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-4', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-4', ), 
    vector=(84.7198, 66.1184, 87.5074))
del mdb.models['Model-1'].rootAssembly.features['fiber-4']
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-4', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-4', ), 
    vector=(84.7198, 66.1184, 87.5074))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(-5.0578, 
    1.6187, -19.2821), axisPoint=(84.7198, 66.1184, 87.5074), instanceList=(
    'fiber-4', ))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='fiber-5', part=
    mdb.models['Model-1'].parts['fiber'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('fiber-5', ), 
    vector=(102.031, 155.939, 100.878))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(-18.7471, 
    -2.233, -6.6), axisPoint=(102.031, 155.939, 100.878), instanceList=(
    'fiber-5', ))
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
    mdb.models['Model-1'].rootAssembly.instances['matrix-1'], ), 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['matrix-1'], 
    name='solidmatrix', originalInstances=SUPPRESS)
mdb.models['Model-1'].rootAssembly.features['matrix-1'].resume()
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
    mdb.models['Model-1'].rootAssembly.instances['fiber-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-5']), instanceToBeCut=
    mdb.models['Model-1'].rootAssembly.instances['matrix-1'], name=
    'emptymatrix', originalInstances=SUPPRESS)
del mdb.models['Model-1'].rootAssembly.features['solidmatrix-1']
mdb.models['Model-1'].rootAssembly.resumeFeatures(('fiber-1', 'fiber-2', 
    'fiber-3', 'fiber-4', 'fiber-5'))
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['fiber-1'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-2'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-3'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-4'], 
    mdb.models['Model-1'].rootAssembly.instances['fiber-5'], 
    mdb.models['Model-1'].rootAssembly.instances['emptymatrix-1']), 
    keepIntersections=ON, name='fullmatrix', originalInstances=SUPPRESS)
