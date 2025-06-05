from mf6Voronoi.geoVoronoi import createVoronoi
from mf6Voronoi.utils import getVoronoiAsShp
import os, json

with open('meshGenerationDict') as jsonFile:
    meshGenerationDict = json.load(jsonFile)

for meshName, meshDict in meshGenerationDict.values:
    #Create mesh object specifying the coarse mesh and the multiplier
    vorMesh = createVoronoi(meshName=meshName,maxRef = meshDict["maxRef"], multiplier=meshDict["multiplier"])

    datasetPath = os.path.join('../../examples/datasets',meshName)

    #Open limit layers and refinement definition layers
    vorMesh.addLimit(meshDict["limitLayer"]["limitName"], os.path.join(datasetPath,meshDict["limitLayer"]["limitShp"]))

    for layerList in  meshDict["layerLayer"]:
        vorMesh.addLayer(layerList[0],os.path.join(datasetPath,layerList[1]),layerList[2])

#Generate point pair array
vorMesh.generateOrgDistVertices()

#Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

# Export generated voronoi mesh
getVoronoiAsShp(vorMesh.modelDis, shapePath='output/'+vorMesh.modelDis['meshName']+'.shp')

#check mesh generation
from mf6Voronoi.meshProperties import meshShape
import os

# open the mesh file
mesh=meshShape('output/'+vorMesh.modelDis['meshName']+'.shp')

# get the list of vertices and cell2d data
gridprops=mesh.get_gridprops_disv()

cell2d = gridprops['cell2d']           #cellid, cell centroid xy, vertex number and vertex id list
vertices = gridprops['vertices']       #vertex id and xy coordinates
ncpl = gridprops['ncpl']               #number of cells per layer
nvert = gridprops['nvert']             #number of verts
centroids=gridprops['centroids']

#check or create an output folder
jsonPath = ('output/'+vorMesh.modelDis['meshName'])
if os.path.isdir(jsonPath):
    print('The output folder %s exists'%jsonPath)
else:
    os.mkdir(jsonPath)
    print('The output folder %s has been generated.'%jsonPath)
    
mesh.save_properties(os.path.join('../checkFiles/meshGeneration/json','%s_disvDict.json'%meshName))