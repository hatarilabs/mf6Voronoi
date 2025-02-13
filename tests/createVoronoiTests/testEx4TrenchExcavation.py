from mf6Voronoi.geoVoronoi import createVoronoi

# from voronoi utils
from mf6Voronoi.utils import getVoronoiAsShp

#Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='trenchExcavation',maxRef = 50, multiplier=1.5)

#Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/trenchExcavation/Shp/modelAoi.shp')
vorMesh.addLayer('ghb','../../examples/trenchExcavation/Shp/compoundGhb.shp',20)
vorMesh.addLayer('wel','../../examples/trenchExcavation/Shp/pumpingWells.shp',5)
vorMesh.addLayer('drn','../../examples/trenchExcavation/Shp/trenchExcavationDissolved.shp',2)

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
    
mesh.save_properties(os.path.join(jsonPath,'disvDict.json'))