from mf6Voronoi.geoVoronoi import createVoronoi

#Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='riverAquifer',maxRef = 100, multiplier=1.5)

#Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/riverAquifer/Shp/ModelLimit1.shp')
vorMesh.addLayer('ghb','../../examples/riverAquifer/Shp/ModelGHB1.shp',20)
vorMesh.addLayer('riv','../../examples/riverAquifer/Shp/ModelRiver2.shp',10)
vorMesh.addLayer('wel','../../examples/riverAquifer/Shp/ModelWell2.shp',2)

#Generate point pair array
vorMesh.generateOrgDistVertices()

#Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

#Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')

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