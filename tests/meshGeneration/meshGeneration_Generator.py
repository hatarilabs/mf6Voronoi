from mf6Voronoi.geoVoronoi import createVoronoi
from mf6Voronoi.utils import getVoronoiAsShp
import matplotlib.pyplot as plt
import geopandas as gpd
import os, json

with open('meshGenerationData.json') as jsonFile:
    meshGenerationDict = json.load(jsonFile)
    
for meshName, meshDict in meshGenerationDict.items():
    #Create mesh object specifying the coarse mesh and the multiplier
    vorMesh = createVoronoi(meshName=meshName,maxRef = meshDict["maxRef"], multiplier=meshDict["multiplier"])

    datasetPath = os.path.join('../../examples/datasets',meshName,'shp')

    #Open limit layers and refinement definition layers
    vorMesh.addLimit(meshDict["limitLayer"]["limitName"], os.path.join(datasetPath,meshDict["limitLayer"]["limitShape"]))

    for layerList in  meshDict["layerLayer"]:
        vorMesh.addLayer(layerList[0],os.path.join(datasetPath,layerList[1]),layerList[2])

    #Generate point pair array
    vorMesh.generateOrgDistVertices()

    #Generate the point cloud and voronoi
    vorMesh.createPointCloud()
    vorMesh.generateVoronoi()

    # Export generated voronoi mesh
    shapePath='../checkFiles/meshGeneration/shp/'+meshName+'.shp'
    getVoronoiAsShp(vorMesh.modelDis, shapePath=shapePath)

    #check mesh generation
    from mf6Voronoi.meshProperties import meshShape
    import os

    # open the mesh file
    mesh=meshShape(shapePath)

    # get the list of vertices and cell2d data
    gridprops=mesh.get_gridprops_disv()

    cell2d = gridprops['cell2d']           #cellid, cell centroid xy, vertex number and vertex id list
    vertices = gridprops['vertices']       #vertex id and xy coordinates
    ncpl = gridprops['ncpl']               #number of cells per layer
    nvert = gridprops['nvert']             #number of verts
    centroids=gridprops['centroids']
        
    mesh.save_properties(os.path.join('../checkFiles/meshGeneration/json','%s_disvDict.json'%meshName))

    # Load the shapefile (replace 'your_shapefile.shp' with the path to your file)
    gdf = gpd.read_file(shapePath)

    # Create a plot
    fig, ax = plt.subplots(figsize=(10, 10))
    gdf.plot(figsize=(35,25), fc='crimson', alpha=0.3, ec='teal', ax=ax)

    # Remove axis for cleaner image (optional)
    ax.set_axis_off()

    # Save the figure as a PNG
    plt.savefig("../checkFiles/meshGeneration/png/%s.png"%meshName, bbox_inches="tight", dpi=300)

    plt.close()