from mf6Voronoi.geoVoronoi import createVoronoi

#Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='regionalModel',maxRef = 200, multiplier=1.5)

#Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/regionalModel/shp/Angascancha_Basin_Extension.shp')
vorMesh.addLayer('river','../../examples/regionalModel/shp/rios.shp',50)

#Generate point pair array
vorMesh.generateOrgDistVertices()

#Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

#Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')