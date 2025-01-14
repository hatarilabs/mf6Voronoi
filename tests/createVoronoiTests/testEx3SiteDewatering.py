from mf6Voronoi.geoVoronoi import createVoronoi

#Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='siteDewatering',maxRef = 10, multiplier=1.5)

#Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/siteDewatering/Shp/modelAoi2.shp')
vorMesh.addLayer('ghb','../../examples/siteDewatering/Shp/modelGhb.shp',2)
vorMesh.addLayer('drn','../../examples/siteDewatering/Shp/drains.shp',1)
vorMesh.addLayer('wel','../../examples/siteDewatering/Shp/pointWell.shp',1)

#Generate point pair array
vorMesh.generateOrgDistVertices()

#Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

#Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')