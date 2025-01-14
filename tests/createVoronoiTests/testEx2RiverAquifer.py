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