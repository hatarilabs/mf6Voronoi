from mf6Voronoi.geoVoronoi import createVoronoi

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

#Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')