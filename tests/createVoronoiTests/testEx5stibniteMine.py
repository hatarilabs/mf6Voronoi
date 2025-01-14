from mf6Voronoi.geoVoronoi import createVoronoi

#Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='stibniteMine',maxRef = 500, multiplier=2.5)

#Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/stibniteMine/Shp/catchment_11N.shp')
vorMesh.addLayer('water','../../examples/stibniteMine/Shp/river_basin_11N.shp',20)
vorMesh.addLayer('pit','../../examples/stibniteMine/Shp/minePits_11N.shp',50)
vorMesh.addLayer('dump','../../examples/stibniteMine/Shp/mineDumps_11N.shp',50)
vorMesh.addLayer('faults','../../examples/stibniteMine/Shp/faults.shp',5)

#Generate point pair array
vorMesh.generateOrgDistVertices()

#Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

#Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')
