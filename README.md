mf6Voronoi
==========================
_The friendly way to create awesome Voronoi meshes for MODFLOW6 DISV_

<img src="https://raw.githubusercontent.com/hatarilabs/mf6Voronoi/refs/heads/main/examples/figures/voronoiMeshinModflow6Disv.png" alt="flopy3" style="width:50;height:20">


## Introduction
Groundwater modeling with several boundary conditions and complex hydrogeological setups require advanced tools for mesh discretizacion that ensures adequate refinement in the zone of interest while preserving a minimal cell account. Type of mesh has to be engineered in a way to preserve computational resources and represent adequately the groundwater flow regime. 

## Package
This python package creates a Voronoi mesh for MODFLOW6 with the DISV (discretized by vertices) option. The package work with geospatial files and has options for selective refinement based on the boundary condition.

These are the main python package characteristics:
- Works with geospatial files on ESRI Shapefile format
- Progressive refinement can be modified with a multiplier
- Summary of the point cloud generated for the Voronoi meshing
- Tested on more than 5 groundwater model datasets
- Output as polygon ESRI Shapefile
- Few steps and arguments for mesh generation

## Requirements
There are few requirements for the package. The most important one is that all the input files has to be in the same system of reference (CRS) and the CRS length unit has to be in meters or feet.

## Example

The package has been designed with a simple and user friendly approach allowing to create awesome meshes on a short amount of steps.


```
# Import the mf6Voronoi package
from mf6Voronoi.geoVoronoi import createVoronoi

# Create mesh object specifying the coarse mesh and the multiplier
vorMesh = createVoronoi(meshName='regionalModel',maxRef = 200, multiplier=1.5)

# Open limit layers and refinement definition layers
vorMesh.addLimit('basin','../../examples/regionalModel/shp/Angascancha_Basin_Extension.shp')
vorMesh.addLayer('river','../../examples/regionalModel/shp/rios.shp',50)

# Generate point pair array
vorMesh.generateOrgDistVertices()

# Generate the point cloud and voronoi
vorMesh.createPointCloud()
vorMesh.generateVoronoi()

# Export generated voronoi mesh
vorMesh.getVoronoiAsShp(outputPath='output')
```
