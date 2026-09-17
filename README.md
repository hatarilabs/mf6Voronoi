# mf6Voronoi

_The friendly way to create awesome Voronoi meshes for MODFLOW6 DISV_

<img src="https://raw.githubusercontent.com/hatarilabs/mf6Voronoi/refs/heads/main/examples/figures/voronoiMeshinModflow6Disv.png" alt="flopy3" style="width:50;height:20">

## Introduction

Groundwater modeling with several boundary conditions and complex hydrogeological setups require advanced tools for mesh discretizacion that ensures adequate refinement in the zone of interest while preserving a minimal cell account. Type of mesh has to be engineered in a way to preserve computational resources and represent adequately the groundwater flow regime.

## Package

This python package creates a Voronoi mesh for MODFLOW6 with the DISV (discretized by vertices) option. The package work with geospatial files and has options for selective refinement based on the boundary condition.

These are the main python package characteristics:

- Works with geospatial files on ESRI Shapefile format and GeoJSON
- Progressive refinement can be modified with a multiplier
- Summary of the point cloud generated for the Voronoi meshing
- Tested on more than 5 groundwater model datasets
- Output as polygon ESRI Shapefile
- Few steps and arguments for mesh generation

## Key Features & Recent Updates

- **GeoPandas Integration:** Complete transition from `fiona` to `geopandas` for faster, vector-optimized spatial data processing.
- **Native GeoJSON Support:** Direct ingestion of GeoJSON datasets alongside traditional Shapefiles.
- **Parallel Processing with Dask:** Scalable mesh generation and spatial computations designed for regional-scale river basins.
- **Vector Flow Visualization:** Advanced tools for extracting and plotting flow vectors in both plan view (2D) and vertical cross-sections.

For a detailed list of changes and release history, see the [Changelog](CHANGELOG.md).

## Benchmarks & Testing

- **Spatial Compatibility:** Validated against 20 distinct spatial datasets.
- **Scalability:** Dask parallelization stress-tested across 6 regional river basins.
- **End-to-End Stability:** Full MODFLOW 6 / Voronoi pipeline verified on 5 complete model scenarios.
- **Flow Vector Tools:** Plan-view and cross-sectional vector generation tested across 5 production cases.

## Requirements

There are few requirements for the package. The most important one is that all the input files has to be in the same system of reference (CRS) and the CRS length unit has to be in meters or feet.

## Tutorials & Learning Resources

For hands-on guides covering mesh generation, GeoJSON integration, Dask parallel processing, and flow vector visualization, visit the official tutorial series:
Complete Guide Collection: [Hatari Labs - Comprehensive MODFLOW 6 & mf6Voronoi Tutorials](https://hatarilabs.com/ih-en/comprehensive-modflow6-and-mf6voronoi-tutorials-collection)

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
