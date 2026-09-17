# Changelog

All notable changes to the `mf6Voronoi` project will be documented in this file.

## [0.0.38] - 2026-09-17

### Added

- **Dask Parallelization:** Native support for parallel computation on regional-scale basins.
- **GeoJSON Support:** Direct ingestion of GeoJSON format for spatial inputs.
- **Flow Vector Visualization:** Extracted cross-sectional and plan-view flow vector generators (`graphs2d_6`).

### Changed

- **GIS Engine Migration:** Replaced `fiona` dependency with `geopandas` for improved performance and spatial operations.
- **Enhanced Testing Suite:** Added full test coverage across 20 spatial datasets, 6 Dask regional basins, 5 complete MODFLOW 6 cases, and 5 vector visualization cases.

### Fixed

- Fixed `ZeroDivisionError` in vertical cross-sectional flow vector generation.
- Fixed visibility issues in quiver plots using `scale_units='width'`.
