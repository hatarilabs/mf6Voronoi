import os
import sys
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon, LineString, Point
import pytest

# Add repo root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from mf6Voronoi.voronoi import VoronoiGrid


def create_synthetic_data(tmp_path):
    # 1. Boundary (1000x1000 box)
    poly = Polygon([(0, 0), (1000, 0), (1000, 1000), (0, 1000), (0, 0)])
    gdf_bound = gpd.GeoDataFrame({"geometry": [poly]}, crs="EPSG:32618")
    bound_path = tmp_path / "boundary.shp"
    gdf_bound.to_file(bound_path)

    # 2. Refinement Area (Polygon)
    poly_ref = Polygon([(100, 100), (300, 100), (300, 300), (100, 300), (100, 100)])
    gdf_ref = gpd.GeoDataFrame({"geometry": [poly_ref]}, crs="EPSG:32618")
    ref_path = tmp_path / "refinement.shp"
    gdf_ref.to_file(ref_path)

    # 3. Refinement Line
    line = LineString([(500, 0), (500, 1000)])
    gdf_line = gpd.GeoDataFrame({"geometry": [line]}, crs="EPSG:32618")
    line_path = tmp_path / "river.shp"
    gdf_line.to_file(line_path)

    # 4. Fixed Points (Wells)
    p1 = Point(800, 800)
    p2 = Point(850, 850)
    gdf_pts = gpd.GeoDataFrame({"geometry": [p1, p2]}, crs="EPSG:32618")
    pts_path = tmp_path / "wells.shp"
    gdf_pts.to_file(pts_path)

    # 5. DEM
    dem_path = tmp_path / "dem.tif"
    res = 10
    rows, cols = 100, 100
    transform = from_origin(0, 1000, res, res)
    data = np.zeros((rows, cols), dtype=np.float32)
    # create a slope
    for r in range(rows):
        for c in range(cols):
            data[r, c] = r  # simple ramp

    with rasterio.open(
        dem_path,
        "w",
        driver="GTiff",
        height=rows,
        width=cols,
        count=1,
        dtype=data.dtype,
        crs="EPSG:32618",
        transform=transform,
    ) as dst:
        dst.write(data, 1)

    return bound_path, ref_path, line_path, pts_path, dem_path


def test_voronoi_workflow(tmp_path):
    bound_path, ref_path, line_path, pts_path, dem_path = create_synthetic_data(
        tmp_path
    )

    # Initialize
    vor = VoronoiGrid(bound_path)

    # Add features
    vor.add_refinement_area(ref_path, cell_size=50)
    vor.add_refinement_line(line_path, cell_size=20)
    vor.add_fixed_points(pts_path)

    # Add DEM
    # min_slope=0, max_slope=10. max_size=100, min_size=10
    vor.add_dem(
        dem_path, min_slope=0, max_slope=10, min_cell_size=10, max_cell_size=100
    )

    # Build
    # Use a coarse background
    vor.build(global_cell_size=200)

    # Check outputs
    assert vor.mesh_gdf is not None
    assert len(vor.mesh_gdf) > 0
    assert "geometry" in vor.mesh_gdf.columns

    # Export
    out_shp = tmp_path / "mesh.shp"
    vor.export_mesh(out_shp)
    assert out_shp.exists()

    # DISV Props
    props = vor.get_disv_properties()
    assert "vertices" in props
    assert "cell2d" in props
    assert len(props["vertices"]) > 0
    assert len(props["cell2d"]) == len(vor.mesh_gdf)


if __name__ == "__main__":
    # Manual run
    from pathlib import Path
    import shutil

    tmp = Path("tests/test_data_manual")
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    try:
        test_voronoi_workflow(tmp)
        print("Test passed!")
    except Exception as e:
        print(f"Test failed: {e}")
        import traceback

        traceback.print_exc()
