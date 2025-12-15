import os
import json
import warnings
import numpy as np
import geopandas as gpd
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.features import rasterize
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, box
from shapely.ops import unary_union, voronoi_diagram
from scipy.spatial import Voronoi, cKDTree
import matplotlib.pyplot as plt

class VoronoiGrid:
    """
    A class to generate a Voronoi grid for MODFLOW 6 DISV packages.
    """

    def __init__(self, boundary_path, min_cell_size=None, max_cell_size=None, crs=None):
        """
        Initialize the VoronoiGrid.

        Parameters
        ----------
        boundary_path : str
            Path to the boundary vector file (SHP, GPKG).
            The first polygon in this file defines the model domain.
        min_cell_size : float, optional
            Minimum allowed cell size (point spacing).
            Used for proximity checks and warnings.
        max_cell_size : float, optional
            Maximum allowed cell size.
            Used as the default background grid spacing.
        crs : str or dict, optional
            Coordinate reference system (e.g., 'EPSG:32618').
            If None, it is read from the file.
        """
        self.boundary_gdf = gpd.read_file(boundary_path)
        if crs:
            self.boundary_gdf = self.boundary_gdf.to_crs(crs)

        # Ensure boundary is a single geometry (union if multiple)
        self.boundary_geom = unary_union(self.boundary_gdf.geometry)
        if self.boundary_geom.geom_type == 'MultiPolygon':
            warnings.warn("Boundary is a MultiPolygon. This might cause issues with simple containment checks.")

        self.min_cell_size = min_cell_size
        self.max_cell_size = max_cell_size

        self.points = []
        self.refinements = []
        self.dem_config = None
        self.mesh_gdf = None
        self.disv_props = None

    def add_refinement_area(self, path, cell_size):
        """
        Add a polygon refinement area.

        Parameters
        ----------
        path : str
            Path to vector file defining the area.
        cell_size : float
            Target cell size (approximate spacing of points) inside this area.
        """
        gdf = gpd.read_file(path)
        if self.boundary_gdf.crs and gdf.crs != self.boundary_gdf.crs:
            gdf = gdf.to_crs(self.boundary_gdf.crs)

        self.refinements.append({
            'type': 'area',
            'geometry': unary_union(gdf.geometry),
            'cell_size': cell_size
        })

    def add_refinement_line(self, path, cell_size):
        """
        Add a line refinement.

        Parameters
        ----------
        path : str
            Path to vector file defining lines (e.g., rivers, faults).
        cell_size : float
            Target point spacing along the line.
        """
        gdf = gpd.read_file(path)
        if self.boundary_gdf.crs and gdf.crs != self.boundary_gdf.crs:
            gdf = gdf.to_crs(self.boundary_gdf.crs)

        self.refinements.append({
            'type': 'line',
            'geometry': unary_union(gdf.geometry),
            'cell_size': cell_size
        })

    def add_fixed_points(self, points, crs=None):
        """
        Add fixed points that must be included as generating seeds (nodes).

        Parameters
        ----------
        points : str, list of tuple, or list of Point
            - If str: Path to vector file defining points.
            - If list of tuple: [(x1, y1), (x2, y2), ...]
            - If list of Point: [Point(x1, y1), ...]
        crs : str, optional
            CRS of the input points if provided as list.
            If None, assumes same as boundary.
            If points is a file, CRS is read from file.
        """
        new_points = []

        if isinstance(points, str):
            # File path
            gdf = gpd.read_file(points)
            if self.boundary_gdf.crs and gdf.crs != self.boundary_gdf.crs:
                gdf = gdf.to_crs(self.boundary_gdf.crs)

            for geom in gdf.geometry:
                if geom.geom_type == 'Point':
                    new_points.append((geom.x, geom.y))
                elif geom.geom_type == 'MultiPoint':
                    for p in geom.geoms:
                        new_points.append((p.x, p.y))

        elif isinstance(points, list):
            # List of points
            # Check type of first element
            if not points:
                return

            first = points[0]
            if isinstance(first, Point):
                # Convert to tuples
                # Handle CRS if needed (basic support)
                if crs and self.boundary_gdf.crs and crs != self.boundary_gdf.crs:
                    # Creating a GDF to reproject is easiest
                    gdf = gpd.GeoDataFrame(geometry=points, crs=crs)
                    gdf = gdf.to_crs(self.boundary_gdf.crs)
                    new_points = [(p.x, p.y) for p in gdf.geometry]
                else:
                    new_points = [(p.x, p.y) for p in points]
            else:
                # tuples/lists
                # Assume provided in correct CRS if crs is None, or reproject manually?
                # For simplicity, if tuple, assumes project CRS unless crs is provided
                if crs and self.boundary_gdf.crs and crs != self.boundary_gdf.crs:
                     pts_geom = [Point(x, y) for x, y in points]
                     gdf = gpd.GeoDataFrame(geometry=pts_geom, crs=crs)
                     gdf = gdf.to_crs(self.boundary_gdf.crs)
                     new_points = [(p.x, p.y) for p in gdf.geometry]
                else:
                    new_points = [(x, y) for x, y in points]

        self.points.extend(new_points)

    def add_dem(self, path, min_slope, max_slope, min_cell_size, max_cell_size):
        """
        Configure DEM-based refinement.

        Parameters
        ----------
        path : str
            Path to the DEM raster (GeoTIFF).
        min_slope : float
            Slope threshold for maximum cell size (flat areas).
        max_slope : float
            Slope threshold for minimum cell size (steep areas).
        min_cell_size : float
            Cell size at max_slope (steepest).
        max_cell_size : float
            Cell size at min_slope (flattest).
        """
        self.dem_config = {
            'path': path,
            'min_slope': min_slope,
            'max_slope': max_slope,
            'min_cell_size': min_cell_size,
            'max_cell_size': max_cell_size
        }

    def _generate_points_in_polygon(self, poly, spacing):
        """Generate a regular grid of points inside a polygon."""
        minx, miny, maxx, maxy = poly.bounds
        x_coords = np.arange(minx, maxx, spacing)
        y_coords = np.arange(miny, maxy, spacing)

        xx, yy = np.meshgrid(x_coords, y_coords)
        pts = np.vstack([xx.ravel(), yy.ravel()]).T

        valid_points = []
        # Use shapely.prepared for efficiency
        from shapely.prepared import prep
        prepared_poly = prep(poly)

        for x, y in pts:
            if prepared_poly.contains(Point(x, y)):
                valid_points.append((x, y))
        return valid_points

    def _generate_points_along_line(self, line, spacing):
        """Generate points along a line string."""
        length = line.length
        num_points = int(length / spacing)
        if num_points == 0:
            return []

        points = [line.interpolate(d) for d in np.linspace(0, length, num_points)]
        return [(p.x, p.y) for p in points]

    def _generate_dem_points(self):
        """Generate points based on DEM slope."""
        if not self.dem_config:
            return []

        print("Generating points from DEM...")
        with rasterio.open(self.dem_config['path']) as src:
            try:
                out_image, out_transform = mask(src, [self.boundary_geom], crop=True)
                out_image = out_image[0]
            except ValueError:
                warnings.warn("DEM does not overlap with boundary.")
                return []

            dx = out_transform.a
            dy = -out_transform.e
            grad_y, grad_x = np.gradient(out_image, dy, dx)
            slope = np.sqrt(grad_x**2 + grad_y**2)

            min_s, max_s = self.dem_config['min_slope'], self.dem_config['max_slope']
            min_sz, max_sz = self.dem_config['min_cell_size'], self.dem_config['max_cell_size']

            slope_clipped = np.clip(slope, min_s, max_s)

            if max_s == min_s:
                target_size = np.full(slope_clipped.shape, max_sz)
            else:
                target_size = max_sz - (slope_clipped - min_s) / (max_s - min_s) * (max_sz - min_sz)

            pixel_area = abs(dx * dy)
            prob = pixel_area / (target_size**2)

            rand_grid = np.random.random(target_size.shape)

            valid_mask = (out_image != src.nodata) if src.nodata is not None else np.ones(out_image.shape, dtype=bool)

            selected = np.where((rand_grid < prob) & valid_mask)
            rows, cols = selected
            xs, ys = rasterio.transform.xy(out_transform, rows, cols, offset='center')

            return list(zip(xs, ys))

    def build(self, global_cell_size=None):
        """
        Generate the Voronoi mesh.

        Parameters
        ----------
        global_cell_size : float, optional
            Background cell size for the whole domain.
            If None, uses self.max_cell_size from init.
            If both are None, no background grid is generated.
        """
        print("Building Voronoi Grid...")

        if global_cell_size is None:
            global_cell_size = self.max_cell_size

        # 1. Fixed Points
        fixed_points = list(set(self.points)) # Remove dupes

        # Check proximity of fixed points against min_cell_size
        if self.min_cell_size and len(fixed_points) > 1:
            tree = cKDTree(fixed_points)
            # Query nearest neighbor (k=2 because the nearest is the point itself)
            dists, _ = tree.query(fixed_points, k=2)
            # dists[:, 1] is the distance to the nearest neighbor
            min_d = np.min(dists[:, 1])
            if min_d < self.min_cell_size:
                warnings.warn(f"Some fixed points are closer than min_cell_size ({self.min_cell_size}). Closest pair distance: {min_d:.2f}")

        # 2. Generated Points
        generated_points = set()

        # Refinements
        for ref in self.refinements:
            geom = ref['geometry']
            spacing = ref['cell_size']

            if ref['type'] == 'area':
                pts = self._generate_points_in_polygon(geom, spacing)
                generated_points.update(pts)
            elif ref['type'] == 'line':
                if geom.geom_type == 'MultiLineString':
                    for line in geom.geoms:
                        pts = self._generate_points_along_line(line, spacing)
                        generated_points.update(pts)
                else:
                    pts = self._generate_points_along_line(geom, spacing)
                    generated_points.update(pts)

        # Background
        if global_cell_size:
            bg_points = self._generate_points_in_polygon(self.boundary_geom, global_cell_size)
            generated_points.update(bg_points)

        # DEM
        dem_pts = self._generate_dem_points()
        generated_points.update(dem_pts)

        # Filter generated points against Fixed Points
        final_points = list(fixed_points)

        if generated_points:
            gen_pts_list = list(generated_points)
            if fixed_points and self.min_cell_size:
                print("Filtering generated points against fixed points...")
                tree = cKDTree(fixed_points)
                # query generated points against fixed points tree
                # if dist < min_cell_size, discard
                dists, _ = tree.query(gen_pts_list)

                valid_gen_pts = []
                for i, p in enumerate(gen_pts_list):
                    if dists[i] >= self.min_cell_size:
                        valid_gen_pts.append(p)

                final_points.extend(valid_gen_pts)
                print(f"Filtered {len(gen_pts_list) - len(valid_gen_pts)} generated points due to proximity.")
            else:
                final_points.extend(gen_pts_list)

        # Remove duplicates again just in case
        final_points = list(set(final_points))

        if len(final_points) < 3:
            raise ValueError("Not enough points to generate Voronoi diagram (need at least 3).")

        print(f"Generating Voronoi for {len(final_points)} points...")

        # Voronoi
        vor = Voronoi(final_points)
        from shapely.geometry import MultiPoint
        mp = MultiPoint(final_points)

        # Use shapely's voronoi_diagram which is robust
        # envelope needs to be large enough
        envelope_buffer = max(global_cell_size if global_cell_size else 1000, 1000)
        envelope = self.boundary_geom.buffer(envelope_buffer)
        regions = voronoi_diagram(mp, envelope=envelope)

        # Clip to boundary
        clipped_regions = []
        for region in regions.geoms:
            clipped = region.intersection(self.boundary_geom)
            if not clipped.is_empty:
                if clipped.geom_type == 'Polygon':
                    clipped_regions.append(clipped)
                elif clipped.geom_type == 'MultiPolygon':
                    clipped_regions.extend(clipped.geoms)

        self.mesh_gdf = gpd.GeoDataFrame(geometry=clipped_regions, crs=self.boundary_gdf.crs)
        self.mesh_gdf['id'] = range(len(self.mesh_gdf))

        print(f"Mesh generated with {len(self.mesh_gdf)} cells.")

    def export_mesh(self, path, driver=None):
        """Export mesh to file."""
        if self.mesh_gdf is None:
            raise ValueError("Mesh not generated. Call build() first.")

        path_str = str(path)
        if driver is None:
            if path_str.lower().endswith('.gpkg'):
                driver = 'GPKG'
            else:
                driver = 'ESRI Shapefile'

        self.mesh_gdf.to_file(path, driver=driver)
        print(f"Mesh exported to {path}")

    def get_disv_properties(self):
        """
        Calculate DISV properties (vertices, cell2d, etc.).
        Returns a dictionary suitable for MF6 DISV.
        """
        if self.mesh_gdf is None:
            raise ValueError("Mesh not generated.")

        print("Calculating DISV properties...")

        all_coords = []
        poly_coords_list = []

        for poly in self.mesh_gdf.geometry:
            coords = list(poly.exterior.coords)[:-1]
            all_coords.extend(coords)
            poly_coords_list.append(coords)

        coords_arr = np.array(all_coords)
        unique_coords, indices = np.unique(coords_arr, axis=0, return_inverse=True)

        vertices = []
        for i, (x, y) in enumerate(unique_coords):
            vertices.append([i, x, y])

        cell2d = []
        current_idx = 0

        for i, poly in enumerate(self.mesh_gdf.geometry):
            n_vert = len(poly_coords_list[i])
            cell_verts = indices[current_idx : current_idx + n_vert]
            current_idx += n_vert

            centroid = poly.centroid
            xc, yc = centroid.x, centroid.y

            cell_record = [i, xc, yc, n_vert] + cell_verts.tolist()
            cell2d.append(cell_record)

        self.disv_props = {
            'nvert': len(vertices),
            'vertices': vertices,
            'cell2d': cell2d,
            'ncpl': len(cell2d)
        }

        return self.disv_props

    def plot(self, interactive=False, **kwargs):
        """
        Plot the mesh.

        Parameters
        ----------
        interactive : bool
            If True, use geopandas.explore() (requires folium).
            If False, use matplotlib plot.
        **kwargs :
            Arguments passed to the plotting function.
        """
        if self.mesh_gdf is None:
            print("Mesh not built.")
            return

        if interactive:
            # Create base map with boundary
            m = self.boundary_gdf.explore(
                name="Boundary", style_kwds={"fill": False, "color": "black", "weight": 2}
            )

            # Add Mesh
            self.mesh_gdf.explore(
                m=m,
                name="Voronoi Mesh",
                style_kwds={"fill": False, "color": "blue", "weight": 0.5},
                **kwargs,
            )

            # Add Fixed Points
            if self.points:
                pts_geom = [Point(x, y) for x, y in self.points]
                pts_gdf = gpd.GeoDataFrame(geometry=pts_geom, crs=self.boundary_gdf.crs)
                pts_gdf.explore(
                    m=m,
                    name="Fixed Points",
                    color="red",
                    marker_kwds={"radius": 2},
                )

            # Add Refinements
            # We construct a GDF for all refinements to show them
            ref_geoms = []
            ref_types = []
            for ref in self.refinements:
                ref_geoms.append(ref["geometry"])
                ref_types.append(f"Refinement ({ref['type']}) - {ref['cell_size']}")

            if ref_geoms:
                ref_gdf = gpd.GeoDataFrame(
                    {"type": ref_types, "geometry": ref_geoms},
                    crs=self.boundary_gdf.crs,
                )
                ref_gdf.explore(
                    m=m,
                    name="Refinements",
                    style_kwds={"color": "green", "weight": 1, "dashArray": "5, 5"},
                )

            # Add Layer Control
            import folium

            folium.LayerControl().add_to(m)

            return m
        else:
            fig, ax = plt.subplots(figsize=(10, 10))
            self.boundary_gdf.plot(
                ax=ax,
                facecolor="none",
                edgecolor="black",
                linewidth=2,
                label="Boundary",
            )
            self.mesh_gdf.plot(
                ax=ax,
                facecolor="none",
                edgecolor="blue",
                linewidth=0.5,
                alpha=0.5,
                **kwargs,
            )

            if self.points:
                x, y = zip(*self.points)
                ax.scatter(x, y, c="red", s=5, label="Fixed Points")

            # Simple plot for refinements in static mode too?
            for ref in self.refinements:
                 gpd.GeoSeries([ref["geometry"]]).plot(ax=ax, color='green', alpha=0.3, linestyle='--')

            plt.title(f"Voronoi Mesh ({len(self.mesh_gdf)} cells)")
            plt.show()
