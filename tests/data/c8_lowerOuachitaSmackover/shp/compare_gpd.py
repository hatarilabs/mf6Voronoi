# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 16:43:37 2026

@author: kmbefus
"""

import os,sys
from osgeo import gdal

vor_dir = r'F:\OneDrive - University of Arkansas\research\python\git\mf6Voronoi'
sys.path.insert(1,vor_dir)
from mf6Voronoi.geoVoronoi import createVoronoi

#%%
wdir = r'D:\data\ar_flood\mesh'
test_limit_fname = os.path.join(vor_dir,'mf6Voronoi','examples','extent_08040201.shp')
test_layer_fname = os.path.join(vor_dir,'mf6Voronoi','examples','riparian_08040201.shp')

verbose=True
max_cell_size = 1000 # meters
min_cell_size = 200

for gpd_on in [True,False][:1]: # not checking for dask speed ups here

    if gpd_on:
        scenario_name = 'gpd_on'
    else:
        scenario_name = 'original'
    
    out_dir = os.path.join(wdir,'vor_test','{}'.format(scenario_name))
    
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    
    vorMesh = createVoronoi(meshName='regionalModel',maxRef=max_cell_size, multiplier=1.1,
                            use_gpd=gpd_on)
    
    vor_shp = os.path.join(out_dir,'{}.shp'.format(vorMesh.modelDis['meshName']))
    
    # Open limit layers and refinement definition layers
    vorMesh.addLimit('basin',test_limit_fname)
    vorMesh.addLayer('river',test_layer_fname,min_cell_size)
    
    # Generate point pair array
    print("Building interpolated vertex lists")
    vorMesh.generateOrgDistVertices()
    
    # Try turning off gpd_on here to see if it helps pinpoint the gpd issue
    vorMesh.settings['use_gpd'] = False
    
    # Generate the point cloud and voronoi mesh
    vorMesh.createPointCloud(verbose=verbose,use_banner=False)
    vorMesh.generateVoronoi(shapePath=vor_shp)