import geopandas as gpd
import os
import folium
import io
import fiona
import json
import pathlib
import numpy as np
import rasterio
import rasterio.mask
import tempfile

import shutil


def read_shp_from_zip(file):
    zipshp = io.BytesIO(open(file, 'rb').read())
    with fiona.BytesCollection(zipshp.read()) as src:
        crs = src.crs
        gdf = gpd.GeoDataFrame.from_features(src, crs=crs)
    return gdf


def write_shp_as_zip(zipLoc,zipDest,baseName):
    shutil.make_archive(base_dir=zipLoc,
        root_dir=zipDest,
        format='zip',
        base_name=baseName)

def shpFromZipAsFiona(file):
    zipshp = io.BytesIO(open(file, 'rb').read())
    fionaObj = fiona.BytesCollection(zipshp.read())
    return fionaObj


def remove_files_and_folder(path_to_file, folder=True):
    folder=os.path.dirname(path_to_file)

    if os.path.isfile(path_to_file):
        os.remove(path_to_file)
        print("File has been deleted")
    else:
        print("File does not exist")
    #for filename in os.listdir(folder):
    #    file_path=os.path.join(folder,filename)
        """
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    if folder:
        os.rmdir(folder)
    """
