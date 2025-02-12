import geopandas as gpd
import os, shutil
import folium
import io
import fiona
import json
import pathlib
import numpy as np
import rasterio
import rasterio.mask
import tempfile
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from shapely.ops import split, unary_union, cascaded_union, voronoi_diagram
import shutil
from collections import OrderedDict


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

def isMultiGeometry(geom):
    return isinstance(geom, (MultiPoint, MultiLineString, MultiPolygon))

#auxiliary funtion to intersect:
def intersectLimitLayer(discLayerGeom, modelDis):
    discGeomList = []  
    #generic 
    if isMultiGeometry(discLayerGeom):
        for partGeom in discLayerGeom.geoms:
            discGeomClip =  modelDis['limitGeometry'].intersection(partGeom)
            if not discGeomClip.is_empty:
                discGeomList.append(discGeomClip)
    else:
        discGeomClip =  modelDis['limitGeometry'].intersection(discLayerGeom)
        if not discGeomClip.is_empty:
            discGeomList.append(discGeomClip)

    unaryGeom = unary_union(discGeomList)

    if isMultiGeometry(unaryGeom):
        unaryFilter = [geom for geom in unaryGeom.geoms]
    else:
        unaryFilter = [unaryGeom]

    return unaryFilter    

def filterVertexCloseLimit(layerRef,layerGeom,pointList,modelDis):
    filterPointList = []
    pointObject = layerGeom.exterior.coords.xy
    pointList = list(zip(pointObject[0],pointObject[1]))

    if layerGeom.buffer(layerRef).within(modelDis['limitGeometry']):
        filterPointList = pointList
    else:
        for point in pointList:
            pointPoint = Point(point)
            if pointPoint.buffer(layerRef).within(modelDis['limitGeometry']):
                filterPointList.append(point)

    if len(filterPointList) > 1:
        filterPointListGeom = LineString(filterPointList)
    elif len(filterPointList) == 1:
        filterPointListGeom = Point(filterPointList)
    else:
        filterPointListGeom = None

    return filterPointList, filterPointListGeom 

def getFionaDictPoly(polyGeom, index):
    polyCoordList = []
    x,y = polyGeom.exterior.coords.xy
    polyCoordList.append(list(zip(x,y)))
    if polyGeom.interiors[:] != []:
        interiorList = []
        for interior in polyGeom.interiors:
            polyCoordList.append(interior.coords[:])
    feature = {
        "geometry": {'type':'Polygon',
                    'coordinates':polyCoordList},
        "properties": OrderedDict([("id",index)]),
    }
    return feature

def getFionaDictPoint(pointGeom, index):
    if isinstance(pointGeom[0], float):
        feature = {
                "geometry": {'type':'Point',
                            'coordinates':(pointGeom[0],pointGeom[1])},
                "properties": OrderedDict([("id",index)]),
            }
        return feature

def initiateOutputFolder(outputPath):
    if os.path.isdir(outputPath):
        print('The output folder %s exists and has been cleared'%outputPath)
        shutil.rmtree(outputPath)
        os.mkdir(outputPath)
    else:
        os.mkdir(outputPath)
        print('The output folder %s has been generated.'%outputPath)