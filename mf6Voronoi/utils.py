import geopandas as gpd
import pandas as pd
import os, shutil, time, json
from pathlib import Path
import io
import fiona
# /--------------- try to import dask_geopandas
try:
  import dask_geopandas as dgpd
  HAS_DASK = True
except ImportError:
  dgpd = None
  HAS_DASK = False
# ----------------/
import numpy as np
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from shapely.ops import unary_union
import shutil
from collections import OrderedDict

# _____________________________%% KMB ADDITIONS

def unique_rows(a,sort=True,return_inverse=False):
    '''
    Find unique rows and return indexes of unique rows
    '''
    a = np.ascontiguousarray(a)
    unique_a,uind,uinv = np.unique(a.view([('', a.dtype)]*a.shape[1]),return_index=True,return_inverse=True)
    if sort:    
        uord = [(uind==utemp).nonzero()[0][0] for utemp in np.sort(uind)]
        outorder = uind[uord]
    else:
        outorder = uind
    if return_inverse:
        return unique_a,uind,uinv
    else:
        return outorder

def xy_in_poly(xy,poly,use_dask=False,nproc=10,return_inds=False):
    
    points_array = np.array(xy)
    point_df = gpd.GeoDataFrame(gpd.pd.DataFrame(np.arange(len(xy)),columns=['id']),
                                         geometry=gpd.points_from_xy(points_array[:,0],points_array[:,1]))
    if use_dask:
        innerpoint_bool = dgpd.from_geopandas(point_df,nproc).within(poly).compute().values
    else:
        innerpoint_bool = point_df.within(poly).values
                 
    filterPointList = points_array[innerpoint_bool].tolist()
    
    if return_inds:
        return filterPointList,innerpoint_bool
    else:
        return filterPointList
        

def findIndex(var, coordArray, intervalNumber=10):
    for interval in range(intervalNumber):
        if var >= coordArray[interval] and var < coordArray[interval+1]:
            return interval
            break

            
def processVertexWithFilterCloseLimitDf(layerRef,layerGeoms,modelDis,use_dask=False,nproc=10):
    ###### For original ######
    #create a temporal dataframe for the current layers
    orgLayerDf = gpd.GeoDataFrame(np.arange(len(layerGeoms)),columns=['id'],geometry=layerGeoms)
    orgCoords = orgLayerDf.get_coordinates()
    orgLayerPtsDf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(orgCoords.x, orgCoords.y))

    # First collect all geometries within limitGeometry
    if use_dask:
        #create a temporal dataframe for the org layers
        orgLayerPtsInsideLimitBool = dgpd.from_geopandas(orgLayerPtsDf, nproc).buffer(layerRef).within(modelDis['limitGeometry'])
    else:
        orgLayerPtsInsideLimitBool = orgLayerPtsDf.buffer(layerRef).within(modelDis['limitGeometry'])

    #these are the points of the geometry that are inside limit geometry 
    orgPointMask = pd.Series(orgLayerPtsInsideLimitBool).values
    orgPointList =  orgLayerPtsDf.loc[orgPointMask].get_coordinates().values.tolist()    

    ###### For distributed ######
    # distLayerDf = orgLayerDf.segmentize(layerRef)
    # distCoords = distLayerDf.get_coordinates()
    # distLayerPtsDf = gpd.GeoDataFrame(geometry=gpd.points_from_xy(distCoords.x, distCoords.y))    
    # --- Distributed Points (Fixed Step Vectorized Interpolation) ---
    dist_raw_points = []
    for geom in layerGeoms:
        if geom.geom_type in ['Polygon', 'MultiPolygon']:
            # Use boundary/exterior length for polygons
            polys = geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]
            for poly in polys:
                # Collect exterior ring AND all interior hole rings
                rings = [poly.exterior] + list(poly.interiors)
                for ring in rings:
                    prog = np.arange(0, ring.length, layerRef)
                    dist_raw_points.extend(ring.interpolate(prog))

        elif geom.geom_type in ['LineString', 'MultiLineString']:
            lines = geom.geoms if geom.geom_type == 'MultiLineString' else [geom]
            for line in lines:
                prog = np.arange(0, line.length, layerRef)
                dist_raw_points.extend(line.interpolate(prog))

        elif geom.geom_type in ['Point', 'MultiPoint']:
            pts = geom.geoms if geom.geom_type == 'MultiPoint' else [geom]
            dist_raw_points.extend(pts)

    distLayerPtsDf = gpd.GeoDataFrame(geometry=dist_raw_points, crs=modelDis.get('crs'))

    # First collect all geometries within limitGeometry
    if use_dask:
        #create a temporal dataframe for the dist layers
        distLayerPtsInsideLimitBool = dgpd.from_geopandas(distLayerPtsDf, nproc).buffer(layerRef).within(modelDis['limitGeometry'])
    else:
        distLayerPtsInsideLimitBool = distLayerPtsDf.buffer(layerRef).within(modelDis['limitGeometry'])

    #these are the points of the geometry that are inside limit geometry 
    distPointMask = pd.Series(distLayerPtsInsideLimitBool).values
    distPointList =  distLayerPtsDf.loc[distPointMask].get_coordinates().values.tolist() 

    # Return list of Shapely Point geometries
    distPointList_asGeom = distLayerPtsDf.loc[distPointMask, 'geometry'].tolist()
    
    return orgPointList, distPointList, distPointList_asGeom

def save_zip(mesh_obj,out_fname,compresslevel=9):
    # source: https://stackoverflow.com/a/57758563
    with gzip.open(out_fname,'wt',encoding='UTF-8',compresslevel=compresslevel) as fout:
        json.dump(mesh_obj.disvDict,fout)
        

def read_zip(in_fname):
    # source: https://stackoverflow.com/a/57758563
    with gzip.open(in_fname,'rt',encoding='UTF-8') as fin:
        data = json.load(fin)
    
    return data

def getPolygonAndInteriors(polyGeom):
    exteriorInteriorPolys = [polyGeom] + [Polygon(ring) for ring in polyGeom.interiors]
    return exteriorInteriorPolys

def makeGeometry(geom_df,geom_col='geometry',gtype_col='geom_type'):
    
    starting_geom_type = geom_df[gtype_col].iloc[0]
    geom_in = geom_df[geom_col].values
    
    points = [p.coords[0] for p in geom_in]
    
    if starting_geom_type == 'Polygon' and len(points) > 2:
        out_geom = Polygon(points)
    elif len(points) > 1:  #len(filterPointList) > 1:
        out_geom = LineString(points)
    elif len(points) == 1: #len(filterPointList) == 1:
        out_geom = Point(points)
    else:
        out_geom = None
    return out_geom

# ----------------------------------------/

def readShpFromZip(file):
    zipshp = io.BytesIO(open(file, 'rb').read())
    with fiona.BytesCollection(zipshp.read()) as src:
        crs = src.crs
        gdf = gpd.GeoDataFrame.from_features(src, crs=crs)
    return gdf


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
        else: 
            return False

    unaryGeom = unary_union(discGeomList)

    if isMultiGeometry(unaryGeom):
        unaryFilter = [geom for geom in unaryGeom.geoms]
    else:
        unaryFilter = [unaryGeom]

    return unaryFilter    

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

###############
# Output functions
###############

def getVoronoiAsShp(modelDis, shapePath=''):
    print('\n/----Generation of the voronoi shapefile----/')
    start = time.time()
    schema_props = OrderedDict([("id", "int")])
    schema={"geometry": "Polygon", "properties": schema_props}

    outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                        crs = modelDis['crs'], schema=schema)
    
    for index, poly in enumerate(modelDis['voronoiRegions'].geoms):
        polyCoordList = []
        x,y = poly.exterior.coords.xy
        polyCoordList.append(list(zip(x,y)))
        if poly.interiors[:] != []:
            interiorList = []
            for interior in poly.interiors:
                polyCoordList.append(interior.coords[:])
        feature = {
            "geometry": {'type':'Polygon',
                        'coordinates':polyCoordList},
            "properties": OrderedDict([("id",index)]),
        }
        outFile.write(feature)
    outFile.close()

    
    end = time.time()
    print('\nTime required for voronoi shapefile: %.2f seconds \n'%(end - start), flush=True)

def getPolyAsShp(modelDis, circleList, shapePath=''):
  """Exports polygon geometries from modelDis to a Shapefile using GeoPandas.

  Handles Polygons, MultiPolygons, GeometryCollections, and dicts of geometries.
  """
  start = time.time()
  geom_data = modelDis[circleList]

  # 1. Extract raw geometries into a Python list
  geom_list = []
  if isinstance(geom_data, dict):
    for val in geom_data.values():
      if isinstance(val, (list, tuple)):
        geom_list.extend(val)
      else:
        geom_list.append(val)
  elif isinstance(geom_data, (list, tuple)):
    geom_list = list(geom_data)
  else:
    geom_list = [geom_data]

  # 2. Build GeoDataFrame
  gdf = gpd.GeoDataFrame(geometry=geom_list, crs=modelDis.get('crs'))

  # 3. Explode MultiGeometries / GeometryCollections into individual geometries
  gdf = gdf.explode(ignore_index=True)

  # 4. Filter only Polygon geometries (discards Lines/Points resulting from GeometryCollection)
  gdf = gdf[gdf.geometry.geom_type == 'Polygon'].copy()

  # 5. Add ID column and save to Shapefile
  if not gdf.empty:
    gdf['id'] = gdf.index.astype(str)
    gdf.to_file(shapePath, driver='ESRI Shapefile')
  else:
    print(f'[WARNING] No valid Polygon geometries found for {circleList}.')

  end = time.time()
  print(
      f'\nTime required for polygon shapefile ({shapePath}):'
      f' {end - start:.2f} seconds \n',
      flush=True,
  )

def getPointsAsShp(modelDis,pointList,shapePath=''):
    schema_props = OrderedDict([("id", "str")])
    schema={"geometry": "Point", "properties": schema_props}
    if shapePath != '':
        outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                            crs = modelDis['crs'], schema=schema)
        if isinstance(modelDis[pointList], dict):
            #print(modelDis[pointList].keys())
            for key, value in modelDis[pointList].items():
                for index, point in enumerate(value):
                    feature = getFionaDictPoint(point, index)
                    if feature != None:
                        outFile.write(feature)
                    else:
                        print('Something went wrong with %s'%point)
        else:
            for index, point in enumerate(modelDis[pointList]):
                feature = getFionaDictPoint(point, index)
                if feature != None:
                    outFile.write(feature)
                else:
                    print('Something went wrong with %s'%point)
        outFile.close()

def exportMeshBuildFeaturesToShp(modelDis, out_dir="debug_point_cloud"):
  """Exporta los puntos originales y distribuidos a archivos Shapefile."""
  if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

  print(f"\n[DEBUG] Exporting all point cloud categories to: {out_dir}")
    
  if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

  # 2. Puntos Distribuidos por capa (vertexDist)
  if "vertexDist" in modelDis and len(modelDis["vertexDist"]) > 0:
    dist_path = os.path.join(out_dir, "p1_vertexDist.shp")
    getPointsAsShp(modelDis, "vertexDist", dist_path)
    print(f"  - Puntos distribuidos exportados: {dist_path}")

  # 1. Puntos Originales (vertexOrg)
  if "vertexOrg" in modelDis and len(modelDis["vertexOrg"]) > 0:
    org_path = os.path.join(out_dir, "p2_vertexOrg.shp")
    getPointsAsShp(modelDis, "vertexOrg", org_path)
    print(f"  - Puntos originales exportados: {org_path}")

  # 2. Buffer Points (vertexBuffer)
  if "vertexBuffer" in modelDis and len(modelDis["vertexBuffer"]) > 0:
    buf_pts_path = os.path.join(out_dir, "p3_vertexBuffer.shp")
    getPointsAsShp(modelDis, "vertexBuffer", buf_pts_path)
    print(f"  - Buffer points exported: {buf_pts_path}")

  # 6. Combined Circle Buffers (circleUnion)
  if "circleUnion" in modelDis and modelDis["circleUnion"] is not None:
    buf_poly_path = os.path.join(out_dir, "p4_circleUnion.shp")
    getPolyAsShp(modelDis, "circleUnion", buf_poly_path)
    print(f"  - Circle buffer polygons exported: {buf_poly_path}")

  # 6. Combined Circle Buffers (circleUnion)
  if "circleUnionByStep" in modelDis and modelDis["circleUnionByStep"] is not None:
    buf_poly_path = os.path.join(out_dir, "p4x_circleUnionByStep.shp")
    getPolyAsShp(modelDis, "circleUnionByStep", buf_poly_path)
    print(f"  - Circle buffer polygons by step exported: {buf_poly_path}")

  # 6. Combined Circle Buffers (circleUnion)
  if "bufferPerCellSize" in modelDis and modelDis["bufferPerCellSize"] is not None:
    buf_poly_path = os.path.join(out_dir, "p4x_bufferPerCellSize.shp")
    getPolyAsShp(modelDis, "bufferPerCellSize", buf_poly_path)
    print(f"  - buffer polygons before unary union: {buf_poly_path}")

  # 7. Circle Buffers with Interiors (circleUnionInteriors)
  if "circleUnionInteriors" in modelDis and modelDis["circleUnionInteriors"] is not None:
    buf_int_path = os.path.join(out_dir, "p5_circleUnionInteriors.shp")
    getPolyAsShp(modelDis, "circleUnionInteriors", buf_int_path)
    print(f"  - Circle buffer interior polygons exported: {buf_int_path}")
    
  # 3. Max Refinement Coarse Grid Points (vertexMaxRef)
  if "vertexMaxRef" in modelDis and len(modelDis["vertexMaxRef"]) > 0:
    max_path = os.path.join(out_dir, "p6_vertexMaxRef.shp")
    getPointsAsShp(modelDis, "vertexMaxRef", max_path)
    print(f"  - Max refinement points exported: {max_path}")

  # 4. Min Refinement Points (vertexMinRef)
  if "vertexMinRef" in modelDis and len(modelDis["vertexMinRef"]) > 0:
    min_path = os.path.join(out_dir, "p7_vertexMinRef.shp")
    getPointsAsShp(modelDis, "vertexMinRef", min_path)
    print(f"  - Min refinement points exported: {min_path}")

  # 6. Combined Circle Buffers (circleUnion)
  if "pointsMaxRefPoly" in modelDis and modelDis["pointsMaxRefPoly"] is not None:
    buf_poly_path = os.path.join(out_dir, "p8_pointsMaxRefPoly.shp")
    getPolyAsShp(modelDis, "pointsMaxRefPoly", buf_poly_path)
    print(f"  - zona donde aplica el refinamiento máximo: {buf_poly_path}")

  # 1. Total Combined Points (vertexTotal)
  if "vertexTotal" in modelDis and len(modelDis["vertexTotal"]) > 0:
    pts_path = os.path.join(out_dir, "p9_vertexTotal.shp")
    getPointsAsShp(modelDis, "vertexTotal", pts_path)
    print(f"  - Total points exported: {pts_path}")

#########
# miscelaneous functions
#########


    
def copyTemplate(templateType, prefix = ''):
    utilsPath = os.path.realpath(__file__)
    examplePath = os.path.join(os.path.dirname(utilsPath),'examples','notebooks')
    jsonPath = os.path.join(os.path.dirname(utilsPath),'examples','json','templates.json')
    workingPath = os.getcwd()
    
    with open(jsonPath, 'r') as file:
        templateDict = json.load(file)

    try:
        tempDict = templateDict[templateType]
        srcPath = str(Path(os.path.join(examplePath, tempDict["template"])))
        if prefix != '':
            dstPath = str(Path(os.path.join(workingPath, prefix+'_'+tempDict["template"])))
        else:
            dstPath = str(Path(os.path.join(workingPath, tempDict["template"])))
        shutil.copy2(srcPath,dstPath)
    except KeyError:
        print("The template: %s doesn't exists capullo"%templateType)

def listTemplates():
    utilsPath = os.path.realpath(__file__)
    jsonPath = os.path.join(os.path.dirname(utilsPath),'examples','json','templates.json')
    with open(jsonPath, 'r') as file:
        templateDict = json.load(file)

    print("/-------- List of available mf6Voronoi templates --------/\n")

    for key in templateDict.keys():
        print("Nr %d: %s"%(templateDict[key]["index"],key))
        print("    File: %s"%(templateDict[key]["template"]))
        print("    Description: %s\n"%(templateDict[key]["desc"]))

def isRunningInJupyter():
    try:
        from IPython import get_ipython
        shell = get_ipython().__class__.__name__
        return shell == 'ZMQInteractiveShell'
    except (NameError, ImportError):
        return False
    
def printBannerHtml():
    from IPython.display import display, HTML

    html_content = """
    <link href="https://fonts.googleapis.com/css2?family=Anton&display=swap" rel="stylesheet">

    <style>
        .styled-text {
        font-family: 'Anton', Impact, sans-serif;
        font-size: 32px;
        font-weight: bold;
        font-style: italic;
        }
    </style>

    <div>
    <a href="https://hatarilabs.com" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/png/hatarilabs.png" alt="Hatarilabs" width="200" height="200"></a> 
            <p class="styled-text">mf6Voronoi will have a web version in 2028</p>
    </div>

    <table border="0px">
    <tbody>
    <tr>
        <td><h3>Follow us:</h3></td>
        <td><a href="https://www.linkedin.com/company/hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-linkedin.svg" alt="Hatarilabs"></a></td>
        <td><a href="https://www.facebook.com/hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-facebook.svg" alt="Hatarilabs"></a></td>
        <td><a href="https://www.instagram.com/hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-instagram.svg" alt="Hatarilabs"></a></td>
        <td><a href="https://www.youtube.com/hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-youtube.svg" alt="Hatarilabs"></a></td>
        <td><a href="https://www.tiktok.com/@_hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-tiktok.svg" alt="Hatarilabs"></a></td>
        <td><a href="https://x.com/hatarilabs" target="_blank">
            <img src="https://olivosbellaterra.com/static/img/svg/icons8-twitterx.svg" alt="Hatarilabs"></a></td>
    </tr>
    </tbody>
    </table>

    """

    display(HTML(html_content))

def printBannerText():
    print('''
                                                                                                    
*mSi                                                                                       
gQQ>                                                                                       
dQU;                                 +|:                                     :v)_          
;PQm'                                %B$s                                    .gQMe          
PYQ7.                               -3QE_                                     <e}'          
c8Qx                                '$RT                                                    
?HM"   )7yw1=       .)r]jJfzi.   `=>!QDuvxxi_   `<s[LCwe<    ,>seua:  ^!C3o' `eur           
oRk= vdZ6qDQE"     ]PF)/+vJBNe`  :l{6Q8!I![s' .ebJ<//%3MD]   )fffQDv.ebPhZY/ QQQ#           
JQS'7b]_  ?Q$r     EWy     3Qg^     pQX       :GWj    _mQ5'     lQ&TT4v   .  rQDl           
5QnCV/    ]QZ<      :"/iiss4Q5:    -GQS        .:|/<v%IPQJ`     !Qk[h;       ?QK"           
.SQXd^    _JQh:    -*53e*ppaDQL.    _&QF       '!6fa{vzzMQa      7Q@q|        1Q@,           
;4QWi     :pQp    _dQY-   xm@Qa     _OQd      ^XQV.   rhHQ!     .wQBo         [QK^           
KkQ6'     '2QO}vc/:&QDa}75L<2Qgx="= .nQMCv)%I"UUQ81}j57>SQhi=|; .dQA'         %QQCi%)_       
/jJ>       82mw[i: /zmVFa|  ;t53j}+  `1mVpn!>_ UuSh21/  =oFy7{; .z#I          .nm57r/.       
                                                                                                                                                                                                                                      
''')
# %%
