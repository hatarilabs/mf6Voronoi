import geopandas as gpd
import os, shutil, time
from pathlib import Path
import io
import fiona
import numpy as np
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from shapely.ops import unary_union
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
        else: 
            return False

    unaryGeom = unary_union(discGeomList)

    if isMultiGeometry(unaryGeom):
        unaryFilter = [geom for geom in unaryGeom.geoms]
    else:
        unaryFilter = [unaryGeom]

    return unaryFilter    

def processVertexFilterCloseLimit(layerRef,layerGeom,modelDis,vertexType):
    #first conditional
    if vertexType == 'Org':
        if layerGeom.geom_type == 'Polygon':
            pointObject = layerGeom.exterior.coords.xy
            pointList = list(zip(pointObject[0],pointObject[1]))
        elif layerGeom.geom_type == 'LineString':
            pointObject = layerGeom.coords.xy
            pointList = list(zip(pointObject[0],pointObject[1]))
        elif layerGeom.geom_type == 'Point':
            #covered in the third conditional
            pass
        else:
            print('Something went wrong with the org vertex')
    elif vertexType == 'Dist':
        pointList = []
        if layerGeom.geom_type == 'Polygon':
            polyLength = layerGeom.exterior.length
            pointProg = np.arange(0,polyLength,layerRef)
            for prog in pointProg:
                pointXY = list(layerGeom.exterior.interpolate(prog).xy)
                pointList.append([pointXY[0][0],pointXY[1][0]])
        elif layerGeom.geom_type == 'LineString':
            lineLength = layerGeom.length
            pointProg = np.arange(0,lineLength,layerRef)
            for prog in pointProg:
                pointXY = list(layerGeom.interpolate(prog).xy)
                pointList.append([pointXY[0][0],pointXY[1][0]])
        elif layerGeom.geom_type == 'Point':
            #covered in the third conditional
            pass
        else:
            print('Something went wrong with the dist vertex')

    #second conditional
    if layerGeom.geom_type == 'Polygon' or layerGeom.geom_type == 'LineString':
        if layerGeom.buffer(layerRef).within(modelDis['limitGeometry']):
            filterPointList = pointList
        else:
            filterPointList = []
            for point in pointList:
                pointPoint = Point(point)
                if pointPoint.buffer(layerRef).within(modelDis['limitGeometry']):
                    filterPointList.append(point)
        
        if layerGeom.geom_type == 'Polygon' and len(filterPointList) > 2:
            filterPointListGeom = Polygon(filterPointList)
        elif len(filterPointList) > 1:
            filterPointListGeom = LineString(filterPointList)
        elif len(filterPointList) == 1:
            filterPointListGeom = Point(filterPointList)
        else:
            filterPointListGeom = None

    #third conditional
    elif layerGeom.geom_type == 'Point':
        pointObject = layerGeom.coords.xy
        point = (pointObject[0][0],pointObject[1][0])
        if layerGeom.buffer(layerRef).within(modelDis['limitGeometry']):
            filterPointList = [point]
            filterPointListGeom = layerGeom
        else:
            filterPointList = []
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

def getPolyAsShp(modelDis,circleList,shapePath=''):
    start = time.time()
    schema_props = OrderedDict([("id", "str")])
    schema={"geometry": "Polygon", "properties": schema_props}
    
    outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                        crs = modelDis['crs'], schema=schema)
    
    if isinstance(modelDis[circleList], dict):
        for key, value in modelDis[circleList].items():
            if isMultiGeometry(value):
                for index, poly in enumerate(value.geoms):
                    feature = getFionaDictPoly(poly, index)
                    outFile.write(feature)

    if isMultiGeometry(modelDis[circleList]):
        for index, poly in enumerate(modelDis[circleList].geoms):
            feature = getFionaDictPoly(poly, index)
            outFile.write(feature)
    else:
        poly = modelDis[circleList]
        feature = getFionaDictPoly(poly, '1')
        outFile.write(feature)
    outFile.close()
    
    end = time.time()
    print('\nTime required for voronoi shapefile: %.2f seconds \n'%(end - start), flush=True)

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

#########
# miscelaneous functions
#########

def printHeader():
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
    
def copyTemplate(templateType, prefix = ''):
    utilsPath = os.path.realpath(__file__)
    examplePath = os.path.join(os.path.dirname(utilsPath),'examples','notebooks')
    workingPath = os.getcwd()

    def generateSrcDstPath():
        srcPath = str(Path(os.path.join(examplePath, templateName)))
        if prefix != '':
            dstPath = str(Path(os.path.join(workingPath, prefix+'_'+templateName)))
        else:
            dstPath = str(Path(os.path.join(workingPath, templateName)))
        return srcPath,dstPath

    if templateType == 'generateVoronoi':
        templateName = 'p1_'+templateType+'.ipynb'
        srcPath,dstPath = generateSrcDstPath()
        shutil.copy2(srcPath,dstPath)
    elif templateType == 'modelCreation':
        templateName = 'p2_'+templateType+'.ipynb'
        srcPath,dstPath = generateSrcDstPath()
        shutil.copy2(srcPath,dstPath)
    elif templateType == 'obsCalculated':
        templateName = 'p3_'+templateType+'.ipynb'
        srcPath,dstPath = generateSrcDstPath()
        shutil.copy2(srcPath,dstPath)
    else:
        print("The template you want doesn't exists capullo")
