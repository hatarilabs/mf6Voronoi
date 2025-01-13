import numpy as np
import copy
import tqdm, time
#from .lloydRelax import Field
#import geopandas as gpd
from scipy.spatial import Voronoi
#import geospatial libraries
import fiona
from tqdm import tqdm
from shapely.ops import split, unary_union, voronoi_diagram
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPoint, MultiPolygon
from collections import OrderedDict
from .utils import shpFromZipAsFiona

class createVoronoi():
    def __init__(self):
        self.discArrays = {}
        self.modelDis = {}

    def addLimit(self, name, shapePath):
        #Create the model limit
        limitShape = shpFromZipAsFiona(shapePath)

        #check if the geometry geometry type is polygon
        if limitShape[0]['geometry']['type'] != 'Polygon':
            print('A polygon layer is needed')
            exit()
        elif len(limitShape) > 1:
            print('Just one polygon is required')
            exit()

        #get all dimensions from the shapefile
        limitGeom = Polygon(limitShape[0]['geometry']['coordinates'][0])
        limitBounds = limitGeom.bounds
        self.modelDis['xMin'], self.modelDis['xMax'] = [limitBounds[i] for i in [0,2]]
        self.modelDis['yMin'], self.modelDis['yMax'] = [limitBounds[i] for i in [1,3]]
        self.modelDis['xDim'] = limitBounds[2] - limitBounds[0]
        self.modelDis['yDim'] = limitBounds[3] - limitBounds[1]
        self.modelDis['limitShape'] = limitShape
        self.modelDis['limitGeometry'] = limitGeom
        self.modelDis['crs'] = limitShape.crs

    def defineParameters(self, maxRef=500, minRef=50, stages=5):
        #Define the refinement sizes
        self.modelDis['refSizeList'] = np.linspace(maxRef,minRef,stages)
        line1 = '/---Cell discretization---/'
        line2 = 'Max. refinement: %d m.'%self.modelDis['refSizeList'].max()
        line3 = 'Min. refinement: %d m.'%self.modelDis['refSizeList'].min()
        line4 = 'Cell size list:'
        line5 = str(self.modelDis['refSizeList'])
        outMessage = line1 + '<br/>' + line2 + '<br/>' + line3 + '<br/>' + line4 + '<br/>' + line5
        return outMessage

    def addLayer(self, name, shapePath):
        #Add layers for mesh definition
        #This feature also clips and store the geometry
        shape = shpFromZipAsFiona(shapePath)
        discGeomList =[]
        if shape.schema['geometry'] == 'Polygon':
            for poly in shape:
                if poly['geometry']['type'] == 'Polygon':
                    for spart in poly['geometry']['coordinates']:
                        polyGeom = Polygon(spart)
                        discGeomClip =  self.modelDis['limitGeometry'].intersection(polyGeom)
                        if not discGeomClip.is_empty:
                            discGeomList.append(discGeomClip)
                elif poly['geometry']['type'] == 'MultiPolygon':
                    for mpart in poly['geometry']['coordinates']:
                        for spart in mpart:
                            polyGeom = Polygon(spart)
                            discGeomClip =  self.modelDis['limitGeometry'].intersection(polyGeom)
                            if not discGeomClip.is_empty:
                                discGeomList.append(discGeomClip)
            unaryPoly = unary_union(discGeomList)
            if unaryPoly.geom_type == 'MultiPolygon':
                unaryMulti = unaryPoly
            elif unaryPoly.geom_type == 'Polygon':
                unaryMulti = MultiPolygon([unaryPoly])
            discGeomList = [poly for poly in unaryMulti.geoms]
            self.discArrays[name] = {'type':'Polygon', 'geomList':discGeomList}
            self.modelDis['discGeomClip'] = unary_union(discGeomList)
        elif shape.schema['geometry'] == 'LineString':
            for line in shape:
                if shape[0]['geometry']['type'] == 'LineString':
                    lineGeom = LineString(line['geometry']['coordinates'])
                    discGeomClip =  self.modelDis['limitGeometry'].intersection(lineGeom)
                    if not discGeomClip.is_empty:
                        discGeomList.append(discGeomClip)
                elif shape[0]['geometry']['type'] == 'MultiLineString':
                    for linePart in line['geometry']['coordinates']:
                        lineGeom = LineString(linePart)
                        discGeomClip =  self.modelDis['limitGeometry'].intersection(lineGeom)
                        if not discGeomClip.is_empty:
                            discGeomList.append(discGeomClip)
            self.discArrays[name] = {'type':'LineString', 'geomList':discGeomList}
        elif shape[0]['geometry']['type'] == 'MultiLineString':
            for line in shape:
                lineGeom = LineString(line['geometry']['coordinates'])
                discGeomClip =  self.modelDis['limitGeometry'].intersection(lineGeom)
                if not discGeomClip.is_empty:
                    discGeomList.append(discGeomClip)
            self.discArrays[name] = {'type':'LineString', 'geomList':discGeomList}
        elif shape.schema['geometry'] == 'Point':
            for point in shape:
                pointGeom = Point(point['geometry']['coordinates'])
                discGeomClip =  self.modelDis['limitGeometry'].intersection(pointGeom)
                if not discGeomClip.is_empty:
                    discGeomList.append(discGeomClip)
            self.discArrays[name] = {'type':'Point', 'geomList':discGeomList}
        else:
            print('You are working with a uncompatible geometry. Remember to use single parts')
            print('Check this file: %s \n'%shapePath)
            exit()

    def orgVertexAsList(self, geomDict):
        #get only the original vertices inside the model limit
        vertexList = []
        if geomDict['type']=='Polygon':
            for poly in geomDict['geomList']:
                pointObject = poly.exterior.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                for index, point in enumerate(pointList):
                    vertexList.append(point)
        elif geomDict['type']=='LineString':
            if geomDict['geomList'][0].geom_type == "LineString":
                for line in geomDict['geomList']:
                    pointObject = line.coords.xy
                    pointList = list(zip(pointObject[0],pointObject[1]))
                    for index, point in enumerate(pointList):
                        vertexList.append(point)
            elif geomDict['geomList'][0].geom_type == "MultiLineString":
                for part in geomDict['geomList']:
                    for line in part:
                        pointObject = line.coords.xy
                        pointList = list(zip(pointObject[0],pointObject[1]))
                        for index, point in enumerate(pointList):
                            vertexList.append(point)
        elif geomDict['type']=='Point':
            for point in geomDict['geomList']:
                pointObject = point.coords.xy
                vertexList.append((pointObject[0][0],pointObject[1][0]))
        else:
            pass
        return vertexList

    def distributedVertexAsList(self, geomDict):
        #distribute vertices along the layer paths
        vertexList = []
        minRef = self.modelDis['refSizeList'].min()
        if geomDict['type']=='Polygon':
            for poly in geomDict['geomList']:
                polyLength = poly.exterior.length
                pointProg = np.arange(0,polyLength,minRef)
                for prog in pointProg:
                    pointXY = list(poly.exterior.interpolate(prog).xy)
                    vertexList.append([pointXY[0][0],pointXY[1][0]])
        elif geomDict['type']=='LineString':
            if geomDict['geomList'][0].geom_type == "LineString":
                for line in geomDict['geomList']:
                    lineLength = line.length
                    pointProg = np.arange(0,lineLength,minRef)
                    for prog in pointProg:
                        pointXY = list(line.interpolate(prog).xy)
                        vertexList.append([pointXY[0][0],pointXY[1][0]])
            elif geomDict['geomList'][0].geom_type == "MultiLineString":
                for part in geomDict['geomList']:
                    for line in part:
                        lineLength = line.length
                        pointProg = np.arange(0,lineLength,minRef)
                        for prog in pointProg:
                            pointXY = list(line.interpolate(prog).xy)
                            vertexList.append([pointXY[0][0],pointXY[1][0]])
        elif geomDict['type']=='Point':
            for point in geomDict['geomList']:
                pointObject = point.coords.xy
                vertexList.append((pointObject[0][0],pointObject[1][0]))
        else:
            pass
        return vertexList

    def extractOrgVertices(self, txtFile=''):
        start = time.time()
        vertexOrgPairList = []
        vertexDistPairList = []
        for key, dictGeom in self.discArrays.items():
            vertexOrgPairList += self.orgVertexAsList(dictGeom)
            vertexDistPairList += self.distributedVertexAsList(dictGeom)
        self.modelDis['vertexOrg'] = vertexOrgPairList
        self.modelDis['vertexDist'] = vertexDistPairList

        self.modelDis['vertexBuffer'] = []
        if txtFile != '':
            np.savetxt(txtFile,self.modelDis['vertexOrg'])

    def generateCirclesAroundRef(self,indexRef,refSize):
        circleList = []
        for point in self.modelDis['vertexDist']:
            circle = Point(point).buffer(refSize)
            circleList.append(circle)
        circleUnion = unary_union(circleList)
        polyPointList = []
        interiorPolygonList = []

        #check type
        if circleUnion.geom_type == 'MultiPolygon':
            circleMulti = circleUnion
        elif circleUnion.geom_type == 'Polygon':
            circleMulti = MultiPolygon([circleUnion])

        for poly in circleMulti.geoms:

            for interior in poly.interiors:
                interiorPolygon =Polygon(interior.coords)
                interiorPolygonList.append(interiorPolygon)

            outerLength = poly.exterior.length
            if indexRef%2 == 0:
                pointProg = np.arange(0,outerLength,np.pi*refSize/3)
            else:
                pointProg = np.arange(np.pi*refSize/6,outerLength+np.pi*refSize/6,np.pi*refSize/3)
            for prog in pointProg:
                pointXY = list(poly.exterior.interpolate(prog).xy)
                polyPointList.append([pointXY[0][0],pointXY[1][0]])

        interiorMulti= MultiPolygon(interiorPolygonList)
        return circleMulti, interiorPolygonList, polyPointList

    def generateAllCircles(self):

        interiorMultiList = []
        for index, ref in enumerate(self.modelDis['refSizeList']):
            circleUnion, interiorPolygonList, polyPointList = self.generateCirclesAroundRef(index,ref)
            #refBuffer = gpd.GeoSeries(circleUnion)
            self.modelDis['vertexBuffer'] += polyPointList
            interiorMultiList+=interiorPolygonList
            if ref == self.modelDis['refSizeList'].max():
                self.modelDis['circleUnion'] = circleUnion
        self.modelDis['interiors'] = unary_union(MultiPolygon(interiorMultiList))

    def getPointsMinMaxRef(self):

        #define refs
        maxRef = self.modelDis['refSizeList'].max()
        minRef = self.modelDis['refSizeList'].min()

        #define objects to store the uniform vertex
        self.modelDis['vertexMaxRef'] =[]
        self.modelDis['vertexMinRef'] =[]

        #for max ref points
        outerPoly = self.modelDis['limitGeometry']
        innerPolys = self.modelDis['circleUnion']
        for poly in innerPolys.geoms:
            initialPoly = copy.copy(outerPoly)
            transPoly = outerPoly.difference(poly)
            if initialPoly.area == transPoly.area:
                outerPoly.geom.interior += poly
            elif initialPoly.area > transPoly.area:
                outerPoly = outerPoly.difference(poly)
        maxRefXList = np.arange(self.modelDis['xMin']+minRef,self.modelDis['xMax'],maxRef)
        maxRefYList = np.arange(self.modelDis['yMin']+minRef,self.modelDis['yMax'],maxRef)

        for xCoord in maxRefXList:
            for yCoord in maxRefYList:
                refPoint = Point(xCoord,yCoord)
                if outerPoly.contains(refPoint):
                    self.modelDis['vertexMaxRef'].append((xCoord,yCoord))

        self.modelDis['pointsMaxRefPoly']=outerPoly

        #for min ref points
        for polyDict in self.discArrays:
            tempDict = self.discArrays[polyDict]
            if tempDict['type']=='Polygon':
                for poly in tempDict['geomList']:
                    bounds = poly.exterior.bounds
                    minRefXList = np.arange(bounds[0]+minRef,bounds[2],minRef)
                    minRefYList = np.arange(bounds[1]+minRef,bounds[3],minRef)

                    for xCoord in minRefXList:
                        for yCoord in minRefYList:
                            refPoint = Point(xCoord,yCoord)
                            if poly.contains(refPoint):
                                self.modelDis['vertexMinRef'].append((xCoord,yCoord))


    def createPointCloud(self):
        start = time.time()
        #Generate all circles and points on circle paths
        self.generateAllCircles()
        #Distribute points over the max and min refinement areas
        self.getPointsMinMaxRef()
        #Compile all points
        totalRawPoints = []
        totalRawPoints += self.modelDis['vertexDist']
        totalRawPoints += self.modelDis['vertexBuffer']
        totalRawPoints += self.modelDis['vertexMaxRef']
        totalRawPoints += self.modelDis['vertexMinRef']
        totalDefPoints = []

        #check if points are inside limit polygon
        for point in totalRawPoints:
            refPoint = Point(point[0],point[1])
            if self.modelDis['limitGeometry'].contains(refPoint):
                totalDefPoints.append(point)
        self.modelDis['vertexTotal'] = totalDefPoints

        line1 = '/---Point cloud---/'
        line2 = 'From layers: %d'%len(self.modelDis['vertexDist'])
        line3 = 'From layer buffers: %d'%len(self.modelDis['vertexBuffer'])
        line4 = 'From max. ref. areas: %d'%len(self.modelDis['vertexMaxRef'])
        line5 = 'From min. ref. areas: %d'%len(self.modelDis['vertexMinRef'])
        line6 = 'Points inside limit: %d'%len(self.modelDis['vertexTotal'])
        end = time.time()
        line7 = 'Time required: %.2f s.'%(end - start)
        outMessage = line1 + '<br/>' + line2 + '<br/>' + line3 + '<br/>' + line4 + '<br/>' + line5 + '<br/>' + line6 + '<br/>' + line7
        return outMessage

    def generateVoronoi(self):
        #create a multipoint object
        pointMulti = MultiPoint(self.modelDis['vertexTotal'])
        #original regions
        regions = voronoi_diagram(pointMulti)
        #object for clipped regions
        clippedRegions = []
        for region in regions.geoms:
            if self.modelDis['limitGeometry'].contains(region):
                clippedRegions.append(region)
            else:
                regionDiff = region.intersection(self.modelDis['limitGeometry'])
                clippedRegions.append(regionDiff)
        for region in clippedRegions:
            if region.area <= 1:
                print(region.area)
        clippedRegionsMulti = MultiPolygon(clippedRegions)
        self.modelDis['voronoiRegions'] = clippedRegionsMulti

    def getPolyAsShp(self,circleList,shapePath=''):
        schema_props = OrderedDict([("id", "int")])
        schema={"geometry": "Polygon", "properties": schema_props}
        if shapePath != '':
            outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                                crs = self.modelDis['crs'], schema=schema)
            if self.modelDis[circleList].geom_type == 'Polygon':
                poly = self.modelDis[circleList]
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
                    "properties": OrderedDict([("id",1)]),
                }
                outFile.write(feature)
            elif self.modelDis[circleList].geom_type == 'MultiPolygon':
                for index, poly in enumerate(self.modelDis[circleList].geoms):
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

    def getPointsAsShp(self,pointList,shapePath=''):
        schema_props = OrderedDict([("id", "int")])
        schema={"geometry": "Point", "properties": schema_props}
        if shapePath != '':
            outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                                crs = self.modelDis['crs'], schema=schema)
            for index, point in enumerate(self.modelDis[pointList]):
                feature = {
                    "geometry": {'type':'Point',
                                'coordinates':(point[0],point[1])},
                    "properties": OrderedDict([("id",index)]),
                }
                outFile.write(feature)
            outFile.close()