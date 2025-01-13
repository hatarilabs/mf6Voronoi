import numpy as np
import copy, sys
import tqdm, time
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi,cKDTree
#import geospatial libraries
import fiona
from tqdm import tqdm
from shapely.ops import split, unary_union, cascaded_union, voronoi_diagram
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from collections import OrderedDict

class createVoronoi():
    def __init__(self):
        self.discGeoms = {}
        self.modelDis = {}
        self.pairArray = None
        self.discLayers = {}

    def isMultiGeometry(self, geom):
        return isinstance(geom, (MultiPoint, MultiLineString, MultiPolygon))

    def addLimit(self, name, shapePath):
        #Create the model limit
        limitShape = fiona.open(shapePath)

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
        self.modelDis['vertexDist'] = {}
        self.modelDis['crs'] = limitShape.crs

    # we shift from minRef to layerRef
    # def defineParameters(self, maxRef=500, minRef=50, multiplier=1):
    def defineParameters(self, maxRef=500, multiplier=1):
        #Define the refinement sizes
        # distList = [minRef]
        # i=1
        # while distList[-1] <= maxRef:
        #     #print(distList)
        #     distValue = distList[-1] + multiplier**i*minRef
        #     if distValue <= maxRef:
        #         distList.append(distValue)       
        #     else:
        #         break
        #     i+=1

        # self.modelDis['refSizeList'] = np.array(distList)
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
        #self.modelDis['minRef'] = minRef
        #self.modelDis['refSizeList'] = [minRef + index*]
        # print('\n/--------Sumary of cell discretization-------/')
        # print('Maximun refinement progressive: %.2f m.'%self.modelDis['refSizeList'].max())
        # print('Maximun refinement coarse areas: %.2f m.'%self.modelDis['maxRef'])
        # print('Minimum refinement: %.2f m.'%self.modelDis['minRef'])
        # print('Cell size list: %s m.'%str(self.modelDis['refSizeList']))
        # print('/--------------------------------------------/\n',flush=True)

    #here we add the layerRef to the function
    #def addLayer(self, name, shapePath):
    def addLayer(self, name, shapePath, layerRef):
        #Add layers for mesh definition
        #This feature also clips and store the geometry
        #geomList is allways a Python of Shapely geometries
        
        spatialDf = gpd.read_file(shapePath)   

        #get the master layer name
        self.discLayers[name] = {'layerRef':layerRef}  
        
        #auxiliary funtion to intersect:
        def intersectGeom(anyGeom, partType, multiPartType):
            discGeomList = []
            #generic 
            
            if self.isMultiGeometry(anyGeom):
                for partGeom in anyGeom.geoms:
                    discGeomClip =  self.modelDis['limitGeometry'].intersection(partGeom)
                    if not discGeomClip.is_empty:
                        discGeomList.append(discGeomClip)
            else:
                discGeomClip =  self.modelDis['limitGeometry'].intersection(anyGeom)
                if not discGeomClip.is_empty:
                    discGeomList.append(discGeomClip)

            unaryPoly = unary_union(discGeomList)
            if unaryPoly.geom_type == partType:
                unaryFilter = [unaryPoly]
            elif unaryPoly.geom_type == multiPartType:
                unaryFilter = [poly for poly in unaryPoly.geoms]
            return unaryFilter

        #looping over the shapefile
        for spatialIndex, spatialRow in spatialDf.iterrows():
            #working for polygons and multipolygons
            if spatialRow.geometry.geom_type == 'Polygon':
                polyGeom = spatialRow.geometry
                unaryFilter = intersectGeom(polyGeom,'Polygon','MultiPolygon')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'Polygon', 
                                                              'geomList':unaryFilter
                                                              'layerName':name}
                
            elif spatialRow.geometry.geom_type == 'MultiPolygon':
                multiPolyGeom = spatialRow.geometry
                unaryFilter = intersectGeom(multiPolyGeom,'Polygon','MultiPolygon')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'Polygon', 
                                                              'geomList':unaryFilter
                                                              'layerName':name}
                
            #working for lines and multilines    
            elif spatialRow.geometry.geom_type == 'LineString':
                lineGeom = spatialRow.geometry
                unaryFilter = intersectGeom(lineGeom,'LineString','MultiLineString')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'LineString', 
                                                              'geomList':unaryFilter,
                                                              'layerName':name}

            elif spatialRow.geometry.geom_type == 'MultiLineString':
                multiLineGeom = spatialRow.geometry
                unaryFilter = intersectGeom(multiLineGeom,'LineString','MultiLineString')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'LineString', 
                                                              'geomList':unaryFilter
                                                              'layerName':name}

            #working for points and multipoints
            elif spatialRow.geometry.geom_type == 'Point':
                pointGeom = spatialRow.geometry
                unaryFilter = intersectGeom(pointGeom,'Point','MultiPoint')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'Point', 
                                                              'geomList':unaryFilter,
                                                              'layerRef':layerRef,
                                                              'layerName':name}

            elif spatialRow.geometry.geom_type == 'MultiPoint':
                pointGeom = spatialRow.geometry
                unaryFilter = intersectGeom(pointGeom,'Point','MultiPoint')
                self.discGeoms[name+'_'+str(spatialIndex)] = {'type':'Point', 
                                                              'geomList':unaryFilter,
                                                              'layerRef':layerRef,
                                                              'layerName':name}
            else:
                print('You are working with a uncompatible geometry. Remember to use single parts')
                print('Check this file: %s \n'%shapePath)
                sys.exit()

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
            for line in geomDict['geomList']:
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
        layerRef = geomDict['layerRef']
        if geomDict['type']=='Polygon':
            for poly in geomDict['geomList']:
                polyLength = poly.exterior.length
                pointProg = np.arange(0,polyLength,layerRef)
                for prog in pointProg:
                    pointXY = list(poly.exterior.interpolate(prog).xy)
                    vertexList.append([pointXY[0][0],pointXY[1][0]])
        elif geomDict['type']=='LineString':
            for line in geomDict['geomList']:
                lineLength = line.length
                pointProg = np.arange(0,lineLength,layerRef)
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

    def extractOrgVertices(self, txtFile='', probIndex=0.01):
        start = time.time()
        vertexOrgPairList = []
        #vertexDistPairList = []
        for key, dictGeom in self.discGeoms.items():
            vertexOrgPairList += self.orgVertexAsList(dictGeom)
            #vertexDistPairList += self.distributedVertexAsList(dictGeom)
            self.modelDis['vertexDist'][key] = self.distributedVertexAsList(dictGeom)
        self.modelDis['vertexOrg'] = vertexOrgPairList

        self.modelDis['vertexBuffer'] = []
        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,key,indexRef,refSize):
        #first we create buffers around points and merge them
        circleList = []
        polyPointList = []
        for point in self.modelDis['vertexDist'][key]:
            circle = Point(point).buffer(refSize)
            circleList.append(circle)
        circleUnion = unary_union(circleList)
        polyPointList = []
        if circleUnion.geom_type == 'MultiPolygon':
            circleMulti = circleUnion
        elif circleUnion.geom_type == 'Polygon':
            circleMulti = MultiPolygon([circleUnion])
        # from the multipolygons 
        for poly in circleMulti.geoms:
            outerLength = poly.exterior.length
            if indexRef%2 == 0:
                pointProg = np.arange(0,outerLength,np.pi*refSize/3)
            else:
                pointProg = np.arange(np.pi*refSize/6,outerLength+np.pi*refSize/6,np.pi*refSize/3)
            for prog in pointProg:
                pointXY = list(poly.exterior.interpolate(prog).xy)
                polyPointList.append([pointXY[0][0],pointXY[1][0]])
        return circleMulti, polyPointList

    def generateAllCircles(self):
        partialCircleUnionList = []

        label = ''

        for layer in self.discLayers:
            distList = [layer['layerRef']]

            i=1
            while distList[-1] <= self.modelDis['maxRef']:
                distValue = distList[-1] + self.modelDis['multiplier']**i*value['layerRef']
                if distValue <= self.modelDis['maxRef']:
                    distList.append(distValue)       
                else:
                    break
                i+=1

            self.discLayers[layer]['distList'] = distList

            print('\n/--------Layer %s discretization-------/'%layer)
            print('Cell size list: %s m.'%str(distList))
            
            #filter discGeom of this layer
            layerDiscGeoms = [ self.discGeoms[x] for x in self.discGeoms if x['layerName'] == layer]

            for layerDiscGeom in layerDiscGeoms:

        for key, value in self.discGeoms.items():
            distList = [value['layerRef']]

            i=1
            while distList[-1] <= self.modelDis['maxRef']:
                #print(distList)
                distValue = distList[-1] + self.modelDis['multiplier']**i*value['layerRef']
                if distValue <= self.modelDis['maxRef']:
                    distList.append(distValue)       
                else:
                    break
                i+=1

            self.discGeoms[key]['distList'] = distList

            if label != key.split('_')[0]:
                label = key.split('_')[0]
                print('\n/--------Layer %s discretization-------/'%label)
                print('Cell size list: %s m.'%str(distList))
                #print('/--------------------------------------------/\n',flush=True)        

            #partialCircleUnionList = []
            for index, ref in enumerate(distList):
                #print(key,index,ref)
                circleUnion, polyPointList = self.circlesAroundRefPoints(key,index,ref)
                refBuffer = gpd.GeoSeries(circleUnion)
                self.modelDis['vertexBuffer'] += polyPointList
                #here we use the maximum progressive refinement
                #if ref == self.modelDis['refSizeList'].max():
                if ref == np.array(distList).max():
                    #self.modelDis['circleUnion'] = circleUnion
                    partialCircleUnionList.append(circleUnion)

        totalCircleUnion = unary_union(partialCircleUnionList)
        self.modelDis['circleUnion'] = totalCircleUnion


        # for index, ref in enumerate(self.modelDis['refSizeList']):
        #     circleUnion, polyPointList = self.circlesAroundRefPoints(index,ref)
        #     refBuffer = gpd.GeoSeries(circleUnion)
        #     self.modelDis['vertexBuffer'] += polyPointList
        #     #here we use the maximum progressive refinement
        #     if ref == self.modelDis['refSizeList'].max():
        #         self.modelDis['circleUnion'] = circleUnion

    def getPointsMinMaxRef(self):

        #define refs
        maxRef = self.modelDis['maxRef']

        refList = []
        for key, value in self.discGeoms.items():
            refList += value['distList']

        #minRef = self.modelDis['minRef']
        minRef = np.array(refList).min()

        #define objects to store the uniform vertex
        self.modelDis['vertexMaxRef'] =[]
        self.modelDis['vertexMinRef'] =[]

        #get the limit geometry where no coarse grid will be generated
        outerPoly = self.modelDis['limitGeometry']
        limitPoly = copy.copy(outerPoly)
        innerPolys = self.modelDis['circleUnion']

        #working with circle unions

        if self.isMultiGeometry(innerPolys):
            for poly in innerPolys.geoms:
                transPoly = outerPoly.difference(poly)
                if limitPoly.area == transPoly.area:
                    outerPoly.geom.interior += poly
                elif limitPoly.area > transPoly.area:
                    outerPoly = outerPoly.difference(poly)
        else:
            transPoly = outerPoly.difference(innerPolys)
            if limitPoly.area == transPoly.area:
                outerPoly.geom.interior += transPoly
            elif limitPoly.area > transPoly.area:
                outerPoly = outerPoly.difference(transPoly)

        #working with mesh disc polys
        for key, value in self.discGeoms.items():
            if value['type'] == 'Polygon':
                for poly in value['geomList']:
                    transPoly = outerPoly.difference(poly)
                    if limitPoly.area == transPoly.area:
                        outerPoly.geom.interior += poly
                    elif limitPoly.area > transPoly.area:
                        outerPoly = outerPoly.difference(poly)
                                 
        #exporting final clipped polygon geometry                         
        self.modelDis['pointsMaxRefPoly']=outerPoly

        #creating points of coarse grid
        maxRefXList = np.arange(self.modelDis['xMin']+minRef,self.modelDis['xMax'],maxRef)
        maxRefYList = np.arange(self.modelDis['yMin']+minRef,self.modelDis['yMax'],maxRef)

        for xCoord in maxRefXList:
            for yCoord in maxRefYList:
                refPoint = Point(xCoord,yCoord)
                if not outerPoly.contains(refPoint):
                    self.modelDis['vertexMaxRef'].append((xCoord,yCoord))

        self.modelDis['pointsMaxRefPoly']=outerPoly

        #for min ref points
        for polyKey in self.discGeoms:
            polyDict = self.discGeoms[polyKey]
            if polyDict['type']=='Polygon':
                for poly in polyDict['geomList']:
                    bounds = poly.exterior.bounds
                    minRefXList = np.arange(bounds[0]+polyDict['layerRef'],bounds[2],polyDict['layerRef'])
                    minRefYList = np.arange(bounds[1]+polyDict['layerRef'],bounds[3],polyDict['layerRef'])

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
        #totalRawPoints += self.modelDis['vertexDist']
        for key in self.modelDis['vertexDist']:
            totalRawPoints += self.modelDis['vertexDist'][key]
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

        print('\n/----Sumary of points for voronoi meshing----/')
        print('Distributed points from layers: %d'%len(self.modelDis['vertexDist']))
        print('Points from layer buffers: %d'%len(self.modelDis['vertexBuffer']))
        print('Points from max refinement areas: %d'%len(self.modelDis['vertexMaxRef']))
        print('Points from min refinement areas: %d'%len(self.modelDis['vertexMinRef']))
        print('Total points inside the limit: %d'%len(self.modelDis['vertexTotal']))
        print('/--------------------------------------------/')
        end = time.time()
        print('\nTime required for point generation: %.2f seconds \n'%(end - start), flush=True)

    def generateVoronoi(self):
        print('\n/----Generation of the voronoi mesh----/')
        start = time.time()
        #create a multipoint object
        pointMulti = MultiPoint(self.modelDis['vertexTotal'])
        #original regions
        regions = voronoi_diagram(pointMulti)
        #object for clipped regions
        clippedRegions = []
        #loop over all polygons
        for region in regions.geoms:
            #for contained polygons
            if self.modelDis['limitGeometry'].contains(region):
                clippedRegions.append(region)
            #for intersected polygons
            else:
                regionDiff = region.intersection(self.modelDis['limitGeometry'])
                #check for clipped region as multipolygon
                if regionDiff.geom_type == 'Polygon':
                    clippedRegions.append(regionDiff)
                elif regionDiff.geom_type == 'MultiPolygon':
                    clippedRegions.extend(list(regionDiff.geoms))
                else: print('Something went wrong')

        clippedRegionsMulti = MultiPolygon(clippedRegions)
        self.modelDis['voronoiRegions'] = clippedRegionsMulti
        end = time.time()
        print('\nTime required for voronoi generation: %.2f seconds \n'%(end - start), flush=True)

    def getPolyAsShp(self,circleList,shapePath=''):
        print('\n/----Generation of the voronoi shapefile----/')
        start = time.time()
        schema_props = OrderedDict([("id", "int")])
        schema={"geometry": "Polygon", "properties": schema_props}
        if shapePath != '':
            outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                                crs = self.modelDis['crs'], schema=schema)
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
        end = time.time()
        print('\nTime required for voronoi shapefile: %.2f seconds \n'%(end - start), flush=True)


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
