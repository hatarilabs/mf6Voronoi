import numpy as np
import copy, sys, os
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
from .utils import (filterVertexCloseLimit, 
                    getFionaDictPoly, 
                    getFionaDictPoint, 
                    intersectLimitLayer, 
                    isMultiGeometry)

class createVoronoi():
    def __init__(self, meshName, maxRef, multiplier):
        #self.discGeoms = {}
        self.modelDis = {}
        self.modelDis['meshName'] = meshName
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
        self.pairArray = None
        self.discLayers = {}

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
        self.modelDis['vertexDistGeoms'] = {}
        self.modelDis['vertexBuffer'] = []
        self.modelDis['crs'] = limitShape.crs

    #here we add the layerRef to the function
    def addLayer(self, layerName, shapePath, layerRef):
        #Add layers for mesh definition
        #This feature also clips and store the geometry
        #geomList is allways a Python of Shapely geometries
        spatialDf = gpd.read_file(shapePath)   

        #get the ref and geoms as a list
        self.discLayers[layerName] = {'layerRef':layerRef,
                                      'layerGeoms':[]}  

        #looping over the shapefile
        i = 1
        for spatialIndex, spatialRow in spatialDf.iterrows():
            if spatialRow.geometry.is_valid:
                geomGeom = spatialRow.geometry
                #get the layer type
                if i==1:
                    self.discLayers[layerName]['layerType'] = geomGeom.geom_type
                    i+=1
                #intersect with the limit layer
                unaryFilter = intersectLimitLayer(geomGeom, self.modelDis)
                self.discLayers[layerName]['layerGeoms'] += unaryFilter
            else:
                print('You are working with a uncompatible geometry. Remember to use single parts')
                print('Check this file: %s \n'%shapePath)
                sys.exit()
            
    #def orgVertexAsList(self, layerGeoms, layerRef):
    def orgVertexAsList(self, layer):
        #get only the original vertices inside the model limit
        vertexList = []
        layerGeoms = self.discLayers[layer]['layerGeoms']
        layerRef = self.discLayers[layer]['layerRef']

        for layerGeom in layerGeoms:
            if layerGeom.geom_type == 'Polygon':
                # pointObject = layerGeom.exterior.coords.xy
                # pointList = list(zip(pointObject[0],pointObject[1]))
                filterPointList = filterVertexCloseLimit(layerRef,layerGeom,self.modelDis)
                if filterPointList != None:
                    vertexList += filterPointList
            elif layerGeom.geom_type == 'LineString':
                # pointObject = layerGeom.coords.xy
                # pointList = list(zip(pointObject[0],pointObject[1]))
                filterPointList  = filterVertexCloseLimit(layerGeom,pointList,self.modelDis)
                if filterPointList != None:
                    vertexList += filterPointList
            elif layerGeom.geom_type == 'Point':
                pointObject = layerGeom.coords.xy
                point = (pointObject[0][0],pointObject[1][0])
                pointPoint = Point(point)
                if pointPoint.buffer(self.modelDis['maxRef']).within(self.modelDis['limitGeometry']):
                    vertexList.append((pointObject[0][0],pointObject[1][0]))
            else:
                print(layerGeom)
                print('/-----Problem has been bound when extracting org vertex-----/')

        return vertexList

    def distributedVertexAsList(self, layerGeoms, layerRef):
        #distribute vertices along the layer paths
        vertexList = []
        vertexGeomList = []

        for layerGeom in layerGeoms:
            if layerGeom.geom_type == 'Polygon':
                polyLength = layerGeom.exterior.length
                pointProg = np.arange(0,polyLength,layerRef)
                pointList = []
                for prog in pointProg:
                    pointXY = list(layerGeom.exterior.interpolate(prog).xy)
                    pointList.append([pointXY[0][0],pointXY[1][0]])
                filterPointList, filterPointGeom = filterVertexCloseLimit(layerGeom,pointList,self.modelDis)
                if filterPointGeom != None:
                    vertexList += filterPointList
                    vertexGeomList.append(filterPointGeom)
            elif layerGeom.geom_type == 'LineString':
                lineLength = layerGeom.length
                pointProg = np.arange(0,lineLength,layerRef)
                pointList = []
                for prog in pointProg:
                    pointXY = list(layerGeom.interpolate(prog).xy)
                    pointList.append([pointXY[0][0],pointXY[1][0]])
                filterPointList, filterPointGeom = filterVertexCloseLimit(layerGeom,pointList,self.modelDis)
                if filterPointGeom != None:
                    vertexList += filterPointList
                    vertexGeomList.append(filterPointGeom)
            elif layerGeom.geom_type == 'Point':
                pointObject = layerGeom.coords.xy
                point = (pointObject[0][0],pointObject[1][0])
                pointPoint = Point(point)
                if pointPoint.buffer(self.modelDis['maxRef']).within(self.modelDis['limitGeometry']):
                    x,y = pointObject[0][0],pointObject[1][0]
                    vertexList.append((x,y))
                    vertexGeomList.append(pointPoint)
            else:
                print('/-----Problem has been bound when extracting dist vertex-----/')

            # # to prevent vertices close to borders
            # for vertex in vertexList:
            #     vertexPoint = Point(vertex)
            #     if vertexPoint.buffer*
            vertexLine = LineString(vertexList)
            # if len(vertexList) > 1:
            #     bufferVertexList = []
            #     vertexLine = LineString(vertexList)
            #     vertexProg = np.arange(0,vertexLine.length,layerRef/4)
            #     for prog in vertexProg:
            #         pointXY = list(vertexLine.interpolate(prog).xy)
            #         bufferVertexList.append([pointXY[0][0],pointXY[1][0]])
            # else:
            #     bufferVertexList = vertexList

        return vertexList, vertexGeomList

    def generateOrgDistVertices(self, txtFile=''):
        start = time.time()
        vertexOrgPairList = []
        #vertexDistPairList = []
        for layer, values in self.discLayers.items():
            layerGeoms = values['layerGeoms']
            layerRef = values['layerRef']
            vertexOrgPairList += self.orgVertexAsList(values['layerGeoms'], values['layerRef'])
            self.modelDis['vertexDist'][layer] = self.distributedVertexAsList(layerGeoms, layerRef)[0]
            self.modelDis['vertexDistGeoms'][layer] = self.distributedVertexAsList(layerGeoms, layerRef)[1]
        self.modelDis['vertexOrg'] = vertexOrgPairList

        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,layer,indexRef,cellSize):
        
        #first we create buffers around points and merge them
        circleList = []
        polyPointList = []
        # for point in self.modelDis['vertexDistBuffer'][layer]:
        #     circle = Point(point).buffer(cellSize)
        #     circleList.append(circle)
        # circleUnions = unary_union(circleList)
        #get the layer space dist
        layerSpaceList = self.discLayers[layer]['layerSpaceList']
        layerSpaceFraction = layerSpaceList.index(cellSize)/len(layerSpaceList)
        firstCellSize = layerSpaceList[0]


        for geom in self.modelDis['vertexDistGeoms'][layer]:
            #fixing for the first cell avoiding long cells
            circle = geom.buffer(cellSize - firstCellSize/2)
            circleList.append(circle)
        circleUnions = unary_union(circleList)
        #to avoid corner
        #circleUnions =  circleUnions.simplify(tolerance=4, preserve_topology=True)
        
        def getPolygonAndInteriors(polyGeom):
            exteriorInteriorPolys = [polyGeom] + [Polygon(ring) for ring in polyGeom.interiors]
            return exteriorInteriorPolys
         
        circleUnionExtIntList = []
        circleUnionExtWithIntList = []
        if circleUnions.geom_type == 'MultiPolygon':
            for circleUnion in circleUnions.geoms:
                circleUnionExtIntList += getPolygonAndInteriors(circleUnion)
                circleUnionExtWithIntList.append(circleUnion)
        elif circleUnions.geom_type == 'Polygon':
            circleUnionExtIntList += getPolygonAndInteriors(circleUnions)
            circleUnionExtWithIntList.append(circleUnions)

        #define opening angle 
        opAngle = 30*np.pi/180
        

        # from the multipolygons 
        polyPointList = []
        for circleUnionExtInt in circleUnionExtIntList:
            outerLength = circleUnionExtInt.exterior.length
            #pointProg = np.arange(0,outerLength,np.sin(np.pi/2 - layerSpaceFraction*np.pi/6)*cellSize)
            pointProg = np.arange(0,outerLength,(0.8 - layerSpaceFraction*0.4)*cellSize)
            # if indexRef%2 == 0:
            #     pointProg = np.arange(0,outerLength,np.sin(opAngle)*cellSize)
            # else:
            #     pointProg = np.arange(np.pi*cellSize/6,outerLength+np.pi*cellSize/6,np.sin(opAngle)*cellSize)
            #     #pointProg = np.arange(np.pi*cellSize/6,outerLength+np.pi*cellSize/6,np.pi*cellSize/3)
            for prog in pointProg:
                pointXY = list(circleUnionExtInt.exterior.interpolate(prog).xy)
                polyPointList.append([pointXY[0][0],pointXY[1][0]])

        circleUnionExtIntMpoly = MultiPolygon(circleUnionExtIntList)
        circleUnionExtWithIntMpoly = MultiPolygon(circleUnionExtWithIntList)

        lastGeometry = self.modelDis['activeArea'][-1]
        partialActiveArea = lastGeometry.difference(circleUnionExtWithIntMpoly)
        self.modelDis['activeArea'].append(partialActiveArea)

        #print(circleUnionExtWithIntList)
        return circleUnionExtWithIntMpoly, circleUnionExtIntMpoly, polyPointList

    def generateAllCircles(self):
        partialCircleUnionList = []
        partialCircleUnionInteriorList = []

        label = ''

        for layer, value in self.discLayers.items():
            cellSizeList = [value['layerRef']]

            i=1
            while cellSizeList[-1] <= self.modelDis['maxRef']:
                cellSize = cellSizeList[-1] + self.modelDis['multiplier']**i*value['layerRef']
                if cellSize <= self.modelDis['maxRef']:
                    cellSizeList.append(cellSize)       
                else:
                    break
                i+=1

            self.discLayers[layer]['layerSpaceList'] = cellSizeList

            print('\n/--------Layer %s discretization-------/'%layer)
            print('Progressive cell size list: %s m.'%str(cellSizeList))

            #starting from the second size if multiple sizes are present
            # if len(cellSizeList) > 1:
            #     cellSizeList = cellSizeList[1:]
            # else:
            #     cellSizeList = cellSizeList

            #fixing distances for the first cell
            #cellSizeList = [size - cellSizeList[0]/2 for size in cellSizeList]
            
            #initiate active area list:
            self.modelDis['activeArea'] = [self.modelDis['limitGeometry']]

            #looping
            for index, cellSize in enumerate(cellSizeList):
                circleUnionInteriors, circleUnion, polyPointList = self.circlesAroundRefPoints(layer,index,cellSize)
                refBuffer = gpd.GeoSeries(circleUnion)
                self.modelDis['vertexBuffer'] += polyPointList
                #here we use the maximum progressive refinement
                #if ref == self.modelDis['refSizeList'].max():
                #for the last discretization
                if cellSize == np.array(cellSizeList).max():
                    #self.modelDis['circleUnion'] = circleUnion
                    partialCircleUnionList.append(circleUnion)
                    partialCircleUnionInteriorList.append(circleUnionInteriors)

        totalCircleUnion = unary_union(partialCircleUnionList)
        totalCircleUnionInteriors = unary_union(partialCircleUnionInteriorList)

        self.modelDis['circleUnion'] = totalCircleUnion
        self.modelDis['circleUnionInteriors'] = totalCircleUnionInteriors

    def getPointsMinMaxRef(self):

        #define refs
        maxRef = self.modelDis['maxRef']

        layerRefList = []
        for key, value in self.discLayers.items():
            layerRefList.append(value['layerRef'])

        #minRef = self.modelDis['minRef']
        minRef = np.array(layerRefList).min()

        #define objects to store the uniform vertex
        self.modelDis['vertexMaxRef'] =[]
        self.modelDis['vertexMinRef'] =[]

        #get the limit geometry where no coarse grid will be generated
        outerPoly = self.modelDis['limitGeometry']
        limitPoly = copy.copy(outerPoly)
        innerPolys = self.modelDis['circleUnionInteriors']

        #working with circle unions
        if isMultiGeometry(innerPolys):
            for poly in innerPolys.geoms:
                transPoly = outerPoly.difference(poly)
                if limitPoly.area == transPoly.area:
                    outerPoly.geom.interior += poly
                elif limitPoly.area > transPoly.area:
                    outerPoly = transPoly
        else:
            transPoly = outerPoly.difference(innerPolys)
            self.modelDis['tempPoly']=transPoly
            if limitPoly.area == transPoly.area:
                outerPoly.geom.interior += transPoly
            elif limitPoly.area > transPoly.area:
                outerPoly = transPoly

        #working with mesh disc polys
        for key, value in self.discLayers.items():
            for layerGeom in value['layerGeoms']:
                if layerGeom.geom_type == 'Polygon':
                    transPoly = outerPoly.difference(layerGeom)
                    if limitPoly.area == transPoly.area:
                        outerPoly.geom.interior += layerGeom
                    elif limitPoly.area > transPoly.area:
                        outerPoly = outerPoly.difference(layerGeom)
                                 
        #exporting final clipped polygon geometry                         
        self.modelDis['pointsMaxRefPoly']=outerPoly

        #creating points of coarse grid
        maxRefXList = np.arange(self.modelDis['xMin']+minRef,self.modelDis['xMax'],maxRef)
        maxRefYList = np.arange(self.modelDis['yMin']+minRef,self.modelDis['yMax'],maxRef)

        for xCoord in maxRefXList:
            for yCoord in maxRefYList:
                refPoint = Point(xCoord,yCoord)
                if outerPoly.contains(refPoint):
                    self.modelDis['vertexMaxRef'].append((xCoord,yCoord))

        self.modelDis['pointsMaxRefPoly']=outerPoly

        #for min ref points
        for key, value in self.discLayers.items():
            for layerGeom in value['layerGeoms']:
                if layerGeom.geom_type == 'Polygon':
                    bounds = layerGeom.exterior.bounds
                    minRefXList = np.arange(bounds[0]+value['layerRef'],bounds[2],value['layerRef'])
                    minRefYList = np.arange(bounds[1]+value['layerRef'],bounds[3],value['layerRef'])

                    for xCoord in minRefXList:
                        for yCoord in minRefYList:
                            refPoint = Point(xCoord,yCoord)
                            if layerGeom.contains(refPoint):
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
        regions = voronoi_diagram(pointMulti, tolerance=0.5)
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

    def getVoronoiAsShp(self, shapePath=''):
        print('\n/----Generation of the voronoi shapefile----/')
        start = time.time()
        schema_props = OrderedDict([("id", "int")])
        schema={"geometry": "Polygon", "properties": schema_props}

        #check or create an output folder
        # if os.path.isdir(outputPath):
        #     print('The output folder %s exists'%outputPath)
        # else:
        #     os.mkdir(outputPath)
        #     print('The output folder %s has been generated.'%outputPath)

        # shapePath = os.path.join(outputPath, self.modelDis['meshName']+'.shp')

        outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                            crs = self.modelDis['crs'], schema=schema)
        
        for index, poly in enumerate(self.modelDis['voronoiRegions'].geoms):
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

    def getPolyAsShp(self,circleList,shapePath=''):
        start = time.time()
        schema_props = OrderedDict([("id", "str")])
        schema={"geometry": "Polygon", "properties": schema_props}

        #check or create an output folder
        # if os.path.isdir(outputPath):
        #     print('The output folder %s exists'%outputPath)
        # else:
        #     os.mkdir(outputPath)
        #     print('The output folder %s has been generated.'%outputPath)

        #shapePath = os.path.join(outputPath, self.modelDis['meshName']+circleList)
        
        outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                            crs = self.modelDis['crs'], schema=schema)
        
        if isinstance(self.modelDis[circleList], dict):
            for key, value in self.modelDis[circleList].items():
                if isMultiGeometry(value):
                    for index, poly in enumerate(value.geoms):
                        feature = getFionaDictPoly(poly, index)
                        outFile.write(feature)

        if isMultiGeometry(self.modelDis[circleList]):
            for index, poly in enumerate(self.modelDis[circleList].geoms):
                feature = getFionaDictPoly(poly, index)
                outFile.write(feature)
        else:
            poly = self.modelDis[circleList]
            feature = getFionaDictPoly(poly, '1')
            outFile.write(feature)
        outFile.close()
        
        end = time.time()
        print('\nTime required for voronoi shapefile: %.2f seconds \n'%(end - start), flush=True)

    def getPointsAsShp(self,pointList,shapePath=''):
        schema_props = OrderedDict([("id", "str")])
        schema={"geometry": "Point", "properties": schema_props}
        if shapePath != '':
            outFile = fiona.open(shapePath,mode = 'w',driver = 'ESRI Shapefile',
                                crs = self.modelDis['crs'], schema=schema)
            if isinstance(self.modelDis[pointList], dict):
                print(self.modelDis[pointList].keys())
                for key, value in self.modelDis[pointList].items():
                    for index, point in enumerate(value):
                        feature = getFionaDictPoint(point, index)
                        if feature != None:
                            outFile.write(feature)
                        else:
                            print('Something went wrong with %s'%point)
            else:
                for index, point in enumerate(self.modelDis[pointList]):
                    feature = getFionaDictPoint(point, index)
                    if feature != None:
                        outFile.write(feature)
                    else:
                        print('Something went wrong with %s'%point)
            outFile.close()
