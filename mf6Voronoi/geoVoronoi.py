import numpy as np
import copy, sys, time
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.spatial import Voronoi,cKDTree
#import geospatial libraries
import fiona
from tqdm import tqdm
from shapely.ops import split, unary_union, voronoi_diagram
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from collections import OrderedDict
from .utils import (processVertexFilterCloseLimit, 
                    intersectLimitLayer, 
                    isMultiGeometry,
                    isRunningInJupyter, 
                    printBannerHtml, 
                    printBannerText)

class createVoronoi():
    def __init__(self, meshName, maxRef, multiplier, overlapping=True):
        #self.discGeoms = {}
        self.modelDis = {}
        self.modelDis['meshName'] = meshName
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
        self.overlapping = overlapping
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
        #initiate active area list:
        self.modelDis['activeArea'] = [self.modelDis['limitGeometry']]

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
                if unaryFilter:
                    #if not unaryFilter.is_empty:
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
            filterPointList = processVertexFilterCloseLimit(layerRef,layerGeom,self.modelDis,'Org')
            if filterPointList != None:
                vertexList += filterPointList
            else:
                print('/-----Problem has been bound when extracting org vertex-----/')

        return vertexList

    def distributedVertexAsList(self, layer):
        #distribute vertices along the layer paths
        vertexList = []
        vertexGeomList = []
        layerGeoms = self.discLayers[layer]['layerGeoms']
        layerRef = self.discLayers[layer]['layerRef']

        for layerGeom in layerGeoms:
            filterPointList, filterPointGeom = processVertexFilterCloseLimit(layerRef,layerGeom,self.modelDis,'Dist')
            if filterPointGeom != None:
                vertexList += filterPointList
                vertexGeomList.append(filterPointGeom)
        return vertexList, vertexGeomList

    def generateOrgDistVertices(self, txtFile=''):
        vertexOrgPairList = []
        for layer, values in self.discLayers.items():
            vertexOrgPairList += self.orgVertexAsList(layer)
            self.modelDis['vertexDist'][layer] = self.distributedVertexAsList(layer)[0]
            self.modelDis['vertexDistGeoms'][layer] = self.distributedVertexAsList(layer)[1]
        self.modelDis['vertexOrg'] = vertexOrgPairList

        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,layer,indexRef,cellSize):
        
        #first we create buffers around points and merge them
        circleList = []
        polyPointList = []
        layerSpaceList = self.discLayers[layer]['layerSpaceList']
        layerSpaceFraction = layerSpaceList.index(cellSize)/len(layerSpaceList)
        firstCellSize = layerSpaceList[0]

        for geom in self.modelDis['vertexDistGeoms'][layer]:
            #fixing for the first cell avoiding long cells
            #circle = geom.buffer(cellSize - firstCellSize/2) #Check this
            circle = geom.buffer(cellSize) #Check this
            circleList.append(circle)
        circleUnions = unary_union(circleList)

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
        
        # from the multipolygons 
        polyPointList = []
        for circleUnionExtInt in circleUnionExtIntList:
            outerLength = circleUnionExtInt.exterior.length
            #pointProg = np.arange(0,outerLength,np.sin(np.pi/2 - layerSpaceFraction*np.pi/6)*cellSize)
            pointProg = np.arange(0,outerLength,(0.8 - layerSpaceFraction*0.4)*cellSize) #To review the cell size
            for prog in pointProg:
                pointXY = list(circleUnionExtInt.exterior.interpolate(prog).xy)
                if self.overlapping:
                    polyPointList.append([pointXY[0][0],pointXY[1][0]])
                else:
                    pointXYPoint = Point(pointXY[0][0],pointXY[1][0])
                    if pointXYPoint.within(self.modelDis['activeArea'][-1]):
                        polyPointList.append([pointXY[0][0],pointXY[1][0]])
                
        circleUnionExtIntMpoly = MultiPolygon(circleUnionExtIntList)
        circleUnionExtWithIntMpoly = MultiPolygon(circleUnionExtWithIntList)

        return circleUnionExtWithIntMpoly, circleUnionExtIntMpoly, polyPointList

    def generateAllCircles(self):
        partialCircleUnionList = []
        partialCircleUnionInteriorList = []    

		#insert banner
        if isRunningInJupyter():
            printBannerHtml()
        else:
            printBannerText()


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

            #looping
            for index, cellSize in enumerate(cellSizeList):
                circleUnionInteriors, circleUnion, polyPointList = self.circlesAroundRefPoints(layer,index,cellSize)
                self.modelDis['vertexBuffer'] += polyPointList
                #for the last discretization
                if cellSize == np.array(cellSizeList).max():
                    #self.modelDis['circleUnion'] = circleUnion
                    partialCircleUnionList.append(circleUnion)
                    partialCircleUnionInteriorList.append(circleUnionInteriors)
                    #working with the final available geometry
                    lastGeometry = self.modelDis['activeArea'][-1]
                    partialActiveArea = lastGeometry.difference(circleUnionInteriors)
                    self.modelDis['activeArea'].append(partialActiveArea)

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

   