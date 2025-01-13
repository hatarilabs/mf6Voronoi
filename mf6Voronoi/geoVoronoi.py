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
    def __init__(self, maxRef, multiplier):
        #self.discGeoms = {}
        self.modelDis = {}
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
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
        self.modelDis['vertexBuffer'] = []
        self.modelDis['crs'] = limitShape.crs

    # we shift from minRef to layerRef
    # def defineParameters(self, maxRef=500, minRef=50, multiplier=1):
    def defineParameters(self, maxRef=500, multiplier=1):
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier

    #here we add the layerRef to the function
    #def addLayer(self, name, shapePath):
    def addLayer(self, layerName, shapePath, layerRef):
        #Add layers for mesh definition
        #This feature also clips and store the geometry
        #geomList is allways a Python of Shapely geometries
      
        spatialDf = gpd.read_file(shapePath)   

        #get the ref and geoms as a list
        self.discLayers[layerName] = {'layerRef':layerRef,
                                      'layerGeoms':[]}  
        
        #auxiliary funtion to intersect:
        def intersectLimitLayer(discLayerGeom):
            discGeomList = []  
            #generic 
            if self.isMultiGeometry(discLayerGeom):
                for partGeom in discLayerGeom.geoms:
                    discGeomClip =  self.modelDis['limitGeometry'].intersection(partGeom)
                    if not discGeomClip.is_empty:
                        discGeomList.append(discGeomClip)
            else:
                discGeomClip =  self.modelDis['limitGeometry'].intersection(discLayerGeom)
                if not discGeomClip.is_empty:
                    discGeomList.append(discGeomClip)

            unaryGeom = unary_union(discGeomList)
            # print('aaa')
            # print(type(unaryGeom))
            # print(isinstance(unaryGeom, (MultiLineString)))
            # print(self.isMultiGeometry(unaryGeom))
            if self.isMultiGeometry(unaryGeom):
                # print('bbb')
                unaryFilter = [geom for geom in unaryGeom.geoms]
            else:
                unaryFilter = [unaryGeom]

            return unaryFilter

        #looping over the shapefile
        for spatialIndex, spatialRow in spatialDf.iterrows():
            if spatialRow.geometry.is_valid:
                geomGeom = spatialRow.geometry
                unaryFilter = intersectLimitLayer(geomGeom)
                self.discLayers[layerName]['layerGeoms'] += unaryFilter
            else:
                print('You are working with a uncompatible geometry. Remember to use single parts')
                print('Check this file: %s \n'%shapePath)
                sys.exit()
            
    def orgVertexAsList(self, layerGeoms):
        #get only the original vertices inside the model limit
        vertexList = []

        for layerGeom in layerGeoms:
            if layerGeom.geom_type == 'Polygon':
                pointObject = layerGeom.exterior.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                for index, point in enumerate(pointList):
                    vertexList.append(point)
            elif layerGeom.geom_type == 'LineString':
                pointObject = layerGeom.coords.xy
                pointList = list(zip(pointObject[0],pointObject[1]))
                for index, point in enumerate(pointList):
                    vertexList.append(point)
            elif layerGeom.geom_type == 'Point':
                pointObject = layerGeom.coords.xy
                vertexList.append((pointObject[0][0],pointObject[1][0]))
            else:
                print(layerGeom)
                print('/-----Problem has been bound when extracting org vertex-----/')

        return vertexList

    def distributedVertexAsList(self, layerGeoms, layerRef):
        #distribute vertices along the layer paths
        vertexList = []

        for layerGeom in layerGeoms:
            if layerGeom.geom_type == 'Polygon':
                polyLength = layerGeom.exterior.length
                pointProg = np.arange(0,polyLength,layerRef)
                for prog in pointProg:
                    pointXY = list(layerGeom.exterior.interpolate(prog).xy)
                    vertexList.append([pointXY[0][0],pointXY[1][0]])
            elif layerGeom.geom_type == 'LineString':
                lineLength = layerGeom.length
                pointProg = np.arange(0,lineLength,layerRef)
                for prog in pointProg:
                    pointXY = list(layerGeom.interpolate(prog).xy)
                    vertexList.append([pointXY[0][0],pointXY[1][0]])
            elif layerGeom.geom_type == 'Point':
                pointObject = layerGeom.coords.xy
                vertexList.append((pointObject[0][0],pointObject[1][0]))
            else:
                print('/-----Problem has been bound when extracting dist vertex-----/')

        return vertexList

    def generateOrgDistVertices(self, txtFile=''):
        start = time.time()
        vertexOrgPairList = []
        #vertexDistPairList = []
        for layer, values in self.discLayers.items():
            vertexOrgPairList += self.orgVertexAsList(values['layerGeoms'])
            layerGeoms = values['layerGeoms']
            layerRef = values['layerRef']
            self.modelDis['vertexDist'][layer] = self.distributedVertexAsList(layerGeoms, layerRef)
        self.modelDis['vertexOrg'] = vertexOrgPairList

        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,layer,indexRef,cellSize):
        #first we create buffers around points and merge them
        circleList = []
        polyPointList = []
        for point in self.modelDis['vertexDist'][layer]:
            circle = Point(point).buffer(cellSize)
            circleList.append(circle)
        circleUnions = unary_union(circleList)
        
        def getPolygonAndInteriors(polyGeom):
            exteriorInteriorPolys = [polyGeom] + [Polygon(ring) for ring in polyGeom.interiors]
            return exteriorInteriorPolys
         
        circleUnionExtIntList = []
        if circleUnions.geom_type == 'MultiPolygon':
            for circleUnion in circleUnions.geoms:
                circleUnionExtIntList += getPolygonAndInteriors(circleUnion)
        elif circleUnions.geom_type == 'Polygon':
            circleUnionExtIntList += getPolygonAndInteriors(circleUnions)

        # from the multipolygons 
        polyPointList = []
        for circleUnionExtInt in circleUnionExtIntList:
            outerLength = circleUnionExtInt.exterior.length
            if indexRef%2 == 0:
                pointProg = np.arange(0,outerLength,np.pi*cellSize/3)
            else:
                pointProg = np.arange(np.pi*cellSize/6,outerLength+np.pi*cellSize/6,np.pi*cellSize/3)
            for prog in pointProg:
                pointXY = list(circleUnionExtInt.exterior.interpolate(prog).xy)
                polyPointList.append([pointXY[0][0],pointXY[1][0]])

        circleUnionExtIntMpoly = MultiPolygon(circleUnionExtIntList)
        return circleUnionExtIntMpoly, polyPointList

    def generateAllCircles(self):
        partialCircleUnionList = []

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

            self.discLayers[layer]['layerSpaceListt'] = cellSizeList

            print('\n/--------Layer %s discretization-------/'%layer)
            print('Progressive cell size list: %s m.'%str(cellSizeList))

            for index, cellSize in enumerate(cellSizeList):
                circleUnion, polyPointList = self.circlesAroundRefPoints(layer,index,cellSize)
                refBuffer = gpd.GeoSeries(circleUnion)
                self.modelDis['vertexBuffer'] += polyPointList
                #here we use the maximum progressive refinement
                #if ref == self.modelDis['refSizeList'].max():
                if cellSize == np.array(cellSizeList).max():
                    #self.modelDis['circleUnion'] = circleUnion
                    partialCircleUnionList.append(circleUnion)

        totalCircleUnion = unary_union(partialCircleUnionList)
        self.modelDis['circleUnion'] = totalCircleUnion

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
        for key, value in self.discLayers.items():
            for layerGeom in value['layerGeoms']:
                if layerGeom.geom_type == 'Polygon':
                    transPoly = outerPoly.difference(layerGeom)
                    if limitPoly.area == transPoly.area:
                        outerPoly.geom.interior += layerGeom
                    elif limitPoly.area > transPoly.area:
                        outerPoly = outerPoly.difference(layerGeom)

            # if value['type'] == 'Polygon':
            #     for poly in value['geomList']:
            #         transPoly = outerPoly.difference(poly)
            #         if limitPoly.area == transPoly.area:
            #             outerPoly.geom.interior += poly
            #         elif limitPoly.area > transPoly.area:
            #             outerPoly = outerPoly.difference(poly)
                                 
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
                
            # polyDict = self.discGeoms[polyKey]
            # if polyDict['type']=='Polygon':
            #     for poly in polyDict['geomList']:
            #         bounds = poly.exterior.bounds
            #         minRefXList = np.arange(bounds[0]+polyDict['layerRef'],bounds[2],polyDict['layerRef'])
            #         minRefYList = np.arange(bounds[1]+polyDict['layerRef'],bounds[3],polyDict['layerRef'])

            #         for xCoord in minRefXList:
            #             for yCoord in minRefYList:
            #                 refPoint = Point(xCoord,yCoord)
            #                 if poly.contains(refPoint):
            #                     self.modelDis['vertexMinRef'].append((xCoord,yCoord))

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
