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
import dask_geopandas as dgpd # KMB
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from collections import OrderedDict
from .utils import (processVertexFilterCloseLimit, 
                    intersectLimitLayer, 
                    isMultiGeometry,
                    isRunningInJupyter, 
                    printBannerHtml, 
                    printBannerText,
                    xy_in_poly, # \/ KMB additions \/
                    processVertexFilterCloseLimit_df,
                    getPolygonAndInteriors)

class createVoronoi():
    def __init__(self, meshName, maxRef, multiplier, overlapping=True, use_gpd=False, use_dask=False, nproc=1): # KMB
        #self.discGeoms = {}
        self.modelDis = {}
        self.modelDis['meshName'] = meshName
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
        self.overlapping = overlapping
        self.discLayers = {}
        self.settings = {'use_gpd':use_gpd,
                         'use_dask': use_dask,
                         'nproc': nproc} # number of processors to use

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
        
        # -------------- KMB --------------
        if self.settings['use_gpd']: #KMB
        
            # find intersection with limit layer
            if self.settings['use_dask']:
                unaryFilter = dgpd.from_geopandas(spatialDf,self.settings['nproc']).intersection(self.modelDis['limitGeometry']).compute()
            else:
                unaryFilter = spatialDf.intersection(self.modelDis['limitGeometry'])
            
            # Remove empty geometries after intersection
            unaryFilter = unaryFilter[~unaryFilter.is_empty]
            
            # keep intersection geometries as a list and explode multipolygons
            self.discLayers[layerName]['layerGeoms'] = unaryFilter.explode().tolist()
            
            # store geometry type
            self.discLayers[layerName]['layerType'] = spatialDf.iloc[0].geometry.geom_type
            
            if len(self.discLayers[layerName]['layerType']) == 0:
                print('You are working with a incompatible geometry. Remember to use single parts')
                print('Check this file: %s \n'%shapePath)
                sys.exit()
        # \-------------- KMB --------------
        else:
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
                    print('You are working with a incompatible geometry. Remember to use single parts')
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
            # -------------- KMB --------------
            if self.settings['use_gpd']: # KMB
                vOPL,vertexList,vertexGeomList = processVertexFilterCloseLimit_df(self.discLayers[layer]['layerRef'],
                                                                              self.discLayers[layer]['layerGeoms'],
                                                                              self.modelDis,
                                                                              use_dask=self.settings['use_dask'],nproc=self.settings['nproc'])
                vertexOrgPairList += vOPL
            else:
                vertexOrgPairList += self.orgVertexAsList(layer)
                vertexList,vertexGeomList = self.distributedVertexAsList(layer) # KMB, run once and extract outputs to dictionary
                #self.modelDis['vertexDist'][layer] = self.distributedVertexAsList(layer)[0] # KMB, this runs distributedVertexAsList first time
                #self.modelDis['vertexDistGeoms'][layer] = self.distributedVertexAsList(layer)[1] # KMB, this runs distributedVertexAsList second time
            
            self.modelDis['vertexDist'][layer] = vertexList # KMB
            self.modelDis['vertexDistGeoms'][layer] = vertexGeomList # KMB
            # /-------------- KMB -------------- 
        self.modelDis['vertexOrg'] = vertexOrgPairList

        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,layer,last_indexBool,cellSize):
        
        #first we create buffers around points and merge them
        circleList = []
        polyPointList = []
        layerSpaceList = self.discLayers[layer]['layerSpaceList']
        layerSpaceFraction = layerSpaceList.index(cellSize)/len(layerSpaceList)
        firstCellSize = layerSpaceList[0]
        
        # -------------- KMB --------------
        if self.settings['use_gpd']:
            layer_df = gpd.GeoDataFrame(geometry=self.modelDis['vertexDistGeoms'][layer])
        
            if self.settings['use_dask']:
                layer_df['geometry'] = dgpd.from_geopandas(layer_df,self.settings['nproc']).buffer(cellSize).compute()
                circleUnions = dgpd.from_geopandas(layer_df,self.settings['nproc']).union_all().compute()
            else:
                layer_df['geometry'] = layer_df.buffer(cellSize)
                circleUnions = layer_df.union_all()
            
            if last_indexBool:
                circle_df = gpd.GeoDataFrame([0],columns=['id'],geometry=[circleUnions]).explode()
                circleUnionExtWithIntList = circle_df.geometry.values.tolist()
            
            # # Add interior ring geometries
            # circle_types = np.unique(circle_df.geom_type)
            # if len(circle_types)==1 and circle_types[0]=='Polygon':
            #     combo_df = circle_df.copy()
            # else:
            #     interiors_df = circle_df.interiors
            #     interiors_df.dropna(inplace=True)
            #     if interiors_df.shape[0] > 0:
            #         interior_inds = [ind for ind in range(circle_df.shape[0]) if len(interiors_df.iloc[ind])>0]
            #         interior_geoms = []
            #         if len(interior_inds) > 0:
            #             for interior_ind in interior_inds:
            #                 # Convert to polygons
            #                 interior_geoms.extend([Polygon(igeom) for igeom in interiors_df.iloc[interior_ind]])
                    
            #             combo_df = gpd.pd.concat([circle_df,gpd.GeoDataFrame(np.arange(len(interior_geoms)),columns=['id'],geometry=interior_geoms)],ignore_index=True)
            #         else:
            #             combo_df = circle_df.copy()
            #     else:
            #         combo_df = circle_df.copy()
                
                
            # circleUnionExtIntList = combo_df.explode().geometry.values.tolist()
            
            circleUnionExtIntList = []
            if circleUnions.geom_type == 'MultiPolygon':
                for circleUnion in circleUnions.geoms:
                    circleUnionExtIntList += getPolygonAndInteriors(circleUnion)
            elif circleUnions.geom_type == 'Polygon':
                circleUnionExtIntList += getPolygonAndInteriors(circleUnions)
            
            combo_in_df = gpd.GeoDataFrame(np.arange(len(circleUnionExtIntList)),
                                    geometry=circleUnionExtIntList)
            
            combo_in_df['geometry'] = combo_in_df.segmentize((0.8 - layerSpaceFraction*0.4)*cellSize)
            if self.overlapping:
                polyPointList = combo_in_df.get_coordinates().values.tolist()
            else:
                temp_points = combo_in_df.get_coordinates().values
                all_points_df = gpd.GeoDataFrame(geometry=gpd.points_from_xy(temp_points[:,0],temp_points[:,1]))
                polyPointList = all_points_df.loc[all_points_df.within(self.modelDis['activeArea'][-1])].get_coordinates().values.tolist()

            
            # if self.settings['use_dask']:
            #     combo_in_df['geometry'] = dgpd.from_geopandas(combo_in_df,self.settings['nproc']).difference(self.modelDis['activeArea'][-1]).compute()
            # else:
            #     combo_in_df['geometry'] = combo_in_df.difference(self.modelDis['activeArea'][-1])

            
            # polyPointList = combo_in_df.get_coordinates().to_numpy().tolist()
        
        # -------------- KMB --------------    
        else:
            for geom in self.modelDis['vertexDistGeoms'][layer]:
                #fixing for the first cell avoiding long cells
                #circle = geom.buffer(cellSize - firstCellSize/2) #Check this
                circle = geom.buffer(cellSize) #Check this
                circleList.append(circle)
            circleUnions = unary_union(circleList)
            
            # KMB - import instead of defining in function call
            # def getPolygonAndInteriors(polyGeom):
            #     exteriorInteriorPolys = [polyGeom] + [Polygon(ring) for ring in polyGeom.interiors]
            #     return exteriorInteriorPolys
             
            circleUnionExtIntList = []
            circleUnionExtWithIntList = []
            if circleUnions.geom_type == 'MultiPolygon':
                for circleUnion in circleUnions.geoms:
                    circleUnionExtIntList += getPolygonAndInteriors(circleUnion)
                    if last_indexBool: # KMB
                        circleUnionExtWithIntList.append(circleUnion)
            elif circleUnions.geom_type == 'Polygon':
                circleUnionExtIntList += getPolygonAndInteriors(circleUnions)
                if last_indexBool: # KMB
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
        
        if last_indexBool: # KMB
            circleUnionExtIntMpoly = MultiPolygon(circleUnionExtIntList)
            circleUnionExtWithIntMpoly = MultiPolygon(circleUnionExtWithIntList)
        else:
            circleUnionExtWithIntMpoly, circleUnionExtIntMpoly = None,None # no need to spend time on these
        
        return circleUnionExtWithIntMpoly, circleUnionExtIntMpoly, polyPointList

    def generateAllCircles(self, use_banner=True, verbose=True): # KMB
        partialCircleUnionList = []
        partialCircleUnionInteriorList = []    

        if use_banner: # KMB :)
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
            
            if verbose:        # KMB     
                print('\n/--------Layer %s discretization-------/'%layer)
                print('Progressive cell size list: %s m.'%str(cellSizeList))

            #looping
            for index, cellSize in enumerate(cellSizeList):
                # KMB - only run circleUnionInterios and circleUnion calculations for last index
                last_index_bool = cellSize == np.array(cellSizeList).max()
                circleUnionInteriors, circleUnion, polyPointList = self.circlesAroundRefPoints(layer,last_index_bool,cellSize)
                self.modelDis['vertexBuffer'] += polyPointList
                #for the last discretization
                if last_index_bool: # KMB
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

    def getPointsMinMaxRef(self, verbose=True): # KMB, could also set verbose at Class level and add to self.settings

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
        
        # -------------- KMB --------------, similiarities in parts, but easier to separate entirely for now
        if self.settings['use_gpd']:
            innerPolys_df = gpd.GeoDataFrame([0],columns=['id'],geometry=[innerPolys]).explode()
            if verbose:
                print("Find maxiumum reference points") 
            
            innerPolys_ext_bool = innerPolys_df.intersects(outerPoly.exterior)
            all_diffs_df = innerPolys_df[innerPolys_ext_bool].copy()
            all_interiors = innerPolys_df[~innerPolys_ext_bool].geometry.values.tolist()
            
            if verbose:
                print("Find maxiumum reference points, disc polys") # slow
            
            #working with mesh disc polys
            for key, value in self.discLayers.items():
                
                layer_df = gpd.GeoDataFrame(np.arange(len(value['layerGeoms'])),columns=['id'],geometry=value['layerGeoms'])
                
                if verbose:
                    print("Collect interior geometries, disc polys")
                
                # Identify interior geometries and add as interiors to outerPoly
                if self.settings['use_dask']:
                    within_bool = dgpd.from_geopandas(layer_df,self.settings['nproc']).within(limitPoly).compute()
                    internal_df = layer_df.loc[within_bool]
                    
                    # Unify overlapping polygons so only largest internal geometries exist
                    internal_geoms = dgpd.from_geopandas(internal_df,self.settings['nproc']).union_all().compute()
                    
                else:
                    within_bool = layer_df.within(limitPoly)
                    internal_geoms = layer_df.loc[within_bool].union_all()
            
            # Need to add interiors
            all_interiors.extend([internal_geoms])
            outerPoly = outerPoly.difference(unary_union(all_interiors)).buffer(0) # buffer to fix invalid issues
            
            if verbose:
                print("Update outerPoly by removing intersecting layer geometries, disc polys")
            
            outer_df = gpd.GeoDataFrame([0],columns=['id'],geometry=[outerPoly]).explode()
            
            # Update outerPoly to have edge geometries cut out of it
            edge_df = layer_df.loc[np.invert(within_bool.values)] # if not within, then must intersect
            
            # Add inner rings that intersect outerPoly to reduce its size
            edge_df = gpd.pd.concat([edge_df,all_diffs_df],ignore_index=True)
            
            if self.settings['use_dask']:
                edge_geom = dgpd.from_geopandas(edge_df,self.settings['nproc']).union_all().compute().buffer(0)
                temp_df = outer_df.difference(edge_geom)
                outerPoly = dgpd.from_geopandas(temp_df,self.settings['nproc']).union_all().compute()
            else:
                edge_geom = edge_df.union_all().buffer(0)
                temp_df = outer_df.difference(edge_geom)
                outerPoly = temp_df.union_all()
                
            self.modelDis['pointsMaxRefPoly']=outerPoly

            #creating points of coarse grid
            maxRefXList = np.arange(self.modelDis['xMin']+minRef,self.modelDis['xMax'],maxRef)
            maxRefYList = np.arange(self.modelDis['yMin']+minRef,self.modelDis['yMax'],maxRef)
            
            maxX,maxY = np.meshgrid(maxRefXList,maxRefYList)
            
            self.modelDis['vertexMaxRef'] = xy_in_poly(np.column_stack([maxX.ravel(),maxY.ravel()]).tolist(),outerPoly,
                                             use_dask=self.settings['use_dask'], nproc=self.settings['nproc'])
            #for min ref points
            if verbose:
                print('Find minimum reference points') # slow
                
            for key, value in self.discLayers.items():
                
                layer_df = gpd.GeoDataFrame(np.arange(len(value['layerGeoms'])),columns=['id'],geometry=value['layerGeoms'])
                all_bound = layer_df.bounds.values
                layerRef = value['layerRef']
                
                multipt_list = [MultiPoint(np.column_stack(list(map(np.ravel,np.meshgrid(np.arange(bounds[0]+layerRef,bounds[2],layerRef),
                                                   np.arange(bounds[1]+layerRef,bounds[3],layerRef))))).tolist()) for bounds in all_bound]
            
                more_points_needed = [not igeom.is_empty for igeom in multipt_list]
                multipt_list = [mp for imp,mp in enumerate(multipt_list) if more_points_needed[imp]] # select only needed multipoints
                
                layer_df = layer_df.iloc[more_points_needed]
                
                multipt_df = gpd.GeoDataFrame(geometry=multipt_list)
                self.modelDis['vertexMinRef'].extend(multipt_df.intersection(layer_df,align=False).get_coordinates().values.tolist())
            
            
        # -------------- KMB --------------    
        else:
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

    def createPointCloud(self, verbose=True, use_banner=True):
        start = time.time()
        #Generate all circles and points on circle paths
        self.generateAllCircles(verbose=verbose, use_banner=use_banner)
        #Distribute points over the max and min refinement areas
        self.getPointsMinMaxRef(verbose=verbose)
        
        # -------------- KMB -------------- 
        if self.settings['use_gpd']:
            totalRawPoints = [np.array(self.modelDis['vertexDist'][key]) for key in self.modelDis['vertexDist'] if len(self.modelDis['vertexDist'][key])>0]

            if len(self.modelDis['vertexBuffer'])>0:
                totalRawPoints.append(np.array(self.modelDis['vertexBuffer']))
            
            if len(self.modelDis['vertexMaxRef'])>0:
                totalRawPoints.append(np.array(self.modelDis['vertexMaxRef']))
            
            if len(self.modelDis['vertexMinRef'])>0:
                totalRawPoints.append(np.array(self.modelDis['vertexMinRef']))

            #check if points are inside limit polygon
            points_array = np.vstack(totalRawPoints)
            totalDefPoints = xy_in_poly(points_array.tolist(),self.modelDis['limitGeometry'],
                                        use_dask=self.settings['use_dask'], nproc=self.settings['nproc'])
            # -------------- KMB -------------- 
        else:
        
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

        if verbose: #KMB
            print('\n/----Sumary of points for voronoi meshing----/')
            print('Distributed points from layers: %d'%len(self.modelDis['vertexDist']))
            print('Points from layer buffers: %d'%len(self.modelDis['vertexBuffer']))
            print('Points from max refinement areas: %d'%len(self.modelDis['vertexMaxRef']))
            print('Points from min refinement areas: %d'%len(self.modelDis['vertexMinRef']))
            print('Total points inside the limit: %d'%len(self.modelDis['vertexTotal']))
            print('/--------------------------------------------/')
            end = time.time()
            print('\nTime required for point generation: %.2f seconds \n'%(end - start), flush=True)

    def generateVoronoi(self, shapePath=None): # KMB, give option to save file when gpd is already loaded
        print('\n/----Generation of the voronoi mesh----/')
        start = time.time()
        #create a multipoint object
        pointMulti = MultiPoint(self.modelDis['vertexTotal'])
        #original regions
        regions = voronoi_diagram(pointMulti)
        
        # -------------- KMB -------------- 
        if self.settings['use_gpd']:
            layer_df = gpd.GeoDataFrame(np.arange(len(regions.geoms)),columns=['id'],geometry=list(regions.geoms),crs=self.modelDis['crs'])
            
            
            # Update geometries intersecting limit polygon
            if self.settings['use_dask']:
                intersecting_bool = dgpd.from_geopandas(layer_df,nproc).intersects(self.modelDis['limitGeometry'].exterior).compute()
                layer_df.loc[intersecting_bool,'geometry'] = dgpd.from_geopandas(layer_df.loc[intersecting_bool],self.settings['nproc']).intersection(self.modelDis['limitGeometry']).compute()
            else:
                intersecting_bool = layer_df.intersects(self.modelDis['limitGeometry'].exterior)
                layer_df.loc[intersecting_bool,'geometry'] = layer_df.loc[intersecting_bool].intersection(self.modelDis['limitGeometry'])
            
            clippedRegions = layer_df.explode().geometry.values.tolist()
        
        # -------------- KMB -------------- 
        else:
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
        
        # -------------- KMB -------------- 
        
        if shapePath is not None: # option to save shapefile directly
            if not self.settings['use_gpd']:
                layer_df = gpd.GeoDataFrame(np.arange(len(clippedRegions)),columns=['id'],geometry=clippedRegions,crs=self.modelDis['crs'])
            
            layer_df.to_file(shapePath)

    def checkVoronoiQuality(self, threshold = 0.001):
        print('\n/----Performing quality verification of voronoi mesh----/')
        self.modelDis['fixPoints'] = []
        # empty list to store distances
        shortSides = []
        
        for index, poly in enumerate(self.modelDis['voronoiRegions'].geoms):
            polyCoordList = []
            x,y = poly.exterior.coords.xy
            polyCoordList.append(list(zip(x,y)))
            if poly.interiors[:] != []:
                for interior in poly.interiors:
                    polyCoordList.append(interior.coords[:])

            # loopo over polygon on polygon list
            for polyCoord in polyCoordList:
                #looping over sides
                for i in range(len(polyCoord) - 1):
                    p1 = polyCoord[i]
                    p2 = polyCoord[i + 1]
                    edge = LineString([p1, p2])
                    length = edge.length
                    if length < threshold:
                        xMean = (p1[0] + p2[0])/2
                        yMean = (p1[1] + p2[1])/2
                        self.modelDis['fixPoints'].append([xMean,yMean])
                        shortSides.append((p1, p2, length))
            
        # Output short sides

        if len(shortSides) == 0:
            print("Your mesh has no edges shorter than your threshold")
        else:
            for side in shortSides:
                print(f"Short side on polygon: {index} with length = {side[2]:.5f}")

    def fixVoronoiShortSides(self):
        self.modelDis['vertexTotal'] = self.modelDis['vertexTotal'] + self.modelDis['fixPoints']