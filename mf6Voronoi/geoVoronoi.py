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
# /--------------- try to import dask_geopandas
try:
  import dask_geopandas as dgpd
  HAS_DASK = True
except ImportError:
  dgpd = None
  HAS_DASK = False
# ----------------/
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon, mapping
from collections import OrderedDict
from .utils import (intersectLimitLayer, 
                    isMultiGeometry,
                    isRunningInJupyter, 
                    printBannerHtml, 
                    printBannerText,
                    xy_in_poly, # \/ KMB additions \/
                    processVertexWithFilterCloseLimitDf,
                    exportMeshBuildFeaturesToShp,
                    getPolygonAndInteriors)

class createVoronoi():
    def __init__(self, meshName, maxRef, multiplier, overlapping=True, use_dask=False, nproc=1):
        #self.discGeoms = {}
        self.modelDis = {}
        self.modelDis['meshName'] = meshName
        self.modelDis['maxRef'] = maxRef
        self.modelDis['multiplier'] = multiplier
        self.overlapping = overlapping
        self.discLayers = {}
        self.settings = {
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

        # // -------------- KMB --------------
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
        # -------------- KMB -------------- //
            
    def generateOrgDistVertices(self, txtFile='', debug=False, out_dir="debug_org_dist"):
        vertexOrgPairList = []
        for layer, values in self.discLayers.items():
            vertexOrgPointList, vertexDistPointList, vertexDistPointList_asGeom = processVertexWithFilterCloseLimitDf(self.discLayers[layer]['layerRef'],
                                                                              self.discLayers[layer]['layerGeoms'],
                                                                              self.modelDis,
                                                                              use_dask=self.settings['use_dask'],
                                                                              nproc=self.settings['nproc'])
            vertexOrgPairList += vertexOrgPointList
            self.modelDis['vertexDist'][layer] = vertexDistPointList
            self.modelDis['vertexDistGeoms'][layer] = vertexDistPointList_asGeom
        self.modelDis['vertexOrg'] = vertexOrgPairList

        if txtFile != '':
            np.savetxt(txtFile+'_org',self.modelDis['vertexOrg'])
            np.savetxt(txtFile+'_dist',self.modelDis['vertexOrg'])

    def circlesAroundRefPoints(self,layer,last_indexBool,cellSize, debug=False):
        #first we create buffers around points and merge them
        #circleList = []
        #polyPointList = []
        crs = self.modelDis.get('crs', None)
        vertexDistGeoms = self.modelDis['vertexDistGeoms'][layer]
        layerSpaceList = self.discLayers[layer]['layerSpaceList']
        layerSpaceFraction = layerSpaceList.index(cellSize)/len(layerSpaceList)
        firstCellSize = layerSpaceList[0]

        # Vectorización del buffer y la unión espacial
        vertexDistGdf = gpd.GeoDataFrame(geometry=vertexDistGeoms, crs=crs)

        if self.settings['use_dask']:
            vertexDistPts = dgpd.from_geopandas(vertexDistGdf, npartitions=self.settings['nproc'])
            vertexDistPtsBuffer = vertexDistPts.buffer(cellSize)
            circleUnions = vertexDistPtsBuffer.union_all().compute()
        else:
            circleUnions = vertexDistGdf.buffer(cellSize).union_all()

        # for geom in self.modelDis['vertexDistGeoms'][layer]:
        #     #fixing for the first cell avoiding long cells
        #     #circle = geom.buffer(cellSize - firstCellSize/2) #Check this
        #     circle = geom.buffer(cellSize) #Check this
        #     circleList.append(circle)
        # circleUnions = unary_union(circleList)

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
        # polyPointList = []
        # for circleUnionExtInt in circleUnionExtIntList:
        #     outerLength = circleUnionExtInt.exterior.length
        #     #pointProg = np.arange(0,outerLength,np.sin(np.pi/2 - layerSpaceFraction*np.pi/6)*cellSize)
        #     pointProg = np.arange(0,outerLength,(0.8 - layerSpaceFraction*0.4)*cellSize) #To review the cell size
        #     for prog in pointProg:
        #         pointXY = list(circleUnionExtInt.exterior.interpolate(prog).xy)
        #         if self.overlapping:
        #             polyPointList.append([pointXY[0][0],pointXY[1][0]])
        #         else:
        #             pointXYPoint = Point(pointXY[0][0],pointXY[1][0])
        #             if pointXYPoint.within(self.modelDis['activeArea'][-1]):
        #                 polyPointList.append([pointXY[0][0],pointXY[1][0]])

        step_dist = (0.8 - layerSpaceFraction * 0.4) * cellSize
        
        raw_points = []
        for circleUnionExtInt in circleUnionExtIntList:
            exterior = circleUnionExtInt.exterior
            outerLength = exterior.length
            
            # Recrea exactamente los mismos pasos de distancia que la versión original
            pointProg = np.arange(0, outerLength, step_dist)
            
            # Interpolación vectorial (mucho más rápida que llamar .interpolate() uno por uno)
            pts = exterior.interpolate(pointProg)
            raw_points.extend(pts)

        # Convertimos los puntos a GeoDataFrame vectorizado
        sampled_pts_gdf = gpd.GeoDataFrame(geometry=raw_points, crs=crs)
        
        # Extraemos coordenadas como arreglo [[x1, y1], [x2, y2], ...]
        coords_array = np.column_stack([sampled_pts_gdf.geometry.x, sampled_pts_gdf.geometry.y])

        # Filtrado espacial vectorizado (idéntico al original)
        if self.overlapping:
            polyPointList = coords_array.tolist()
        else:
            active_area_geom = self.modelDis['activeArea'][-1]
            if self.settings['use_dask']:
                d_sampled = dgpd.from_geopandas(sampled_pts_gdf, npartitions=self.settings['nproc'])
                inside_mask = d_sampled.within(active_area_geom).compute().values
            else:
                inside_mask = sampled_pts_gdf.within(active_area_geom).values
            
            polyPointList = coords_array[inside_mask].tolist()
                
        circleUnionExtIntMpoly = MultiPolygon(circleUnionExtIntList)
        circleUnionExtWithIntMpoly = MultiPolygon(circleUnionExtWithIntList)
        
        return circleUnionExtWithIntMpoly, circleUnionExtIntMpoly, polyPointList


    def generateAllCircles(self, debug=False, verbose=True):
        partialCircleUnionList = []
        partialCircleUnionInteriorList = []
        if debug:
            self.modelDis['circleUnionByStep'] = {}
            self.modelDis['bufferPerCellSize'] = {}
            

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
                circleUnionInteriors, circleUnion, polyPointList = self.circlesAroundRefPoints(layer,last_index_bool,cellSize, debug)
                self.modelDis['vertexBuffer'] += polyPointList

                if debug:
                    # Store the circleUnion for the current step size
                    if cellSize not in self.modelDis['circleUnionByStep']:
                        self.modelDis['circleUnionByStep'][cellSize] = []
                    self.modelDis['circleUnionByStep'][cellSize].append(circleUnion)

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

    def getPointsMinMaxRef(self, verbose=True): # KMB, could also set verbose   at Class level and add to self.settings
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
    # -------------- KMB --------------//    
 
    
    def createPointCloud(self, verbose=True, debug=False, out_dir="debug_point_cloud"): 
        start = time.time()
        #Generate all circles and points on circle paths
        self.generateAllCircles(debug)
        #Distribute points over the max and min refinement areas
        self.getPointsMinMaxRef()
        #Compile all points
        # //------------------------ KMB
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
        # -------------- KMB -------------- //

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

        # --- EXPORTAR A SHAPEFILES SI DEBUG=TRUE ---
            if debug:
                print(
                    f"\n[DEBUG] Exportando todos los sets del point cloud:"
                    f" {out_dir}"
                )
                exportMeshBuildFeaturesToShp(self.modelDis, out_dir=out_dir)


    def generateVoronoi(self, shapePath=None):
        print('\n/----Generation of the voronoi mesh----/')
        start = time.time()
        #create a multipoint object
        pointMulti = MultiPoint(self.modelDis['vertexTotal'])
        #original regions
        regions = voronoi_diagram(pointMulti)

        # // -----------------KMB
        layer_df = gpd.GeoDataFrame(np.arange(len(regions.geoms)),columns=['id'],geometry=list(regions.geoms),crs=self.modelDis['crs'])
        
        # Update geometries intersecting limit polygon
        if self.settings['use_dask']:
            intersecting_bool = dgpd.from_geopandas(layer_df,self.settings['nproc']).intersects(self.modelDis['limitGeometry'].exterior).compute()
            layer_df.loc[intersecting_bool,'geometry'] = dgpd.from_geopandas(layer_df.loc[intersecting_bool],self.settings['nproc']).intersection(self.modelDis['limitGeometry']).compute()
        else:
            intersecting_bool = layer_df.intersects(self.modelDis['limitGeometry'].exterior)
            layer_df.loc[intersecting_bool,'geometry'] = layer_df.loc[intersecting_bool].intersection(self.modelDis['limitGeometry'])

        #clippedRegions = layer_df.explode(index_parts=False).geometry.values.tolist()
        # --- FIX: Explode GeoSeries directly and filter empty geometries ---
        exploded_geoms = layer_df.geometry.explode(ignore_index=True)
        clippedRegions = exploded_geoms[~exploded_geoms.is_empty].tolist()
        
        clippedRegionsMulti = MultiPolygon(clippedRegions)
        self.modelDis['voronoiRegions'] = clippedRegionsMulti
        end = time.time()
        print('\nTime required for voronoi generation: %.2f seconds \n'%(end - start), flush=True)

        if shapePath is not None: # option to save shapefile directly
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