import json
import numpy as np
import fiona
# KMB additions
import gzip
import geopandas as gpd
# --------
from shapely.geometry import Polygon
from tqdm import tqdm
from .utils import isRunningInJupyter, printBannerHtml, printBannerText, findIndex, unique_rows # KMB

class meshShape():
    def __init__(self,path_to_mesh_shp=None,mesh_df=None): # KMB addition
        self.mesh=path_to_mesh_shp
        self.mesh_df = mesh_df # KMB
        self.disvDict = {}
        self.spatialIndexDict = {}
        
    def get_gridprops_disv(self,use_gpd=False, intervalNumber=10,save_spatial_index=True, verbose=False): # KMB additions
    
        if use_gpd: # KMB
            if self.mesh_df is None:
                if self.mesh is None:
                    print('Either a mesh_df or mesh shapefile name is needed. Exiting.')
                    exit()
                    
                self.mesh_df = gpd.read_file(self.mesh)
            
            self.mesh_df = self.mesh_df.explode()
            
            if verbose:
                print('\nCreating a unique list of vertices [[x1,y1],[x2,y2],...]')
            
            allVerticesList = self.mesh_df.get_coordinates().values
        
            if verbose:
                print('\nExtracting cell2d data')
                
            polyLenVerticies = self.mesh_df.count_coordinates().values.tolist()  # number of verticies in each polygon
            cumul_polyLen = np.hstack([0,np.cumsum(polyLenVerticies)])
            polyIndicies = np.column_stack([cumul_polyLen[:-1],cumul_polyLen[1:]])
            polyAllIndicies = [np.arange(istart,iend) for istart,iend in polyIndicies]

            uniqueVerticesArray,uind,uinv = unique_rows(allVerticesList,return_inverse=True)
            uniqueVerticesList = uniqueVerticesArray.tolist()
            centroids = self.mesh_df.geometry.centroid.apply(lambda c: c.coords[0]).tolist() # (x,y) coordinates of mesh centroids
            
            cell2dArrays = [[ind,*c,nvert,*(uinv[poly_inds].ravel()).tolist()] for ind,(c,nvert,poly_inds) in enumerate(zip(centroids,polyLenVerticies,polyAllIndicies))] # index,centroid,len(vertex_indexes), vertex_indexes
            
            self.disvDict['ncpl'] = len(self.mesh_df)
            self.disvDict['nvert'] = len(uniqueVerticesList)
            self.disvDict['uniqueVerticesList'] = uniqueVerticesList
            self.disvDict['vertices'] = np.column_stack([np.arange(len(uniqueVerticesList)),uniqueVerticesList]).tolist()
            self.disvDict['cell2d'] = cell2dArrays
            self.disvDict['centroids'] = centroids
            
            if save_spatial_index:
                if verbose:
                    print('\nExtracting grid spatial index data')
                
                # Get grid index
                meshBounds = self.mesh_df.total_bounds
                gridXarray = np.linspace(meshBounds[0],meshBounds[2],intervalNumber+1)
                gridYarray = np.linspace(meshBounds[1],meshBounds[3],intervalNumber+1)
                
                poly_bounds = self.mesh_df.bounds # bounds of each cell
                xInterBeg = np.searchsorted(gridXarray,poly_bounds['minx'].values)
                xInterEnd = np.searchsorted(gridXarray,poly_bounds['maxx'].values)
                yInterBeg = np.searchsorted(gridYarray,poly_bounds['miny'].values)
                yInterEnd = np.searchsorted(gridYarray,poly_bounds['maxy'].values)
                gridIndexList = np.column_stack([xInterBeg,xInterEnd,yInterBeg,yInterEnd]).reshape([len(self.mesh_df),2,2]).tolist()
            
        else: # original
            vorMesh = fiona.open(self.mesh)
            # Get grid index
            # intervalNumber = 10 # KMB
            meshBounds = vorMesh.bounds
            gridXarray = np.linspace(meshBounds[0],meshBounds[2],intervalNumber+1)
            gridYarray = np.linspace(meshBounds[1],meshBounds[3],intervalNumber+1)

            totalVerticesList = []
            cell2dArrays = []
            polygonCentroidList = []
            gridIndexList =[]
            #defining function
            # KMB - load from utils
            #def findIndex(var, coordArray):
            #    for interval in range(intervalNumber):
            #        if var > coordArray[interval] and var < coordArray[interval+1]:
            #        if var >= coordArray[interval] and var < coordArray[interval+1]:
            #            return interval
            #            break

            # #insert banner
            # if isRunningInJupyter():
            #   printBannerHtml()
            # else:
            #   printBannerText()

            print('\nCreating a unique list of vertices [[x1,y1],[x2,y2],...]')
            for index, row in enumerate(tqdm(vorMesh, total= len(vorMesh))):
                #vertices xy
                if len(row['geometry']['coordinates']) == 1:
                    #print(row['geometry']['coordinates'])
                    xyList = [[i[0],i[1]] for i in row['geometry']['coordinates'][0]]
                    totalVerticesList += xyList
                elif len(row['geometry']['coordinates']) > 1:
                    print(row['geometry']['coordinates'])
                    print(index)
                    for vertexList in row['geometry']['coordinates']:
                        #print(vertexList)
                        xyList = [[i[0],i[1]] for i in vertexList[0]]
                        #print(xyList)
                        totalVerticesList += xyList
                else:
                    pass
            uniqueVerticesArray = np.unique(np.array(totalVerticesList), axis=0)
            uniqueVerticesList = uniqueVerticesArray.tolist()

            vertexIndexDict = {}
            for index, vertex in enumerate(uniqueVerticesList):
                strVertex = str(vertex)
                vertexIndexDict[strVertex]=index
                
                
            print('\nExtracting cell2d data and grid index')
            centroids=[]
            for index,row in enumerate(tqdm(vorMesh)):#.iterrows(), total= vorMesh.shape[0]):
                rowCoords = row['geometry']['coordinates'][0]
                rowPoly = Polygon(rowCoords)
                #print(index)
                #print(rowPoly.bounds)
                #coords = rowGeometry.exterior.coords
                #cell2d array
                cellArray = []
                #add index
                cellArray.append(index)
                #add centroid
                cellArray += list(rowPoly.centroid.coords[0])

                centroids.append(tuple(rowPoly.centroid.coords[0]))
                #working with vertices number and vertex
                vertexIndexList = []
                for vertex in rowCoords:
                    #print(vertex)
                    strVertex = str(list(vertex))
                    #print(vertexIndexDict[strVertex])
                    vertexIndexList.append(vertexIndexDict[strVertex])
                cellArray.append(len(vertexIndexList))
                cellArray += vertexIndexList
                cell2dArrays.append(cellArray)
                #get grid index
                xmin = rowPoly.bounds[0] #min(coords.xy[0])
                xmax = rowPoly.bounds[2] #max(coords.xy[0])
                ymin = rowPoly.bounds[1] #min(coords.xy[1])
                ymax = rowPoly.bounds[3] #max(coords.xy[1])
                xInterBeg = findIndex(xmin,gridXarray)
                xInterEnd = findIndex(xmax,gridXarray)
                yInterBeg = findIndex(ymin,gridYarray)
                yInterEnd = findIndex(ymax,gridYarray)
                gridIndexList.append([[xInterBeg,xInterEnd],[yInterBeg,yInterEnd]])

            # uniqueVerticesArray = np.unique(np.array(totalVerticesList),axis=0) # KMB, repeated from above
            # uniqueVerticesList = uniqueVerticesArray.tolist() # KMB, repeated from above, can remove
            indexedVerticesList = [[index, row[0], row[1]] for index, row in enumerate(uniqueVerticesList)]

            
            self.disvDict['ncpl'] = len(vorMesh)
            self.disvDict['nvert'] = len(uniqueVerticesList)
            self.disvDict['uniqueVerticesList']=uniqueVerticesList
            self.disvDict['vertices']=indexedVerticesList
            self.disvDict['cell2d'] = cell2dArrays
            self.disvDict['centroids'] =centroids

        if save_spatial_index: # KMB
            self.spatialIndexDict['intervalNumber'] = intervalNumber
            self.spatialIndexDict['gridXarray'] = list(gridXarray)
            self.spatialIndexDict['gridYarray'] = list(gridYarray)
            self.spatialIndexDict['gridIndexList'] = gridIndexList

        return self.disvDict

    def save_properties(self,save_path,save_zip=False,compresslevel=9): # KMB, give option to save as zip
        if save_zip:
            with gzip.open(save_path,'wt',encoding='UTF-8',compresslevel=compresslevel) as outf:
                json.dump(self.disvDict,outf)
        else:
            with open(save_path, 'w') as outf:
                json.dump(self.disvDict, outf)


