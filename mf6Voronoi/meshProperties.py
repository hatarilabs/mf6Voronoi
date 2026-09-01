import json
import numpy as np
import fiona
from shapely.geometry import Polygon
from tqdm import tqdm
from .utils import isRunningInJupyter, printBannerHtml, printBannerText

class meshShape:
	def __init__(self,path_to_mesh_shp,path_to_mesh_shp=None,mesh_df=None): # KMB addition
        self.mesh=path_to_mesh_shp
        self.mesh_df = mesh_df # KMB
        self.disvDict = {}
        self.spatialIndexDict = {}

	def get_gridprops_disv(self, intervalNumber=10,save_spatial_index=True, verbose=False):

		# //------------------------KMB 
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

