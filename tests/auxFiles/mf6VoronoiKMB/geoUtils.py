import os, io, fiona, folium
import geopandas as gpd
from scipy.spatial import Voronoi
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from .utils import readShpFromZip


# Select a colormap
cmap = plt.cm.viridis
         
def plotOrgDistPoints(vorMesh):
  limitXY = vorMesh.modelDis['limitGeometry']
  limitDf = gpd.GeoDataFrame([{'geometry':limitXY}])

  orgPoints = np.array(vorMesh.modelDis['vertexOrg'])

  tempDistPoints = []
  for key, value in vorMesh.modelDis['vertexDist'].items()  :
     tempDistPoints += value
  distPoints = np.array(tempDistPoints)
  #print(distPoints)

  fig, ax = plt.subplots(figsize=(36,24))
  ax.scatter(orgPoints[:,0], orgPoints[:,1], 
             label='Original', alpha=0.5, ec='crimson', fc='none')
  ax.scatter(distPoints[:,0], distPoints[:,1], s=2, marker='^',
             label='Distributed', alpha=0.5, ec='slateblue')
  limitDf.plot(ax=ax, label='Limit', alpha=0.5, ec='teal', fc='none', ls='-')
  for key, layer in vorMesh.discLayers.items():
     distLayersList = []
     for layerGeom in layer['layerGeoms']:
        distLayersList.append(layerGeom)
     layerDf = gpd.GeoDataFrame(geometry=distLayersList)
     #layerDf.plot(ax=ax, alpha=0.5, label=key.split('_')[0])
     layerDf.plot(ax=ax, alpha=0.5, label=None)
  ax.legend(loc='upper left')

  plt.show()

def plotCirclesPoints(vorMesh):
  #limitXY = vorMesh.modelDis['limitGeometry']
  circleUnionDf = gpd.GeoDataFrame([{'geometry':vorMesh.modelDis['circleUnion']}])
  polyPointList = vorMesh.modelDis['vertexBuffer']
  orgPoints = np.array(vorMesh.modelDis['vertexOrg'])

  tempDistPoints = []
  for key, value in vorMesh.modelDis['vertexDist'].items():
     tempDistPoints += value
  distPoints = np.array(tempDistPoints)
  #distPoints = np.array(vorMesh.modelDis['vertexDist'])
  polyPoints = np.array(polyPointList)

  fig, ax = plt.subplots(figsize=(36,24))
  ax.scatter(orgPoints[:,0], orgPoints[:,1], 
             label='Original', alpha=0.5, ec='crimson', fc='none')
  ax.scatter(distPoints[:,0], distPoints[:,1], s=5, marker='^',
             label='Distributed', alpha=0.5, ec='slateblue')
  ax.scatter(polyPoints[:,0], polyPoints[:,1], s=5, marker='^',
            label='CirclePoints', alpha=0.5, ec='tan')
  circleUnionDf.plot(ax=ax, label='circleUnion', alpha=0.5, ec='teal', fc='none', ls='-')
  #polyPointDf.plot(ax=ax, label='Limit', alpha=0.5, ec='tan', fc='none', ls='-')

  for key, layer in vorMesh.discLayers.items():
    distLayersList = []
    for layerGeom in layer['layerGeoms']:
      distLayersList.append(layerGeom)
    layerDf = gpd.GeoDataFrame(geometry=distLayersList)
    layerDf.plot(ax=ax, alpha=0.5, label=None)
  ax.legend(loc='upper left')

  plt.show()

def plotKeyList(vorMesh, keyList):
  fig, ax = plt.subplots(figsize=(36,24))
  for key in keyList:
    keyDf = gpd.GeoDataFrame([{'geometry':vorMesh.modelDis[key]}])
    keyDf.plot(ax=ax, label=key, alpha=0.5, ec='teal', fc='none', ls='-')
  
  ax.legend(loc='upper left')

  plt.show()
     
def bcGraph(self, ref_obj, ref_type, *args, **kwargs):
    project=self.getProject(self.request)

    limitLayer= LimitLayer.objects.get(project=project)
    meshLayer= VoronoiMesh.objects.get(project=project)

    gdfLimit=readShpFromZip(limitLayer.limitFile.path)
    gdfMesh=readShpFromZip(meshLayer.meshFile.path)

    if ref_type == "Recharge":
        gdfRef=readShpFromZip(ref_obj.rchFile.path)
    elif ref_type == "Evapotranspiration":
        gdfRef=readShpFromZip(ref_obj.evtFile.path)
    else:
        gdfRef=readShpFromZip(ref_obj.refinementFile.path)

    gdfRef['id'] = np.arange(0,gdfRef.shape[0])

    meshmap = folium.Map(max_zoom=20)

    gridGroup = folium.FeatureGroup(name='Grid')
    limGroup = folium.FeatureGroup(name='Limit Layer')
    refGroup = folium.FeatureGroup(name=ref_type)

    basemap2 = ('https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}')
    folium.WmsTileLayer(
      url=basemap2,
      layers=None,
      name='Google Sattelite',
      attr='Google',
      control=True,
      overlay=False,
    ).add_to(meshmap)

    gridStyle = {'fillColor': '#5AAACD',"weight": 1,"opacity": 0.65,"fillOpacity": 0}
    gridLayer = folium.GeoJson(data=gdfMesh["geometry"], style_function = lambda x: gridStyle)
    gridGroup.add_child(gridLayer)

    if gdfRef["geometry"][0].geom_type == 'Point':
        refLayer = folium.GeoJson(data=gdfRef,popup=folium.GeoJsonPopup(fields=['id']),
        marker = folium.CircleMarker(radius = 1, # Radius in metres
                                       weight = 6, #outline weight
                                       color = '#7B7C7E',
                                       fill_opacity = 1))
        refGroup.add_child(refLayer)
    else:
        refStyle = {'color': '#00CED1',"weight": 4,"opacity": 0.65,"fillOpacity": 0.2,"radius":3}
        refLayer = folium.GeoJson(data=gdfRef, style_function = lambda x: refStyle,
            popup=folium.GeoJsonPopup(fields=['id']))
        refGroup.add_child(refLayer)

    limStyle = {'color': '#ef4310',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
    limLayer = folium.GeoJson(data=gdfLimit["geometry"], style_function = lambda x: limStyle)
    limGroup.add_child(limLayer)

    meshmap.add_child(limGroup)
    meshmap.add_child(gridGroup)
    meshmap.add_child(refGroup)

    bPts = gdfMesh.to_crs("EPSG:4326").total_bounds
    xMinyMin = list(bPts[:2])[::-1]
    xMaxyMax = list(bPts[2:])[::-1]
    meshmap.fit_bounds([xMinyMin,xMaxyMax], padding=(1, 1))
    folium.LayerControl('topleft', collapsed= False).add_to(meshmap)

    meshmap = meshmap._repr_html_()

    return meshmap

def limitGraph(self, *args, **kwargs):
    project=self.getProject(self.request)
    limitLayer= LimitLayer.objects.get(project=project)

    gdfLimit=readShpFromZip(limitLayer.limitFile.path)

    meshmap = folium.Map(max_zoom=20)

    limGroup = folium.FeatureGroup(name='Limit Layer')

    basemap2 = ('https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}')
    folium.WmsTileLayer(
      url=basemap2,
      layers=None,
      name='Google Sattelite',
      attr='Google',
      control=True,
      overlay=False,
    ).add_to(meshmap)

    limStyle = {'color': '#ef4310',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
    limLayer = folium.GeoJson(data=gdfLimit["geometry"], style_function = lambda x: limStyle)
    limGroup.add_child(limLayer)

    meshmap.add_child(limGroup)

    bPts = gdfLimit.to_crs("EPSG:4326").total_bounds
    xMinyMin = list(bPts[:2])[::-1]
    xMaxyMax = list(bPts[2:])[::-1]
    meshmap.fit_bounds([xMinyMin,xMaxyMax], padding=(1, 1))
    folium.LayerControl('topleft', collapsed= False).add_to(meshmap)

    meshmap = meshmap._repr_html_()

    return meshmap

def meshGraph(self,mesh_obj, *args, **kwargs):
    project=self.getProject(self.request)
    limitLayer= LimitLayer.objects.get(project=project)
    refinementLayers=RefinementLayer.objects.filter(project=project)

    gdfLimit=readShpFromZip(limitLayer.limitFile.path)
    gdfRefDict = {}
    for index, refinemet in enumerate(refinementLayers):
        gdfRefDict[index] = readShpFromZip(refinementLayers[index].refinementFile.path)
    gdfMesh = readShpFromZip(mesh_obj.meshFile.path)

    meshmap = folium.Map(max_zoom=20)

    gridGroup = folium.FeatureGroup(name='Grid')
    limGroup = folium.FeatureGroup(name='Limit Layer')
    refGroup = folium.FeatureGroup(name='Refinement Layers')

    basemap2 = ('https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}')
    folium.WmsTileLayer(
      url=basemap2,
      layers=None,
      name='Google Sattelite',
      attr='Google',
      control=True,
      overlay=False,
    ).add_to(meshmap)

    gridStyle = {'fillColor': '#5AAACD',"weight": 1,"opacity": 0.65,"fillOpacity": 0}
    gridLayer = folium.GeoJson(data=gdfMesh["geometry"], style_function = lambda x: gridStyle)
    gridGroup.add_child(gridLayer)

    meshmap.add_child(gridGroup)

    for ref in gdfRefDict:
        if gdfRefDict[ref]["geometry"][0].geom_type == 'Point':
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"],
                marker = folium.CircleMarker(radius = 1, # Radius in metres
                                           weight = 6, #outline weight
                                           color = '#7B7C7E',
                                           fill_opacity = 1))
            refGroup.add_child(refLayer)
        else:
            refStyle = {'color': '#00CED1',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"], style_function = lambda x: refStyle)
            refGroup.add_child(refLayer)

    limStyle = {'color': '#ef4310',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
    limLayer = folium.GeoJson(data=gdfLimit["geometry"], style_function = lambda x: limStyle)
    limGroup.add_child(limLayer)

    meshmap.add_child(limGroup)
    meshmap.add_child(refGroup)

    bPts = gdfMesh.to_crs("EPSG:4326").total_bounds
    xMinyMin = list(bPts[:2])[::-1]
    xMaxyMax = list(bPts[2:])[::-1]
    meshmap.fit_bounds([xMinyMin,xMaxyMax], padding=(1, 1))
    folium.LayerControl('topleft', collapsed= False).add_to(meshmap)

    meshmap = meshmap._repr_html_()

    return meshmap

def pointCloudGraph(self, meshForm, vorMesh, *args, **kwargs):
    project=self.getProject(self.request)
    limitLayer= LimitLayer.objects.get(project=project)
    refinementLayers=RefinementLayer.objects.filter(project=project)

    outNamePointCloud = 'voronoiPointCloud'
    outNameCircleUnion = 'voronoiCircleUnion'

    outNamePointCloudDir = os.path.join(meshForm.getAbsDir,'pointCloudTemp')

    outTempPointCloudPath = os.path.join(outNamePointCloudDir,outNamePointCloud+'.shp')
    outCircleUnionPath = os.path.join(outNamePointCloudDir,'voronoiCircleUnion.shp')
    outVertexBufferPath = os.path.join(outNamePointCloudDir,'vertexBuffer.shp')
    outInteriorsPath = os.path.join(outNamePointCloudDir,'interiors.shp')

    if not os.path.isdir(outNamePointCloudDir):
        os.makedirs(outNamePointCloudDir, exist_ok=True)

    vorMesh.getPointsAsShp('vertexTotal',outTempPointCloudPath)
    vorMesh.getPolyAsShp('circleUnion',outCircleUnionPath)
    vorMesh.getPointsAsShp('vertexBuffer',outVertexBufferPath)
    vorMesh.getPolyAsShp('interiors',outInteriorsPath)

    gdfLimit=readShpFromZip(limitLayer.limitFile.path)
    gdfRefDict = {}
    for index, refinemet in enumerate(refinementLayers):
        gdfRefDict[index] = readShpFromZip(refinementLayers[index].refinementFile.path)
    gdfMesh=readShpFromZip(meshForm.meshFile.path)

    gdfPointCloud=gpd.read_file(outTempPointCloudPath)
    gdfCircleUnion=gpd.read_file(outCircleUnionPath)
    gdfVertexBuffer=gpd.read_file(outVertexBufferPath)
    gdfInteriors=gpd.read_file(outInteriorsPath)

    #add point cloud map
    gdfPointCloud=gdfPointCloud.to_crs("EPSG:4326")
    gdfPtClBds=gdfPointCloud.total_bounds
    ptClMap = folium.Map(max_zoom=20)

    gridGroup = folium.FeatureGroup(name='Grid')
    refGroup = folium.FeatureGroup(name='Refinement Layers')
    limGroup = folium.FeatureGroup(name='Limit Layer')
    pointGroup = folium.FeatureGroup(name='Point Cloud')
    circleUnionGroup = folium.FeatureGroup(name='Circle Union')
    vertexBufferGroup = folium.FeatureGroup(name='Vertex Buffer')
    interiorsGroup = folium.FeatureGroup(name='Interiors')

    pointLayer = folium.GeoJson(data=gdfPointCloud["geometry"],
        marker = folium.CircleMarker(radius = 1, # Radius in metres
                                   weight = 2, #outline weight
                                   fill_color = '#5AAACD',
                                   fill_opacity = 1))
    pointGroup.add_child(pointLayer)

    gridStyle = {'fillColor': '#5AAACD',"weight": 1,"opacity": 0.65,"fillOpacity": 0}
    gridLayer = folium.GeoJson(data=gdfMesh["geometry"], style_function = lambda x: gridStyle)
    gridGroup.add_child(gridLayer)

    for ref in gdfRefDict:
        if gdfRefDict[ref]["geometry"][0].geom_type == 'Point':
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"],
                marker = folium.CircleMarker(radius = 1, # Radius in metres
                                           weight = 6, #outline weight
                                           color = '#7B7C7E',
                                           fill_opacity = 1))
            refGroup.add_child(refLayer)
        else:
            refStyle = {'color': '#00CED1',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"], style_function = lambda x: refStyle)
            refGroup.add_child(refLayer)

    limStyle = {'color': '#ef4310',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
    limLayer = folium.GeoJson(data=gdfLimit["geometry"], style_function = lambda x: limStyle)
    limGroup.add_child(limLayer)

    circleUnionStyle = {'color': '#80BCD8',"weight": 4,"opacity": 0.65,"fillOpacity": 0}
    circleUnionLayer = folium.GeoJson(data=gdfCircleUnion["geometry"], style_function = lambda x: circleUnionStyle)
    circleUnionGroup.add_child(circleUnionLayer)

    vertexBufferLayer = folium.GeoJson(data=gdfVertexBuffer["geometry"],
        marker = folium.CircleMarker(radius = 4,
                                    weight = .5,
                                   color = '#80BCD8',
                                   fill_opacity = 0.5))
    vertexBufferGroup.add_child(vertexBufferLayer)

    if not gdfInteriors.empty:
        interiorsStyle = {'color': '#B22222',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
        interiorsLayer = folium.GeoJson(data=gdfInteriors["geometry"], style_function = lambda x: interiorsStyle)
        interiorsGroup.add_child(interiorsLayer)

    ptClMap.add_child(gridGroup)
    ptClMap.add_child(limGroup)
    ptClMap.add_child(refGroup)
    ptClMap.add_child(pointGroup)
    ptClMap.add_child(circleUnionGroup)
    ptClMap.add_child(vertexBufferGroup)
    if not gdfInteriors.empty:
        ptClMap.add_child(interiorsGroup)

    bPts = gdfMesh.to_crs("EPSG:4326").total_bounds
    xMinyMin = list(bPts[:2])[::-1]
    xMaxyMax = list(bPts[2:])[::-1]
    ptClMap.fit_bounds([xMinyMin,xMaxyMax], padding=(1, 1))
    folium.LayerControl('topleft', collapsed= False).add_to(ptClMap)
    ptClMap=ptClMap._repr_html_() #up

    return ptClMap

def waterTableGraph(self, c_map, Line, *args, **kwargs):
    project=self.getProject(self.request)
    limitLayer= LimitLayer.objects.get(project=project)
    refinementLayers=RefinementLayer.objects.filter(project=project)
    meshLayer = VoronoiMesh.objects.get(project=project)

    gdfLimit=readShpFromZip(limitLayer.limitFile.path)
    gdfRefDict = {}
    for index, refinemet in enumerate(refinementLayers):
        gdfRefDict[index] = readShpFromZip(refinementLayers[index].refinementFile.path)
    gdfMesh = readShpFromZip(meshLayer.meshFile.path)

    meshmap = folium.Map(max_zoom=20)

    gridGroup = folium.FeatureGroup(name='Grid')
    limGroup = folium.FeatureGroup(name='Limit Layer')
    refGroup = folium.FeatureGroup(name='Refinement Layers')
    contourGroup = folium.FeatureGroup(name='Contours')

    basemap2 = ('https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={z}')
    folium.WmsTileLayer(
      url=basemap2,
      layers=None,
      name='Google Sattelite',
      attr='Google',
      control=True,
      overlay=False,
    ).add_to(meshmap)

    gridStyle = {'fillColor': '#5AAACD',"weight": 1,"opacity": 0.65,"fillOpacity": 0}
    gridLayer = folium.GeoJson(data=gdfMesh["geometry"], style_function = lambda x: gridStyle)
    gridGroup.add_child(gridLayer)

    meshmap.add_child(gridGroup)

    for ref in gdfRefDict:
        if gdfRefDict[ref]["geometry"][0].geom_type == 'Point':
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"],
                marker = folium.CircleMarker(radius = 1, # Radius in metres
                                           weight = 6, #outline weight
                                           color = '#7B7C7E',
                                           fill_opacity = 1))
            refGroup.add_child(refLayer)
        else:
            refStyle = {'color': '#00CED1',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
            refLayer = folium.GeoJson(data=gdfRefDict[ref]["geometry"], style_function = lambda x: refStyle)
            refGroup.add_child(refLayer)

    limStyle = {'color': '#ef4310',"weight": 4,"opacity": 0.65,"fillOpacity": 0,"radius":3}
    limLayer = folium.GeoJson(data=gdfLimit["geometry"], style_function = lambda x: limStyle)
    limGroup.add_child(limLayer)

    meshmap.add_child(limGroup)
    meshmap.add_child(refGroup)

    if not c_map.empty:
        red = Color("#97FFFF")
        colors = list(red.range_to(Color("#104E8B"),len(Line)))
        for idx, r in c_map.to_crs("EPSG:4326").iterrows():
             heads_to_map(r,contourGroup,colors[idx])
        meshmap.add_child(contourGroup)

    bPts = gdfMesh.to_crs("EPSG:4326").total_bounds
    xMinyMin = list(bPts[:2])[::-1]
    xMaxyMax = list(bPts[2:])[::-1]
    meshmap.fit_bounds([xMinyMin,xMaxyMax], padding=(1, 1))
    folium.LayerControl('topleft', collapsed= False).add_to(meshmap)

    meshmap = meshmap._repr_html_()

    return meshmap

class Field():
  '''
  Create a Voronoi map that can be used to run Lloyd
  relaxation on an array of 2D points. For background,
  see: https://en.wikipedia.org/wiki/Lloyd%27s_algorithm
  '''

  def __init__(self, *args, **kwargs):
    '''
    Store the points and bounding box of the points to which
    Lloyd relaxation will be applied.
    @param np.array `arr`: a numpy array with shape n, 2, where n
      is the number of 2D points to be moved
    @param float `epsilon`: the delta between the input point
      domain and the pseudo-points used to constrain the points
    '''
    arr = args[0]
    if not isinstance(arr, np.ndarray) or arr.shape[1] != 2:
      raise Exception('Please provide a numpy array with shape n,2')
    self.points = arr
    # find the bounding box of the input data
    self.domains = self.get_domains(arr)
    # ensure no two points have the exact same coords
    self.bb_points = self.get_bb_points(arr)
    self.constrain = kwargs.get('constrain', True)
    self.build_voronoi()

  def constrain_points(self):
    '''
    Update any points that have drifted beyond the boundaries of this space
    '''
    for point in self.points:
      if point[0] < self.domains['x']['min']: point[0] = self.domains['x']['min']
      if point[0] > self.domains['x']['max']: point[0] = self.domains['x']['max']
      if point[1] < self.domains['y']['min']: point[1] = self.domains['y']['min']
      if point[1] > self.domains['y']['max']: point[1] = self.domains['y']['max']


  def get_domains(self, arr):
    '''
    Return an object with the x, y domains of `arr`
    '''
    x = arr[:, 0]
    y = arr[:, 1]
    return {
      'x': {
        'min': min(x),
        'max': max(x),
      },
      'y': {
        'min': min(y),
        'max': max(y),
      }
    }


  def get_bb_points(self, arr):
    '''
    Given an array of 2D points, return the four vertex bounding box
    '''
    return np.array([
      [self.domains['x']['min'], self.domains['y']['min']],
      [self.domains['x']['max'], self.domains['y']['min']],
      [self.domains['x']['min'], self.domains['y']['max']],
      [self.domains['x']['max'], self.domains['y']['max']],
    ])


  def build_voronoi(self):
    '''
    Build a voronoi map from self.points. For background on
    self.voronoi attributes, see: https://docs.scipy.org/doc/scipy/
      reference/generated/scipy.spatial.Voronoi.html
    '''
    # build the voronoi tessellation map
    self.voronoi = Voronoi(self.points, qhull_options='Qbb Qc Qx')

    # constrain voronoi vertices within bounding box
    if self.constrain:
      for idx, vertex in enumerate(self.voronoi.vertices):
        x, y = vertex
        if x < self.domains['x']['min']:
          self.voronoi.vertices[idx][0] = self.domains['x']['min']
        if x > self.domains['x']['max']:
          self.voronoi.vertices[idx][0] = self.domains['x']['max']
        if y < self.domains['y']['min']:
          self.voronoi.vertices[idx][1] = self.domains['y']['min']
        if y > self.domains['y']['max']:
          self.voronoi.vertices[idx][1] = self.domains['y']['max']

  def find_centroid(self, vertices):
    '''
    Find the centroid of a Voroni region described by `vertices`,
    and return a np array with the x and y coords of that centroid.
    The equation for the method used here to find the centroid of a
    2D polygon is given here: https://en.wikipedia.org/wiki/
      Centroid#Of_a_polygon
    @params: np.array `vertices` a numpy array with shape n,2
    @returns np.array a numpy array that defines the x, y coords
      of the centroid described by `vertices`
    '''
    area = 0
    centroid_x = 0
    centroid_y = 0
    for i in range(len(vertices)-1):
      step = (vertices[i  , 0] * vertices[i+1, 1]) - \
             (vertices[i+1, 0] * vertices[i  , 1])
      area += step
      centroid_x += (vertices[i, 0] + vertices[i+1, 0]) * step
      centroid_y += (vertices[i, 1] + vertices[i+1, 1]) * step
    area /= 2
    # prevent division by zero - equation linked above
    if area == 0: area += 0.0000001
    centroid_x = (1.0/(6.0*area)) * centroid_x
    centroid_y = (1.0/(6.0*area)) * centroid_y
    # prevent centroids from escaping bounding box
    if self.constrain:
      if centroid_x < self.domains['x']['min']: centroid_x = self.domains['x']['min']
      if centroid_x > self.domains['x']['max']: centroid_x = self.domains['x']['max']
      if centroid_y < self.domains['y']['min']: centroid_y = self.domains['y']['min']
      if centroid_y > self.domains['y']['max']: centroid_y = self.domains['y']['max']
    return np.array([centroid_x, centroid_y])


  def relax(self):
    '''
    Moves each point to the centroid of its cell in the voronoi
    map to "relax" the points (i.e. jitter the points so as
    to spread them out within the space).
    '''
    centroids = []
    for idx in self.voronoi.point_region:
      # the region is a series of indices into self.voronoi.vertices
      # remove point at infinity, designated by index -1
      region = [i for i in self.voronoi.regions[idx] if i != -1]
      # enclose the polygon
      region = region + [region[0]]
      # get the vertices for this regioncd
      verts = self.voronoi.vertices[region]
      # find the centroid of those vertices
      centroids.append(self.find_centroid(verts))
    self.points = np.array(centroids)
    self.constrain_points()
    self.build_voronoi()


  def get_points(self):
    '''
    Return the input points in the new projected positions
    @returns np.array a numpy array that contains the same number
      of observations in the input points, in identical order
    '''
    return self.points