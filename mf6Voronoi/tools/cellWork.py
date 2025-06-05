import rasterio
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon

def getLayCellElevTupleFromRaster(gwf,interIx,rasterPath,geomPath):
    rasterSrc = rasterio.open(rasterPath)
    geomSrc = gpd.read_file(geomPath)
    insideCellsIds = []
    layCellTupleList = []
    cellElevList = []

    #model parameters
    nlay = gwf.modelgrid.nlay
    xCenter = gwf.modelgrid.xcellcenters
    yCenter = gwf.modelgrid.ycellcenters
    rasterElev = [elev[0] for elev in rasterSrc.sample(zip(xCenter,yCenter))] 
    topBotm = gwf.modelgrid.top_botm

    #working with the cell ids
    #loop over the geometries to get the cellids
    for index, row in geomSrc.iterrows():
        tempCellIds = interIx.intersect(row.geometry).cellids
        for cell in tempCellIds:
            insideCellsIds.append(cell)

    #working with the cell elevations and create laycell tuples
    for cell in insideCellsIds:
        #looping over elevations
        if topBotm[-1, cell] < rasterElev[cell] <= topBotm[0,cell]:
            cellElevList.append(rasterElev[cell])
        else: 
            print('The cell %d has a elevation of %.2f outside the model vertical domain'%(cell,rasterElev[cell]))
        #looping through layers
        for lay in range(nlay):  
            if topBotm[lay+1, cell] < rasterElev[cell] <= topBotm[lay,cell]:
                layCellTupleList.append((lay,cell))

    return layCellTupleList, cellElevList

def getLayCellElevTupleFromElev(gwf,interIx,elevValue,geomPath):
    geomSrc = gpd.read_file(geomPath)
    insideCellsIds = []
    layCellTupleList = []

    #model parameters
    nlay = gwf.modelgrid.nlay
    topBotm = gwf.modelgrid.top_botm

    #working with the cell ids
    #loop over the geometries to get the cellids
    for index, row in geomSrc.iterrows():
        tempCellIds = interIx.intersect(row.geometry).cellids
        for cell in tempCellIds:
            insideCellsIds.append(cell)

    #working with the cell elevations
    for cell in insideCellsIds:
        for lay in range(nlay):  # Loop through layers\n",
            if topBotm[lay+1, cell] < elevValue <= topBotm[lay,cell]:
                layCellTupleList.append((lay,cell))


    return layCellTupleList
        
    



