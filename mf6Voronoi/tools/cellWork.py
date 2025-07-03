import rasterio
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon
from typing import Union

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

def getLayCellElevTupleFromElev(gwf,
                                interIx,
                                elevValue: Union[float,int],
                                geomPath):
    
    if isinstance(elevValue,(int,float)):
        print("You have inserted a fixed elevation")
    else:
        raise TypeError("Elevation value has to be a number or a list of numbers")
        
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

def getLayCellElevTupleFromObs(gwf,
                                interIx,
                                obsPath: str,
                                nameField: str,
                                elevField: str):
    
    obsDf = gpd.read_file(obsPath)
    insideCellsIds = []
    layCellTupleList = []
    nameList = []

    #model parameters
    nlay = gwf.modelgrid.nlay
    topBotm = gwf.modelgrid.top_botm

    #working with the cell ids
    #loop over the geometries to get the cellids
    for index, row in obsDf.iterrows():
        tempCellIds = interIx.intersect(row.geometry).cellids
        for cell in tempCellIds:
            insideCellsIds.append(cell)
        nameList.append(row[nameField])

    #working with the cell elevations
    for index, cell in enumerate(insideCellsIds):
        print('Working for cell %d'%cell)
        layerCell = False
        for lay in range(nlay):  # Loop through layers\n",
            if topBotm[lay+1, cell] < obsDf.iloc[index][elevField] <= topBotm[lay,cell]:
                layCellTupleList.append((lay,cell))
                print('Well screen elev of %.2f found at layer %d'%(obsDf.iloc[index][elevField],lay))
                layerCell = True
        if not layerCell:
            print('No layer was found for screen elevation: %.2f'%obsDf.iloc[index][elevField])

    return nameList, layCellTupleList 
    



