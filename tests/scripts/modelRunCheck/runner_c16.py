from datetime import datetime
from pathlib import Path
import os
import sys
import json
import flopy
import rasterio
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import MultiLineString, Point
from mf6Voronoi.tools.cellWork import (
    getLayCellElevTupleFromRaster,
    getLayCellElevTupleFromElev,
    getLayCellElevTupleFromObs,
)

# --------------------------
# Carga de Configuración de Rutas
# --------------------------

jsonPath = '../../json/modelCases.json'

with open(jsonPath, "r") as f:
    configDict = json.load(f)

# --------------------------
# Argumentos de CLI
# --------------------------

#processing argument
if len(sys.argv) < 2:
    sys.exit(
        "Usage: python script.py <caseName> <modelNumber>\n\n"
        "Arguments:\n"
        "  caseName   : Name of the test case\n"
        "  runType    : normalRun | parallelRun\n"
        "  modelMumber  : int\n"
    )

caseName = sys.argv[1]
runType = sys.argv[2]

if runType == 'normalRun':
    print('Computing with normal CPU')
elif runType == 'parallelRun':
    print('Computing with parallelized computing')
else:
    sys.exit(f"Invalid execution mode: '{runType}'. Expected 'normalRun' or 'parallelRun'.")

modelNumber = sys.argv[3]

# --------------------------
# Defining folders
# --------------------------

caseGenDataFolder = configDict["base_dirs"]["case_gen_data"]
meshGenDataFolder = configDict["base_dirs"]["mesh_gen_data"]

# Get today's current date Format date: %d (day), %b (short month), %y (2-digit year)
todayStr = datetime.now().strftime("%d%b%y")

meshName = configDict["mesh_name"][modelNumber]

meshCaseDir = os.path.join(meshGenDataFolder, f"{caseName}_{todayStr}", sys.argv[2], meshName)
jsonDir = os.path.join(meshCaseDir, 'json')
rstDir = os.path.join(caseGenDataFolder,meshName,'rst')
shpDir = os.path.join(caseGenDataFolder,meshName,'shp')
modelWs = os.path.join(meshCaseDir, 'model')
os.makedirs(modelWs, exist_ok=True)
mfBin = os.path.join(caseGenDataFolder,'modflowBin','mf6')

# --------------------------
# Model folder
# --------------------------

# open the json file
with open(os.path.join(jsonDir,'disvDict.json')) as file: ## Org
    gridProps = json.load(file) ## Org

cell2d = gridProps['cell2d']           #cellid, cell centroid xy, vertex number and vertex id list
vertices = gridProps['vertices']       #vertex id and xy coordinates
ncpl = gridProps['ncpl']               #number of cells per layer
nvert = gridProps['nvert']             #number of verts
centroids=gridProps['centroids']   

##############hasta aqui copiar#############################

#Extract dem values for each centroid of the voronois
src = rasterio.open(os.path.join(rstDir,'n33w111_wgs84_int32_20m.tif'))  ## Org
elevation=[x for x in src.sample(centroids)] ## Org

finalSrc = rasterio.open(os.path.join(rstDir,'demWithReservoir.tif'))
finalElevation=[x for x in finalSrc.sample(centroids)]

#fix elevation to stay allways below final elevation
for index, elev in enumerate(elevation):
    if elev > finalElevation[index]:
        elevation[index] = finalElevation[index] 
    else:
        pass

nlay = 14 #10   ## Org

#original and final surface
orgMtop = np.array([elev[0] for i,elev in enumerate(elevation)]) ## Bottom of layer 1
mtop = np.array([elev[0] for i,elev in enumerate(finalElevation)])

zbot=np.zeros((nlay,ncpl)) ## Org

AcuifInf_Bottom = 700 ## Org

#bottom original layer 1
botLay1 = np.copy(AcuifInf_Bottom + (0.95 * (orgMtop - AcuifInf_Bottom)))

#added layers
zbot[0,] = botLay1 + 0.8 * (mtop - botLay1) 
zbot[1,] = botLay1 + 0.6 * (mtop - botLay1) 
zbot[2,] = botLay1 + 0.4 * (mtop - botLay1) 
zbot[3,] = botLay1 + 0.2 * (mtop - botLay1) 

#original layers but different index
zbot[4,] = AcuifInf_Bottom + (0.95 * (orgMtop - AcuifInf_Bottom)) ## <==== updated
zbot[5,] = AcuifInf_Bottom + (0.90 * (orgMtop - AcuifInf_Bottom)) ## <==== updated
zbot[6,] = AcuifInf_Bottom + (0.85 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 85%
zbot[7,] = AcuifInf_Bottom + (0.78 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 
zbot[8,] = AcuifInf_Bottom + (0.71 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 
zbot[9,] = AcuifInf_Bottom + (0.64 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 
zbot[10,] = AcuifInf_Bottom + (0.57 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 
zbot[11,] = AcuifInf_Bottom + (0.50 * (orgMtop - AcuifInf_Bottom)) ## <==== updated 50%
zbot[12,] = AcuifInf_Bottom + (0.25 * (orgMtop - AcuifInf_Bottom)) ## <==== updated
zbot[13,] = AcuifInf_Bottom ## <==== updated

# create simulation
simName = 'mf6Sim' ## Org
modelName = 'mf6Model' ## Org
sim = flopy.mf6.MFSimulation(sim_name=modelName, version='mf6', ## Org
                             exe_name=mfBin, ## Org
                             sim_ws=modelWs) ## Org

# create tdis package
tdis_rc = [(1000.0, 1, 1.0)] ## Org
tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='SECONDS', ## Org
                             perioddata=tdis_rc) ## Org

# create gwf model
gwf = flopy.mf6.ModflowGwf(sim, ## Org
                           modelname=modelName, ## Org
                           save_flows=True, ## Org
                           newtonoptions="NEWTON UNDER_RELAXATION") ## Org

# create iterative model solution and register the gwf model with it
ims = flopy.mf6.ModflowIms(sim, ## Org
                           complexity='COMPLEX', ## Org
                           outer_maximum=150, ## Org
                           inner_maximum=30, ## Org
                           linear_acceleration='BICGSTAB') ## Org
sim.register_ims_package(ims,[modelName]) ## Org

# disv
disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, ## Org
                                top=mtop, botm=zbot, ## Org
                                nvert=nvert, vertices=vertices, ## Org
                                cell2d=cell2d) ## Org

disv.top.plot(figsize=(12,8), alpha=0.8) ## Org
plt.show()

crossSection = gpd.read_file(os.path.join(shpDir,'crossSectionReservoir.shp')) ## Org
sectionLine =list(crossSection.iloc[0].geometry.coords) ## Org

fig, ax = plt.subplots(figsize=(12,8)) ## Org
modelxsect = flopy.plot.PlotCrossSection(model=gwf, line={'Line': sectionLine}) ## Org
linecollection = modelxsect.plot_grid(lw=0.5) ## Org
ax.set_ylim(900,1300)
ax.set_xlim(1500,2500)
plt.show()

# initial conditions
waterTableSrc = rasterio.open(os.path.join(rstDir,'waterTable.tif'))
regWaterTable = [x for x in waterTableSrc.sample(centroids)]
ic = flopy.mf6.ModflowGwfic(gwf, strt=np.stack([regWaterTable for i in range(nlay)])) ## Org

ic.plot(mflay=0)
plt.show()

kxOnes = np.ones([nlay,ncpl])
#original matrix of Kx values
Kx =[4E-4, 4E-4, 4E-4, 4E-4, 4E-4, 5E-5, 3E-6, 3E-6, 2.5E-6, 2.5E-6, 2.5E-6, 1E-6, 9E-7, 5E-7] ## Org
layerValues = np.array(Kx).reshape(nlay,1)
kxArray = kxOnes*layerValues

icelltype = [1,1,1,1,1,1,1,1,1,1,0,0,0,0] ## Org

# Define intersection object
interIx = flopy.utils.gridintersect.GridIntersect(gwf.modelgrid) ## Org

#reservoir apply kx to cells in layers from 1 to 5 
resDf = gpd.read_file(os.path.join(shpDir,'reservoirEnvelope.shp'))
resCells = interIx.intersect(resDf.iloc[0].geometry).cellids
for cell in resCells:
    for layer in range(5):
        kxArray[layer,cell] = 1

#dam apply kx to cells in layers from 1 to 5 from cells downstream
damDf = gpd.read_file(os.path.join(shpDir,'damEnvelopeDownstream.shp'))
damCells = interIx.intersect(damDf.iloc[0].geometry).cellids
for cell in damCells:
    for layer in range(6):
        kxArray[layer,cell] = 1e-8

#apply dam kx to layers upstream
layCellTupleList, cellElevList = getLayCellElevTupleFromRaster(gwf,
                                                               interIx,
                                                               os.path.join(rstDir,'damElev.tif'),
                                                               os.path.join(shpDir,'damEnvelopeUpstream.shp')) ## # <=====
#looping over celltuples to get cells from intersected layer to layer 5
for layCell in layCellTupleList:
    #filter layers below layer 5
    if layCell[0] <=5:
        for extendedLayer in range(layCell[0],6):
            kxArray[extendedLayer,layCell[1]] = 1e-8

# node property flow
npf = flopy.mf6.ModflowGwfnpf(gwf, ## Org
                              save_specific_discharge=True, ## Org
                              icelltype=icelltype, ## Org
                              k=kxArray) ## Org

## <==== inserted
crossSection = gpd.read_file(os.path.join(shpDir,'crossSectionReservoir.shp')) ## <==== inserted
sectionLine =list(crossSection.iloc[0].geometry.coords) ## <==== inserted

fig, ax = plt.subplots(figsize=(12,8)) ## <==== inserted
modelxsect = flopy.plot.PlotCrossSection(model=gwf, line={'Line': sectionLine}) ## <==== inserted
linecollection = modelxsect.plot_grid(lw=0.5) ## <==== inserted
modelxsect.plot_array(np.log(npf.k.array), alpha=0.5) ## <====== inserted
ax.grid() ## <==== inserted
ax.set_ylim(900,1300)
ax.set_xlim(1600,2600)
plt.show()

# define storage and transient stress periods
sto = flopy.mf6.ModflowGwfsto(gwf, ## Org
                              iconvert=1, ## Org
                              steady_state={ ## Org
                                0:True, ## Org
                              } ## Org
                              ) ## Org

rchr = 0.15/365/86400 ## Org
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rchr) ## Org
evtr = 1.2/365/86400 ## Org
evt = flopy.mf6.ModflowGwfevta(gwf,ievt=1,surface=mtop,rate=evtr,depth=1.0) ## Org

# Define intersection object
interIx = flopy.utils.gridintersect.GridIntersect(gwf.modelgrid) ## Org

#open the river shapefile
rivers =gpd.read_file(os.path.join(shpDir,'river_basin_clip.shp')) ## Org
list_rivers=[] ## Org
for i in range(rivers.shape[0]): ## Org
    list_rivers.append(rivers['geometry'].loc[i]) ## Org
    
riverMls = MultiLineString(lines=list_rivers) ## Org

#intersec rivers with our grid
riverCells=interIx.intersect(riverMls).cellids ## Org
riverCells[:10] ## Org

#river package
riverSpd = {} ## Org
riverSpd[0] = [] ## Org
for cell in riverCells: ## Org
    riverSpd[0].append([(0,cell),mtop[cell],0.01]) ## Org
riv = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=riverSpd) ## Org

#river plot
riv.plot(mflay=0) ## Org
plt.show()

#general boundary condition

#regional flow
layCellTupleList, cellElevList = getLayCellElevTupleFromRaster(gwf,
                                                               interIx,
                                                               os.path.join(rstDir,'waterTable.tif'),
                                                               os.path.join(shpDir,'regionalFlow.shp')) ## # <===== Inserted
regionalSpd = {} ## # <===== Inserted
regionalSpd[0] = [] ## # <===== Inserted
for index, layCellTuple in enumerate(layCellTupleList): ## Org
    regionalSpd[0].append([layCellTuple,cellElevList[index],0.01,'regionalFlow']) # <===== Inserted

layCellTupleList = getLayCellElevTupleFromElev(gwf,
                                               interIx,
                                               1220,
                                               os.path.join(shpDir,'reservoirEnvelope.shp'))
for layCellTuple in layCellTupleList:
    regionalSpd[0].append([layCellTuple, 1220, 10, 'reservoir'])

ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=regionalSpd, boundnames=True) ## <==== modified

# Observation package for Drain
obsDict = { # <===== Inserted 
    "{}.ghb.obs.csv".format(modelName): [ # <===== Inserted 
        ("regional_flow", "ghb", "regionalFlow"), # <===== Inserted 
        ("reservoir", "ghb", "reservoir") # <===== Inserted 
    ] # <===== Inserted 
} # <===== Inserted 

# Attach observation package to DRN package
ghb.obs.initialize( # <===== Inserted 
    filename=gwf.name+".ghb.obs", # <===== Inserted 
    digits=10, # <===== Inserted 
    print_input=True, # <===== Inserted 
    continuous=obsDict # <===== Inserted 
) # <===== Inserted

#ghb plot
ghb.plot(mflay=0, kper=0) # <===== Inserted
plt.show()

#oc
head_filerecord = f"{gwf.name}.hds" ## Org
budget_filerecord = f"{gwf.name}.cbc" ## Org
oc = flopy.mf6.ModflowGwfoc(gwf, ## Org
                            head_filerecord=head_filerecord, ## Org
                            budget_filerecord = budget_filerecord, ## Org
                            saverecord=[("HEAD", "LAST"),("BUDGET","LAST")]) ## Org

# Run the simulation
sim.write_simulation() ## Org
success, buff = sim.run_simulation() ## Orgclear