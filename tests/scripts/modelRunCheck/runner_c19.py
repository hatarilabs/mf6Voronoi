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

#### Part 2b: Model construction and simulation
#Extract dem values for each centroid of the voronois
from shapely.geometry import Point

src = rasterio.open(os.path.join(rstDir, 'modelDem.tif'))  ## Org
seaDf = gpd.read_file(os.path.join(shpDir, 'seaBeach.shp')) #<==== updated

elevation = [] #initialize list

for centroid in centroids: #fix values of raster inside the sea polygon
    centPoint = Point(centroid) #get geometry of the cell centroid

    if centPoint.within(seaDf.iloc[0].geometry): #if it is inside the sea
        elevation.append(0)
    else:
        elevation += [x[0].item() for x in src.sample([centroid])] #assing raster value

nlay = 8 ## Org

mtop=np.array(elevation) #[elev[0] for i,elev in enumerate(elevation)]) ## Org
zbot=np.zeros((nlay,ncpl)) ## Org


AcuifInf_Bottom = -150# <==== update
zbot[0,] = AcuifInf_Bottom + (0.875 * (mtop - AcuifInf_Bottom)) ## Org
zbot[1,] = AcuifInf_Bottom + (0.75 * (mtop - AcuifInf_Bottom)) ## Org
zbot[2,] = AcuifInf_Bottom + (0.625 * (mtop - AcuifInf_Bottom)) ## Org
zbot[3,] = AcuifInf_Bottom + (0.5 * (mtop - AcuifInf_Bottom)) ## Org
zbot[4,] = AcuifInf_Bottom + (0.375 * (mtop - AcuifInf_Bottom)) ## Org
zbot[5,] = AcuifInf_Bottom + (0.25 * (mtop - AcuifInf_Bottom)) ## Org
zbot[6,] = AcuifInf_Bottom + (0.125 * (mtop - AcuifInf_Bottom)) ## Org
zbot[7,] = AcuifInf_Bottom ## Org

#### Create simulation and model

# create simulation
simName = 'mf6Sim' ## Org
modelName = 'mf6Model' ## Org
sim = flopy.mf6.MFSimulation(sim_name=simName, version='mf6', ## Org
                             exe_name=mfBin, ## Org
                             continue_=True,
                             sim_ws=modelWs) ## Org

# create tdis package
tdis_rc = [(86400.0*365*50, 1, 1.0)] + [(86400*365*30, 1, 1.0)] ## 30 years, 15 years, 15 years
#tdis_rc = [(86400.0*365*50, 1, 1.0)] + [(86400*365*30, 1, 1.0) for level in range(2)] ## 30 years, 15 years, 15 years
print(tdis_rc[:3]) ## Org

tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='SECONDS', ## Org
                             perioddata=tdis_rc, ## Org
                            nper=2) ## Org

# create gwf model
gwf = flopy.mf6.ModflowGwf(sim, ## Org
                           modelname=modelName, ## Org
                           save_flows=True, ## Org
                           newtonoptions="NEWTON UNDER_RELAXATION") ## Org

# create iterative model solution and register the gwf model with it
imsGwf = flopy.mf6.ModflowIms(sim, ## Org
                              pname='ims_gwf',
                           complexity='COMPLEX', ## Org
                           outer_maximum=150, ## Org
                           inner_maximum=50, ## Org
                           outer_dvclose=0.1, ## Org
                           inner_dvclose=0.0001, ## Org
                           backtracking_number=20, ## Org
                           linear_acceleration='BICGSTAB') ## Org
sim.register_ims_package(imsGwf,[modelName]) ## Org

# disv
disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, ## Org
                                top=mtop, botm=zbot, ## Org
                                nvert=nvert, vertices=vertices, ## Org
                                cell2d=cell2d) ## Org

plt.figure(figsize=(12, 8))
disv.top.plot(figsize=(12,8), alpha=0.8) ## Org
plt.show()

crossSection = gpd.read_file(os.path.join(shpDir,'crossSection.shp')) ## Org
sectionLine =list(crossSection.iloc[0].geometry.coords) ## Org

fig, ax = plt.subplots(figsize=(12,8)) ## Org
modelxsect = flopy.plot.PlotCrossSection(model=gwf, line={'Line': sectionLine}) ## Org
linecollection = modelxsect.plot_grid(lw=0.5) ## Org
ax.set_ylim(-150,100)
ax.grid() ## Org
plt.show()

# initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=np.stack([mtop for i in range(nlay)])) ## Org
Kx =[4E-3 for x in range(3)] + [1E-3 for x in range(3)] + [4E-4 for x in range(2)] ## <=== updated
icelltype = [1 for x in range(5)] + [0 for x in range(nlay - 5)] ## Org

# node property flow
npf = flopy.mf6.ModflowGwfnpf(gwf, ## Org
                              save_specific_discharge=True, ## Org
                              icelltype=icelltype, ## Org
                              k=Kx, ## Org
                              k33=Kx) ## <== updated

# define storage and transient stress periods
sto = flopy.mf6.ModflowGwfsto(gwf, ## Org
                              iconvert=1, ## Org
                              steady_state={ ## Org
                                0:True, ## Org
                              },
                              transient={
                                  1:True, ## Org
                                  2:True, ## Org
                              },
                              ss=1e-06,
                              sy=0.001,
                              ) ## Org

rchr = 0.2/365/86400 ## Org
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rchr) ## Org
evtr = 1.2/365/86400 ## Org
evt = flopy.mf6.ModflowGwfevta(gwf,ievt=1,surface=mtop,rate=evtr,depth=1.0) ## Org

# Define intersection object
interIx = flopy.utils.gridintersect.GridIntersect(gwf.modelgrid) ## Org

#general boundary condition

ghbSpd = {} ## # <===== Inserted
ghbSpd[0] = [] ## # <===== Inserted

# #regional flow
# layCellTupleList, cellElevList = getLayCellElevTupleFromRaster(gwf,
#                                                                interIx,
#                                                                '../rst/waterTable.tif',
#                                                                '../shp/regGhb.shp') ## # <===== Inserted

# for index, layCellTuple in enumerate(layCellTupleList): ## Org
#     ghbSpd[0].append([layCellTuple,cellElevList[index],0.01, 0, 'regflow']) # <===== Inserted


layCellTupleList = getLayCellElevTupleFromElev(gwf,
                                               interIx,
                                               0,
                                               os.path.join(shpDir,'seaGhb.shp'))
for layCellTuple in layCellTupleList:
    ghbSpd[0].append([layCellTuple, 0, 0.20, 0, 'sea'])

ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghbSpd, auxiliary=['CONCENTRATION'], boundnames=True) ## <==== modified

# Observation package for Drain
obsDict = { # <===== Inserted 
    "{}.ghb.obs.csv".format(modelName): [ # <===== Inserted 
        ("regflow", "ghb", "regionalFlow"), # <===== Inserted 
        ("sea", "ghb", "sea") # <===== Inserted 
    ] # <===== Inserted 
} # <===== Inserted 

# Attach observation package to DRN package
ghb.obs.initialize( # <===== Inserted 
    filename=gwf.name+".ghb.obs", # <===== Inserted 
    digits=10, # <===== Inserted 
    print_input=True, # <===== Inserted 
    continuous=obsDict # <===== Inserted 
) # <===== Inserted

#define buy package
buyModName = 'modelBuy'
Csalt = 35.
Cfresh = 0.
densesalt = 1025.
densefresh = 1000.
denseslp = (densesalt - densefresh) / (Csalt - Cfresh)

pd = [(0, denseslp, 0., buyModName, 'CONCENTRATION')]
buy = flopy.mf6.ModflowGwfbuy(gwf, denseref=1000., nrhospecies=1,packagedata=pd)

#ghb plot
ghb.plot(mflay=0, kper=0) # <===== Inserted
plt.show()

from copy import copy
#well bc

wellSpd = {} ## # <===== Inserted
wellSpd[0] = [] ## # <===== Inserted
wellSpd[1] = []

#regional flow
# from raster
# layCellTupleList, cellElevList = getLayCellElevTupleFromRaster(gwf,
#                                                                interIx,
#                                                                '../rst/modelDemMinus45.tif',
#                                                                '../shp/wellsStage1.shp') ## # <===== Inserted

# for index, layCellTuple in enumerate(layCellTupleList): ## Org
#     wellSpd[1].append([layCellTuple,-0.01,'wellsStage1']) # <===== Inserted

# from elevation
layCellTupleList = getLayCellElevTupleFromElev(gwf,
                                               interIx,
                                               -20,
                                               os.path.join(shpDir,'pumpingWells.shp'))
for layCellTuple in layCellTupleList:
    wellSpd[1].append([layCellTuple, -0.01, 'wellsStage1'])

wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=wellSpd, boundnames=True) ## <==== modified

# Observation package for Drain
obsDict = { # <===== Inserted 
    "{}.wel.obs.csv".format(modelName): [ # <===== Inserted 
        ("wellsStage1", "wel", "wellsStage1") # <===== Inserted 
    ] # <===== Inserted 
} # <===== Inserted 

# Attach observation package to DRN package
wel.obs.initialize( # <===== Inserted 
    filename=gwf.name+".wel.obs", # <===== Inserted 
    digits=10, # <===== Inserted 
    print_input=True, # <===== Inserted 
    continuous=obsDict # <===== Inserted 
) # <===== Inserted

refBounds = gpd.read_file(os.path.join(shpDir,'modelRef.shp')).total_bounds

#ghb plot
fig, ax = plt.subplots()

#
mV = flopy.plot.PlotMapView(gwf)
mV.plot_bc('WEL', kper=1, ax=ax, plotAll=True)
mV.plot_bc('GHB', kper=1, ax=ax)
mV.plot_grid(alpha=0.5, lw=0.2)
ax.set_xlim(refBounds[0]-1000,refBounds[2]+1000)
ax.set_ylim(refBounds[1]-1000,refBounds[3]+1000)

#open the river shapefile
rivers =gpd.read_file(os.path.join(shpDir, 'riverBasin.shp')) ## Org

list_rivers=[] ## Org
for i in range(rivers.shape[0]): ## Org
    list_rivers.append(rivers['geometry'].loc[i]) ## Org
    
riverMls = MultiLineString(lines=list_rivers) ## Org

#intersec rivers with our grid
riverCells=interIx.intersect(riverMls).cellids ## Org

#river package
riverSpd = {} ## Org
riverSpd[0] = [] ## Org
for cell in riverCells: ## Org
    riverSpd[0].append([(0,cell),mtop[cell],0.01]) ## Org

riv = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=riverSpd) ## Org
#river plot
riv.plot(mflay=0) ## Org

#create transport package
gwt = flopy.mf6.ModflowGwt(sim, modelname=buyModName)

#register solver for transport model
imsGwt = flopy.mf6.ModflowIms(sim,
                              pname='ims_gwt', 
                              #print_option='SUMMARY', ## Org 
                              outer_dvclose=2e-4, ## Org
                              inner_dvclose=3e-4, ## Org
                              linear_acceleration='BICGSTAB') ## Org
sim.register_ims_package(imsGwt,[gwt.name])

#define spatial discretization
gwtDisv = flopy.mf6.ModflowGwtdisv(gwt, nlay=disv.nlay.data,
                                   ncpl=disv.ncpl.data,
                                   nvert=disv.nvert.data,
                                   top=disv.top.data,
                                   botm=disv.botm.data,
                                   vertices=disv.vertices.array.tolist(),
                                   cell2d=disv.cell2d.array.tolist(),
                                  )

## Org

sim.register_ims_package(imsGwf,[modelName])

#define starting concentrations
strtConc = np.zeros((disv.nlay.data, disv.ncpl.data), dtype=np.float32)

ghbList = ghb.stress_period_data.array[0].tolist()

for ghbItem in ghbList:
    if ghbItem[4] == 'sea':
        strtConc[:,ghbItem[0][1]] = 35 #apply for all layers below the ghb
gwtIc = flopy.mf6.ModflowGwtic(gwt, strt=strtConc)

# create plot of initial concentratios
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(1, 1, 1, aspect = 'equal')
mapview = flopy.plot.PlotMapView(model=gwf,layer = 1)

plot_array = mapview.plot_array(strtConc,masked_values=[-1e+30], cmap=plt.cm.summer)
plt.colorbar(plot_array, shrink=0.75,orientation='horizontal', pad=0.08, aspect=50)
plt.show()

#define advection
adv = flopy.mf6.ModflowGwtadv(gwt, scheme='UPSTREAM')
#define dispersion
dsp = flopy.mf6.ModflowGwtdsp(gwt,alh=10,ath1=10)
#define mobile storage and transfer
porosity = 0.30
sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
#define sink and source package
sourcerecarray = ['GHB_0','AUX','CONCENTRATION']
ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

#define constant concentration package
cncSp = []
for row in ghb.stress_period_data.array[0]:
    if row['boundname'] == 'sea':
        cncSp.append([row[0],35])

cncSpd = {0:cncSp,1:cncSp}
cnc = flopy.mf6.ModflowGwtcnc(gwt,stress_period_data=cncSpd)
# cnc.plot(mflay=0, lw=0.1, figsize=(12,12))

#working with observation points 
obsList = []
nameList, obsLayCellList = getLayCellElevTupleFromObs(gwf, ## Org
                  interIx, ## Org
                  os.path.join(shpDir,'obsPoints.shp'), ## Org
                  'Name', ## Org
                  'Elev') ## Org

for obsName, obsLayCell in zip(nameList, obsLayCellList): ## Org
    obsList.append((obsName,'concentration',obsLayCell[0]+1,obsLayCell[1]+1)) ## Org


obs = flopy.mf6.ModflowUtlobs( ## Org
    gwt,
    filename=gwt.name+'.obs', ## Org
    digits=10, ## Org
    print_input=True, ## Org
    continuous={gwt.name+'.obs.csv': obsList} ## Org
)

#oc for flow 
head_filerecord = f"{gwf.name}.hds" ## Org
budget_filerecord = f"{gwf.name}.cbc" ## Org
oc = flopy.mf6.ModflowGwfoc(gwf, ## Org
                            head_filerecord=head_filerecord, ## Org
                            budget_filerecord = budget_filerecord, ## Org
                            saverecord=[("HEAD", "LAST"),("BUDGET","LAST")]) ## Org

#oc for transport
oc = flopy.mf6.ModflowGwtoc(gwt,
                            concentration_filerecord=buyModName+'.ucn',
                            saverecord=[('CONCENTRATION', 'ALL')])

#define model flow and transport exchange
name = 'modelExchange'
gwfgwt = flopy.mf6.ModflowGwfgwt(sim, exgtype='GWF6-GWT6',
                                 exgmnamea=gwf.name, exgmnameb=buyModName,
                                 filename='{}.gwfgwt'.format(name))

# Run the simulation
#sim.write_simulation() ## Org
success, buff = sim.run_simulation() ## Org