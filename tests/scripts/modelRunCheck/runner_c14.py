from datetime import datetime
from pathlib import Path
import os
import sys
import json
import flopy
import rasterio
import platform
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
mf_ext = ".exe" if platform.system() == "Windows" else ""
mfBin = os.path.join(caseGenDataFolder,'modflowBin',f'mf6{mf_ext}')

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

nlay = 3 ## Org

mtop=np.array([50 for i in range(ncpl)]) ## <=== updated
zbot=np.zeros((nlay,ncpl)) ## <=== updated

zbot[0,] = [30 for i in range(ncpl)] ## <=== updated
zbot[1,] = [10 for i in range(ncpl)] ## <=== updated
zbot[2,] = [-10 for i in range(ncpl)] ## <=== updated

# create simulation
simName = 'mf6Sim' ## Org
modelName = 'mf6Model' ## Org
sim = flopy.mf6.MFSimulation(sim_name=modelName, version='mf6', ## Org
                             exe_name=mfBin, ## Org
                             sim_ws=modelWs) ## Org

# create tdis package
tdis_rc = [(1.0, 1, 1.0)] + [(86400*10, 5, 1.5) for level in range(2)] ## Org
print(tdis_rc[:3]) ## Org

tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='SECONDS', ## Org
                             perioddata=tdis_rc, ## Org
                            nper=3) ## Org

# create gwf model
gwf = flopy.mf6.ModflowGwf(sim, ## Org
                           modelname=modelName, ## Org
                           save_flows=True, ## Org
                           newtonoptions="NEWTON UNDER_RELAXATION") ## Org

# create iterative model solution and register the gwf model with it
ims = flopy.mf6.ModflowIms(sim, ## Org
                           complexity='COMPLEX', ## Org
                           outer_maximum=150, ## Org
                           inner_maximum=50, ## Org
                           outer_dvclose=0.1, ## Org
                           inner_dvclose=0.0001, ## Org
                           backtracking_number=20, ## Org
                           linear_acceleration='BICGSTAB') ## Org
sim.register_ims_package(ims,[modelName]) ## Org

# disv
disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, ## Org
                                top=mtop, botm=zbot, ## Org
                                nvert=nvert, vertices=vertices, ## Org
                                cell2d=cell2d) ## Org

disv.top.plot(figsize=(12,8), alpha=0.8) ## Org
plt.show()

crossSection = gpd.read_file(os.path.join(shpDir,'crossSectionTotal.shp')) ## Org
sectionLine =list(crossSection.iloc[0].geometry.coords) ## Org

fig, ax = plt.subplots(figsize=(12,8)) ## Org
modelxsect = flopy.plot.PlotCrossSection(model=gwf, line={'Line': sectionLine}) ## Org
linecollection = modelxsect.plot_grid(lw=0.5) ## Org
ax.grid() ## Org

plt.show()

# initial conditions

ic = flopy.mf6.ModflowGwfic(gwf, strt=np.stack([40 for i in range(nlay)])) ## Org
#headsInitial = np.load('npy/headCalibInitial.npy')
#ic = flopy.mf6.ModflowGwfic(gwf, strt=headsInitial)

Kx =[3E-5 for x in range(3)] ## <=== updated
icelltype = [1, 0, 0] ## <=== updated

# node property flow
npf = flopy.mf6.ModflowGwfnpf(gwf, ## Org
                              save_specific_discharge=True, ## Org
                              icelltype=icelltype, ## Org
                              k=Kx, ## Org
                              k33=Kx) ## Org

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

#river package
layCellTupleList = getLayCellElevTupleFromElev(gwf,interIx,39.5,os.path.join(shpDir,'river.shp')) ## <=== updated
riverSpd = {} ## Org
riverSpd[0] = [] ## Org
for index, layCellTuple in enumerate(layCellTupleList): ## Org
    riverSpd[0].append([layCellTuple,39.5,0.01,38]) ## Org
riverSpd[0][:5]

riv = flopy.mf6.ModflowGwfriv(gwf, stress_period_data=riverSpd) ## Org

#river plot
riv.plot(mflay=0, kper=1) ## Org
plt.show()

crossSection = gpd.read_file(os.path.join(shpDir,'crossSectionTotal.shp')) ## Org
sectionLine =list(crossSection.iloc[0].geometry.coords) ## Org

fig, ax = plt.subplots(figsize=(12,8)) ## Org
xsect = flopy.plot.PlotCrossSection(model=gwf, line={'Line': sectionLine}) ## Org
lc = xsect.plot_grid(lw=0.5) ## Org
xsect.plot_bc('RIV',kper=2) ## Org
ax.grid() ## Org
plt.show()

#regional flow package
layCellTupleList = getLayCellElevTupleFromElev(gwf,interIx,40,os.path.join(shpDir,'regionalFlow.shp')) ## <=== updated
ghbSpd = {} ## Org
ghbSpd[0] = [] ## Org
for index, layCellTuple in enumerate(layCellTupleList): ## <=== updated
    ghbSpd[0].append([layCellTuple,40,0.01]) ## <=== updated

ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghbSpd)
#regional flow plot
ghb.plot(mflay=0, kper=0) ## <===== modified
plt.show()

#well package
layCellTupleList = getLayCellElevTupleFromElev(gwf,interIx,20,os.path.join(shpDir,'wells.shp')) ## <=== updated
welSpd = {} ## Org
welSpd[0] = [] ## Org
for index, layCellTuple in enumerate(layCellTupleList): ## <=== updated
    welSpd[0].append([layCellTuple,-0.001]) ## <=== updated


wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=welSpd)
#regional flow plot
#wel.plot(mflay=1, kper=0, ec='crimson') ## <===== modified

#well package
layCellTupleListUzf = getLayCellElevTupleFromElev(gwf,interIx,40,os.path.join(shpDir,'piezometer.shp')) ## <=== updated
print(layCellTupleListUzf)
piezoCell = layCellTupleListUzf[0][1]

layCellTupleList = getLayCellElevTupleFromElev(gwf,interIx,50,os.path.join(shpDir,'infiltrationPond.shp')) ## <=== updated

packageData = []
for index, layCellTuple in enumerate(layCellTupleList): ## <=== updated
    if layCellTuple[1] == piezoCell:
        packageData.append([index, layCellTuple, 1, 0, 0.1, 3e-5, 0.1, 0.35, 0.15, 3.5,'surfRate'+str(piezoCell)]) ## <=== updated
    else:
        packageData.append([index, layCellTuple, 1, 0, 0.1, 3e-5, 0.1, 0.35, 0.15, 3.5,'surfRate']) ## <=== updated

periodData = {}
periodData[0] = []
periodData[1] = []
periodData[2] = []
for index, layCellTuple in enumerate(layCellTupleList): ## <=== updated
    periodData[0].append([index, 1.11e-06, 5e-08, 1.5, 0.1, 0, 0, 0])
    periodData[1].append([index, 1.67e-06, 5e-08, 1.5, 0.1, 0, 0, 0])
    periodData[2].append([index, 2.68e-06, 5e-08, 1.5, 0.1, 0, 0, 0])

# Create UZF package and parameters
uzf = flopy.mf6.ModflowGwfuzf(gwf,
                              ntrailwaves=7,
                              nwavesets=40,
                              simulate_et=True,
                              simulate_gwseep=False,
                              packagedata=packageData,
                              perioddata=periodData,
                              wc_filerecord='uzf.out',
                              boundnames=True)

# Observation package for Drain
uzfDict = { # <===== Inserted 
    "{}.uzf.obs.csv".format(modelName): [ # <===== Inserted 
        ("wc_0.2", "water-content", "surfRate"+str(piezoCell), 0.2), # <===== Inserted 
        ("wc_0.5", "water-content", "surfRate"+str(piezoCell), 0.5), # <===== Inserted 
        ("wc_1", "water-content", "surfRate"+str(piezoCell), 1), # <===== Inserted 
        ("wc_1.5", "water-content", "surfRate"+str(piezoCell), 1.5), # <===== Inserted 
        ("wc_2", "water-content", "surfRate"+str(piezoCell), 2), # <===== Inserted 
        ("wc_3", "water-content", "surfRate"+str(piezoCell), 3), # <===== Inserted 
        ("wc_4", "water-content", "surfRate"+str(piezoCell), 4), # <===== Inserted 
        ("wc_5", "water-content", "surfRate"+str(piezoCell), 5), # <===== Inserted 
        ("wc_6", "water-content", "surfRate"+str(piezoCell), 6), # <===== Inserted 
        ("wc_7", "water-content", "surfRate"+str(piezoCell), 7), # <===== Inserted 
        ("wc_8", "water-content", "surfRate"+str(piezoCell), 8), # <===== Inserted 
        ("wc_9", "water-content", "surfRate"+str(piezoCell), 9), # <===== Inserted 
    ] # <===== Inserted 
} # <===== Inserted 

# Attach observation package to DRN package
uzf.obs.initialize( # <===== Inserted 
    filename=gwf.name+".uzf.obs", # <===== Inserted 
    digits=10, # <===== Inserted 
    print_input=True, # <===== Inserted 
    continuous=uzfDict # <===== Inserted 
) # <===== Inserted 

#oc
head_filerecord = f"{gwf.name}.hds" ## Org
budget_filerecord = f"{gwf.name}.cbc" ## Org
oc = flopy.mf6.ModflowGwfoc(gwf, ## Org
                            head_filerecord=head_filerecord, ## Org
                            budget_filerecord = budget_filerecord, ## Org
                            saverecord=[("HEAD", "LAST"),("BUDGET","LAST")]) ## Org

# Run the simulation
sim.write_simulation() ## Org
success, buff = sim.run_simulation() ## Org