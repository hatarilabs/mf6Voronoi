from datetime import datetime
from pathlib import Path
import os
import sys
import json
import flopy
import platform
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

#Extract dem values for each centroid of the voronois
src = rasterio.open(os.path.join(rstDir,'asterDem18S.tif'))
elevation=[x for x in src.sample(centroids)]
nlay = 5

mtop=np.array([elev[0] for i,elev in enumerate(elevation)])
zbot=np.zeros((nlay,ncpl))


AcuifInf_Bottom = 2800
zbot[0,] = mtop - 30
zbot[1,] = AcuifInf_Bottom + (0.85 * (mtop - AcuifInf_Bottom))
zbot[2,] = AcuifInf_Bottom + (0.70 * (mtop - AcuifInf_Bottom))
zbot[3,] = AcuifInf_Bottom + (0.50 * (mtop - AcuifInf_Bottom))
zbot[4,] = AcuifInf_Bottom

#### Create simulation and model
# create simulation
simName = 'mf6Sim'
modelName = 'mf6Model'
sim = flopy.mf6.MFSimulation(sim_name=modelName, version='mf6', 
                             exe_name=mfBin, 
                             sim_ws=modelWs)
# create tdis package
tdis_rc = [(1000.0, 1, 1.0)]
tdis = flopy.mf6.ModflowTdis(sim, pname='tdis', time_units='DAYS', 
                             perioddata=tdis_rc)
# create gwf model
gwf = flopy.mf6.ModflowGwf(sim, modelname=modelName, save_flows=True)
# create iterative model solution and register the gwf model with it
ims = flopy.mf6.ModflowIms(sim,
                           complexity='COMPLEX',
                           outer_maximum=100,
                           inner_maximum=100, 
                           linear_acceleration='BICGSTAB',
                          )
sim.register_ims_package(ims,[modelName])
# disv
disv = flopy.mf6.ModflowGwfdisv(gwf, nlay=nlay, ncpl=ncpl, 
                                top=mtop, botm=zbot, 
                                nvert=nvert, vertices=vertices, 
                                cell2d=cell2d)
# initial conditions
ic = flopy.mf6.ModflowGwfic(gwf, strt=np.stack([mtop for i in range(nlay)]))
Kx =[4E-4,5E-6,1E-6,9E-7,5E-7]
icelltype = [1,1,0,0,0]

# node property flow
npf = flopy.mf6.ModflowGwfnpf(gwf, xt3doptions=[('xt3d')],
                              save_specific_discharge=True,
                              icelltype=icelltype, 
                              k=Kx)
# define storage and transient stress periods
sto = flopy.mf6.ModflowGwfsto(gwf,
                              iconvert=1,
                              steady_state={
                                0:True,
                              }
                              )

#### Working with rechage, evapotranspiration
rchr = 0.15/365/86400
rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rchr)
evtr = 1.2/365/86400
evt = flopy.mf6.ModflowGwfevta(gwf,ievt=1,surface=mtop,rate=evtr,depth=1.0)

#### Definition of the intersect object
#### For the manipulation of spatial data to determine hydraulic parameters or boundary conditions

# Define intersection object
interIx = flopy.utils.gridintersect.GridIntersect(gwf.modelgrid)
from shapely.geometry import MultiLineString

#open the river shapefile
rivers =gpd.read_file(os.path.join(shpDir,'river_basin.shp'))
list_rivers=[]
for i in range(rivers.shape[0]):
    list_rivers.append(rivers['geometry'].loc[i])
    
riverMls = MultiLineString(lines=list_rivers)

#intersec rivers with our grid
riverCells=interIx.intersect(riverMls).cellids

#river package
riverSpd = {}
riverSpd[0] = []
for cell in riverCells:
    riverSpd[0].append([(0,cell),mtop[cell],0.01]) 
riv = flopy.mf6.ModflowGwfdrn(gwf, stress_period_data=riverSpd)
#river plot
riv.plot(mflay=0)
plt.show()

#oc
head_filerecord = f"{gwf.name}.hds" ## Org
budget_filerecord = f"{gwf.name}.cbc" ## Org
oc = flopy.mf6.ModflowGwfoc(gwf, ## Org
                            head_filerecord=head_filerecord, ## Org
                            budget_filerecord = budget_filerecord, ## Org
                            saverecord=[("HEAD", "LAST"),("BUDGET","LAST")]) ## Org

# Run the simulation
sim.write_simulation()
success, buff = sim.run_simulation()