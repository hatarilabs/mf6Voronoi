from mf6Voronoi.geoVoronoi import createVoronoi
from mf6Voronoi.utils import getVoronoiAsShp
import matplotlib.pyplot as plt
from datetime import datetime
import geopandas as gpd
from pathlib import Path
import os, json, sys

# --------------------------
# Processing CLI arguments
# --------------------------
#processing argument
if len(sys.argv) < 4:
    sys.exit(
        "Usage: python script.py <caseName> <runType> <caseType> [<fromCase>] [<fileType>]\n\n"
        "Arguments:\n"
        "  caseName   : Name of the test case\n"
        "  runType    : normalRun | parallelRun\n"
        "  caseType   : casesDask | casesNormal\n"
        "  fromCase   : Optional start index (default: 0)\n"
        "  fileType   : shp | geojson (default: Shapefile)" \
        "  debug      : True | False (default: False)"
    )

caseName = sys.argv[1]

if sys.argv[2] == 'normalRun':
    useDask = False
    nproc = 1
    print('Computing with normal CPU')
elif sys.argv[2] == 'parallelRun':
    useDask = True
    nproc = 4      
    print('Computing with parallelized computing')
else:
    sys.exit(
        f"Invalid execution mode: '{sys.argv[2]}'. Expected 'normalRun' or 'parallelRun'."
    )

# Case file setup
if sys.argv[3] == "casesDask":
    json_path = "../../json/meshCasesDask.json"
    print("Working with Dask cases")
elif sys.argv[3] == "casesNormal":
    json_path = "../../json/meshCasesNormal.json"
    print("Working with normal cases")
else:
    sys.exit("The given mesh case is wrong capullo!")

with open(json_path) as jsonFile:
    meshGenerationDict = json.load(jsonFile)

# Processing fromCase 
try :
    fromCase = int(sys.argv[4])
except IndexError:
    fromCase = 0

# Processing fileType 
try :
    fileType = sys.argv[5]
except IndexError:
    fileType = 'shp'

# Processing debugModel 
try :
    debugMode = sys.argv[6]
    if debugMode == 'True':
        debug = True
    else:
        debug = False
except IndexError:
    debug = False

# --------------------------
# Defining folders
# --------------------------
testDataFolder = "/home/hatari/projects/mf6Voronoi/tests/data"
outputDataFolder = "/home/hatari/projects/mf6Voronoi/tests/output"

# Get today's current date Format date: %d (day), %b (short month), %y (2-digit year)
todayStr = datetime.now().strftime("%d%b%y")

# Process mesh cases
for meshName, meshDict in list(meshGenerationDict.items())[fromCase:]:
    if fileType == 'shp':
        datasetPath = os.path.join(testDataFolder, meshName, "shp")
    elif fileType == 'geojson':
        datasetPath = os.path.join(testDataFolder, meshName, "geojson")
    else:
        print('El tipo de archivo espacial no existe capullo')
    verifDir = os.path.join(outputDataFolder, f"{caseName}_{todayStr}", sys.argv[2], meshName,'shp')

    outputShape = os.path.join(verifDir, f"{meshName}.shp")
    os.makedirs(os.path.dirname(outputShape), exist_ok=True)

    if 'overlapping' in meshDict.keys():
        overlapping = meshDict['overlapping']
    else:
        overlapping = True


    if useDask:

        #Create mesh object specifying the coarse mesh and the multiplier
        vorMesh = createVoronoi(meshName=meshName,
                                maxRef = meshDict["maxRef"], 
                                multiplier=meshDict["multiplier"],
                                use_dask=useDask, 
                                nproc=nproc,
                                overlapping=overlapping)

        #Open limit layers and refinement definition layers
        vorMesh.addLimit(meshDict["limitLayer"]["limitName"],
                         os.path.join(datasetPath, meshDict["limitLayer"]["limitShape"]+'.'+fileType))

        for layerList in  meshDict["layerLayer"]:
            vorMesh.addLayer(layerList[0],
                             os.path.join(datasetPath,layerList[1]+'.'+fileType),
                             layerList[2])

        vorMesh.generateOrgDistVertices(debug=debug, out_dir=verifDir)
        vorMesh.createPointCloud(debug=debug, out_dir=verifDir)
        vorMesh.generateVoronoi(shapePath=outputShape)

    else:
        #Create mesh object specifying the coarse mesh and the multiplier
        vorMesh = createVoronoi(meshName=meshName,
                                maxRef = meshDict["maxRef"], 
                                multiplier=meshDict["multiplier"],
                                overlapping=overlapping)

        #Open limit layers and refinement definition layers
        vorMesh.addLimit(meshDict["limitLayer"]["limitName"], 
                         os.path.join(datasetPath,meshDict["limitLayer"]["limitShape"]+'.'+fileType))

        for layerList in  meshDict["layerLayer"]:
            vorMesh.addLayer(layerList[0],
                             os.path.join(datasetPath,layerList[1]+'.'+fileType),
                             layerList[2])
        if debug:
            vorMesh.generateOrgDistVertices(debug=True, out_dir=verifDir)
            vorMesh.createPointCloud(debug=True, out_dir=verifDir)
            vorMesh.generateVoronoi(shapePath=outputShape)
        else:
            vorMesh.generateOrgDistVertices()
            vorMesh.createPointCloud()
            vorMesh.generateVoronoi()
            getVoronoiAsShp(vorMesh.modelDis, shapePath=outputShape)
    