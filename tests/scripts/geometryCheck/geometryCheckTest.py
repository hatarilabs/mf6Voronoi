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
        "Usage: python script.py <caseName> <runType: normalRun|parallelRun> <caseType: casesDask|casesNormal> <fromCase: 0>"
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
    json_path = "meshCasesDask.json"
    print("Working with Dask cases")
elif sys.argv[3] == "casesNormal":
    json_path = "meshCasesNormal.json"
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

# --------------------------
# Defining folders
# --------------------------
testDataFolder = "/home/hatari/projects/mf6Voronoi/tests/data"
outputDataFolder = "/home/hatari/projects/mf6Voronoi/tests/output"

# Get today's current date Format date: %d (day), %b (short month), %y (2-digit year)
todayStr = datetime.now().strftime("%d%b%y")

# Process mesh cases
for meshName, meshDict in list(meshGenerationDict.items())[fromCase:]:
    datasetPath = os.path.join(testDataFolder, meshName, "shp")
    verifDir = os.path.join(outputDataFolder, f"{caseName}_{todayStr}", sys.argv[2], meshName, )
    outputShape = os.path.join(verifDir, f"{meshName}.shp")

    os.makedirs(os.path.dirname(outputShape), exist_ok=True)

    if useDask:

        #Create mesh object specifying the coarse mesh and the multiplier
        vorMesh = createVoronoi(meshName=meshName,
                                maxRef = meshDict["maxRef"], 
                                multiplier=meshDict["multiplier"],
                                use_dask=useDask, 
                                nproc=nproc)

        #Open limit layers and refinement definition layers
        vorMesh.addLimit(meshDict["limitLayer"]["limitName"],
                         os.path.join(datasetPath, meshDict["limitLayer"]["limitShape"]))

        for layerList in  meshDict["layerLayer"]:
            vorMesh.addLayer(layerList[0],
                             os.path.join(datasetPath,layerList[1]),
                             layerList[2])

        vorMesh.generateOrgDistVertices(debug=True, out_dir=verifDir)
        vorMesh.createPointCloud(debug=True, out_dir=verifDir)
        vorMesh.generateVoronoi(shapePath=outputShape)

    else:
        #Create mesh object specifying the coarse mesh and the multiplier
        vorMesh = createVoronoi(meshName=meshName,
                                maxRef = meshDict["maxRef"], 
                                multiplier=meshDict["multiplier"])

        #Open limit layers and refinement definition layers
        vorMesh.addLimit(meshDict["limitLayer"]["limitName"], 
                         os.path.join(datasetPath,meshDict["limitLayer"]["limitShape"]))

        for layerList in  meshDict["layerLayer"]:
            vorMesh.addLayer(layerList[0],
                             os.path.join(datasetPath,layerList[1]),
                             layerList[2])

        #Generate point pair array
        vorMesh.generateOrgDistVertices()

        #Generate the point cloud 
        vorMesh.createPointCloud()

        #generate voronoi and export directly the shapefile
        vorMesh.generateVoronoi()

        getVoronoiAsShp(vorMesh.modelDis, shapePath=outputShape)
    