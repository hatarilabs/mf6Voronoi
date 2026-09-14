from mf6Voronoi.geoVoronoi import createVoronoi
from mf6Voronoi.utils import getVoronoiAsShp
from datetime import datetime
import os, json, sys, time

# --------------------------
# Processing CLI arguments
# --------------------------
if len(sys.argv) < 3:
    sys.exit(
        "Usage: python daskBenchmark.py <caseName> <caseType>\n\n"
        "Arguments:\n"
        "  caseName   : Name of the test case\n"
        "  caseType   : casesSingle | casesSingleDask"
    )

caseName = sys.argv[1]
caseType = sys.argv[2]

# Case file setup
if caseType == "casesSingle":
    json_path = "../../json/meshCasesSingle.json"
    print("Working with single normal cases")
elif caseType == "casesSingleDask":
    json_path = "../../json/meshCasesSingleDask.json"
    print("Working with single Dask cases")
else:
    sys.exit(f"Invalid caseType: '{caseType}'. Expected 'casesSingle'.")

with open(json_path) as jsonFile:
    meshGenerationDict = json.load(jsonFile)

# --------------------------
# Defining folders and Benchmark File
# --------------------------
testDataFolder = "/home/hatari/projects/mf6Voronoi/tests/data"
outputDataFolder = "/home/hatari/projects/mf6Voronoi/tests/output"

todayStr = datetime.now().strftime("%d%b%y")
benchmarkDir = os.path.join(outputDataFolder, f"{caseName}_{todayStr}")
os.makedirs(benchmarkDir, exist_ok=True)

log_file_path = os.path.join(benchmarkDir, "benchmark_results.txt")

# Write CSV header if file does not exist yet
if not os.path.exists(log_file_path):
    with open(log_file_path, "w") as log_file:
        log_file.write("MeshName,RunType,NProc,TimePointGen,TimeVoronoiGen,TimeShapefileGen,TotalTime\n")

# Process mesh cases
for meshName, meshDict in meshGenerationDict.items():
    datasetPath = os.path.join(testDataFolder, meshName, "shp")
    overlapping = meshDict.get('overlapping', True)

    # Sequence of runs: normalRun (n=1) then parallelRun from n=2 to n=12
    runs = [('normalRun', False, 1)] + [('parallelRun', True, n) for n in range(2, 13)]

    for runType, useDask, nproc in runs:
        print(f"\n--- Running {meshName} | Mode: {runType} | Cores: {nproc} ---")
        
        verifDir = os.path.join(benchmarkDir, runType, f"proc_{nproc}", meshName, 'shp')
        outputShape = os.path.join(verifDir, f"{meshName}.shp")
        os.makedirs(os.path.dirname(outputShape), exist_ok=True)

        if useDask:
            vorMesh = createVoronoi(
                meshName=meshName,
                maxRef=meshDict["maxRef"], 
                multiplier=meshDict["multiplier"],
                use_dask=useDask, 
                nproc=nproc,
                overlapping=overlapping
            )

            vorMesh.addLimit(
                meshDict["limitLayer"]["limitName"],
                os.path.join(datasetPath, meshDict["limitLayer"]["limitShape"] + '.shp')
            )
            for layerList in meshDict["layerLayer"]:
                vorMesh.addLayer(
                    layerList[0],
                    os.path.join(datasetPath, layerList[1] + '.shp'),
                    layerList[2]
                )

            # 1. Time required for point generation
            t0 = time.perf_counter()
            vorMesh.generateOrgDistVertices()
            vorMesh.createPointCloud()
            t_point_gen = time.perf_counter() - t0

            # 2. Time required for voronoi generation (Dask handles direct export inside generateVoronoi)
            t1 = time.perf_counter()
            vorMesh.generateVoronoi(shapePath=outputShape)
            t_vor_gen = time.perf_counter() - t1

            # 3. Time required for voronoi shapefile (Included in voronoi gen step for Dask)
            t_shp_gen = 0.0

        else:
            vorMesh = createVoronoi(
                meshName=meshName,
                maxRef=meshDict["maxRef"], 
                multiplier=meshDict["multiplier"],
                overlapping=overlapping
            )

            vorMesh.addLimit(
                meshDict["limitLayer"]["limitName"], 
                os.path.join(datasetPath, meshDict["limitLayer"]["limitShape"] + '.shp')
            )
            for layerList in meshDict["layerLayer"]:
                vorMesh.addLayer(
                    layerList[0],
                    os.path.join(datasetPath, layerList[1] + '.shp'),
                    layerList[2]
                )

            # 1. Time required for point generation
            t0 = time.perf_counter()
            vorMesh.generateOrgDistVertices()
            vorMesh.createPointCloud()
            t_point_gen = time.perf_counter() - t0

            # 2. Time required for voronoi generation
            t1 = time.perf_counter()
            vorMesh.generateVoronoi()
            t_vor_gen = time.perf_counter() - t1

            # 3. Time required for voronoi shapefile
            t2 = time.perf_counter()
            getVoronoiAsShp(vorMesh.modelDis, shapePath=outputShape)
            t_shp_gen = time.perf_counter() - t2

        t_total = t_point_gen + t_vor_gen + t_shp_gen

        # Write results row to TXT file
        log_line = f"{meshName},{runType},{nproc},{t_point_gen:.4f},{t_vor_gen:.4f},{t_shp_gen:.4f},{t_total:.4f}\n"
        with open(log_file_path, "a") as log_file:
            log_file.write(log_line)

        print(f"Completed in {t_total:.2f}s | PointGen: {t_point_gen:.2f}s | VoronoiGen: {t_vor_gen:.2f}s | ShpGen: {t_shp_gen:.2f}s")