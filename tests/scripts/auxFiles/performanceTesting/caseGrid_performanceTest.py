import os
import json
import sys
import time
from datetime import datetime
from mf6Voronoi.geoVoronoi import createVoronoi

# --------------------------
# Processing CLI arguments
# --------------------------
if len(sys.argv) < 2:
    sys.exit(
        "Usage: python cellScalingBenchmark.py <caseName> [jsonFileName]\n\n"
        "Arguments:\n"
        "  caseName     : Name of the test case\n"
        "  jsonFileName : Optional JSON file name (default: meshCasesSingleDask.json)"
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

log_file_path = os.path.join(benchmarkDir, "cell_benchmark_results.txt")

# Write header if file does not exist
if not os.path.exists(log_file_path):
    with open(log_file_path, "w") as log_file:
        log_file.write("MeshName,Iteration,TotalPoints,TimePointGen,TimeVoronoiGen,TimeShapefileGen,TotalTime\n")

# Process mesh cases
for meshName, meshDict in meshGenerationDict.items():
    datasetPath = os.path.join(testDataFolder, meshName, "shp")
    overlapping = meshDict.get('overlapping', True)
    
    # Store original refinement values
    base_maxRef = meshDict["maxRef"]
    base_layer_refs = [layer[2] for layer in meshDict["layerLayer"]]

    # Run original case (Iteration 0) + 10 iterations decreasing refinement by 10% each time (factor = 0.9^k)
    for k in range(11):
        factor = (0.9) ** k
        current_maxRef = base_maxRef * factor

        print(f"\n--- Running {meshName} | Iteration {k} (Refinement Factor: {factor:.4f}) ---")
        
        verifDir = os.path.join(benchmarkDir, f"iter_{k}", meshName, 'shp')
        outputShape = os.path.join(verifDir, f"{meshName}.shp")
        os.makedirs(os.path.dirname(outputShape), exist_ok=True)

        # Initialize Voronoi instance always with Dask
        vorMesh = createVoronoi(
            meshName=meshName,
            maxRef=current_maxRef, 
            multiplier=meshDict["multiplier"],
            use_dask=True,
            nproc=6,
            overlapping=overlapping
        )

        vorMesh.addLimit(
            meshDict["limitLayer"]["limitName"],
            os.path.join(datasetPath, meshDict["limitLayer"]["limitShape"] + '.shp')
        )
        for idx, layerList in enumerate(meshDict["layerLayer"]):
            current_layerRef = base_layer_refs[idx] * factor
            vorMesh.addLayer(
                layerList[0],
                os.path.join(datasetPath, layerList[1] + '.shp'),
                current_layerRef
            )

        # 1. Time required for point generation
        t0 = time.perf_counter()
        vorMesh.generateOrgDistVertices()
        vorMesh.createPointCloud()
        t_point_gen = time.perf_counter() - t0

        # Get total points generated inside the limit
        total_points = len(vorMesh.modelDis['vertexTotal'])

        # 2. Time required for voronoi generation
        t1 = time.perf_counter()
        vorMesh.generateVoronoi(shapePath=outputShape)
        t_vor_gen = time.perf_counter() - t1

        # 3. Time required for voronoi shapefile (0.0 when using Dask)
        t_shp_gen = 0.0

        t_total = t_point_gen + t_vor_gen + t_shp_gen

        # Log row to TXT file
        log_line = f"{meshName},{k},{total_points},{t_point_gen:.4f},{t_vor_gen:.4f},{t_shp_gen:.4f},{t_total:.4f}\n"
        with open(log_file_path, "a") as log_file:
            log_file.write(log_line)

        print(f"Points: {total_points:,} | Total Time: {t_total:.2f}s | PointGen: {t_point_gen:.2f}s | VoronoiGen: {t_vor_gen:.2f}s | ShpGen: {t_shp_gen:.2f}s")