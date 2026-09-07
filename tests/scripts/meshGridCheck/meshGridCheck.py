from mf6Voronoi.meshProperties import meshShape
from datetime import datetime
import os, sys, json, flopy 
import matplotlib.pyplot as plt

#processing argument
if len(sys.argv) < 3:
    sys.exit(
        "Usage: python script.py <caseName> <runType> <caseType> [<fromCase>] [<fileType>]\n\n"
        "Arguments:\n"
        "  caseName   : Name of the test case\n"
        "  runType    : normalRun | parallelRun\n"
        "  caseType   : casesDask | casesNormal\n"
    )

caseName = sys.argv[1] 

if sys.argv[2] == 'normalRun':
    print('Computing with normal CPU')
elif sys.argv[2] == 'parallelRun':
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

# --------------------------
# Defining folders
# --------------------------

outputDataFolder = "/home/hatari/projects/mf6Voronoi/tests/output" 

# Get today's current date Format date: %d (day), %b (short month), %y (2-digit year)
todayStr = datetime.now().strftime("%d%b%y")

for meshName, meshDict in list(meshGenerationDict.items()):
    outDir = os.path.join(outputDataFolder, f"{caseName}_{todayStr}", sys.argv[2], meshName)
    outputShape = os.path.join(outDir, 'shp',f"{meshName}.shp")

    # open the mesh file
    mesh=meshShape(outputShape) ## Org
    # get the list of vertices and cell2d data
    gridprops=mesh.get_gridprops_disv() ## Org

    jsonDir = os.path.join(outDir, 'json')
    imgDir = os.path.join(outDir, 'jpg')
    os.makedirs(jsonDir, exist_ok=True)
    os.makedirs(imgDir, exist_ok=True)

    #export disv
    mesh.save_properties(os.path.join(jsonDir,'disvDict.json')) ## Org

    # Cargar el diccionario que acabas de exportar
    disv_json_path = os.path.join(jsonDir, 'disvDict.json')
    with open(disv_json_path, 'r') as f:
        disvData = json.load(f)

    # Crear modelo FloPy temporal
    sim = flopy.mf6.MFSimulation(sim_name="temp_sim", version="mf6")
    gwf = flopy.mf6.ModflowGwf(sim, modelname="temp_model")
    
    flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=1,
        ncpl=disvData['ncpl'],
        nvert=disvData['nvert'],
        vertices=disvData['vertices'],
        cell2d=disvData['cell2d']
    )

    # Mostrar/Guardar con PlotMapView
    fig, ax = plt.subplots(figsize=(10, 10))
    pmv = flopy.plot.PlotMapView(model=gwf, ax=ax)
    pmv.plot_grid(edgecolor="black", linewidth=0.4)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(f"Malla DISV: {meshName}")
    
    # Si estás corriendo en un servidor o script automatizado, guarda la figura:
    plot_path = os.path.join(imgDir, f"{meshName}_disv_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)