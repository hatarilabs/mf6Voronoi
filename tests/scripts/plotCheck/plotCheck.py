import os
import sys
import json
from datetime import datetime
import matplotlib.pyplot as plt
import geopandas as gpd
import flopy

from mf6Voronoi.tools.graphs2d import FlowVectorGenerator, crossSectionFlowVectorGenerator

# Processing arguments
if len(sys.argv) < 3:
    sys.exit(
        "Usage: python script.py <caseName> <runType>\n\n"
        "Arguments:\n"
        "  caseName   : Name of the test case\n"
        "  runType    : normalRun | parallelRun\n"
    )

caseName = sys.argv[1] 
runType = sys.argv[2]

if runType == 'normalRun':
    print('Computing with normal CPU')
elif runType == 'parallelRun':
    print('Computing with parallelized computing')
else:
    sys.exit(
        f"Invalid execution mode: '{runType}'. Expected 'normalRun' or 'parallelRun'."
    )

# --------------------------
# 1. Configuración de Argumentos y JSONs
# --------------------------
jsonPath = '../../json/modelCases.json'
csJsonPath = '../../json/crossSections.json'

with open(jsonPath, "r") as f:
    configDict = json.load(f)

with open(csJsonPath, "r") as f:
    csDict = json.load(f)["cross_sections"]

caseGenDataFolder = configDict["base_dirs"]["case_gen_data"]
meshGenDataFolder = configDict["base_dirs"]["mesh_gen_data"]
todayStr = datetime.now().strftime("%d%b%y")

# --------------------------
# 2. Procesamiento Automatizado de Todos los Casos
# --------------------------
for modelNumber, meshName in configDict["mesh_name"].items():
    print(f"\n==========================================")
    print(f"Procesando Caso [{modelNumber}]: {meshName}")
    print(f"==========================================")

    # Construcción de rutas según la estructura del proyecto
    meshCaseDir = os.path.join(meshGenDataFolder, f"{caseName}_{todayStr}", runType, meshName)
    modelWs = os.path.join(meshCaseDir, 'model')
    shpDir = os.path.join(caseGenDataFolder, meshName, 'shp')
    imgDir = os.path.join(meshCaseDir, 'img')

    os.makedirs(imgDir, exist_ok=True)

    if meshName not in csDict:
        print(f"Omitiendo: '{meshName}' no tiene configuración en crossSections.json.")
        continue

    shpName = csDict[meshName]["shpName"]
    csShpPath = os.path.join(shpDir, shpName)

    if not os.path.exists(csShpPath):
        print(f"Error: No existe el shapefile en {csShpPath}")
        continue

    csGdf = gpd.read_file(csShpPath)
    csGeometry = csGdf.geometry.iloc[0]

    simName = 'mf6Sim'
    modelName = 'mf6Model'
    mfBin = os.path.join(caseGenDataFolder, 'modflowBin', 'mf6')

    hdsFile = os.path.join(modelWs, f"{modelName}.hds")
    cbcFile = os.path.join(modelWs, f"{modelName}.cbc")

    if not (os.path.exists(hdsFile) and os.path.exists(cbcFile)):
        print(f"Error: No se encontraron los archivos de salida binarios en {modelWs}")
        continue

    sim = flopy.mf6.MFSimulation.load(
        sim_name=simName,
        sim_ws=modelWs,
        exe_name=mfBin
    )
    gwf = sim.get_model(modelName)

    # Obtenemos la última capa disponible en el modelo
    lastLayer = gwf.modelgrid.nlay - 1
    print(f" Graficando última capa: Capa {lastLayer + 1} (índice {lastLayer})")

    # --------------------------
    # 3. Generación y Guardado: Planta (Última Capa)
    # --------------------------
    print(" Generando figura de flujo en planta...")
    fig = FlowVectorGenerator(
        gwf=gwf,
        layer=lastLayer,
        plotGrid=True,
        plotContour=True,
        scale=20
    )
    
    if fig is not None:
        planImgPath = os.path.join(imgDir, "flow_vectors_plan.png")
        fig.savefig(planImgPath, dpi=300, bbox_inches='tight')
        plt.close(fig)

    # --------------------------
    # 4. Generación y Guardado: Sección Transversal Completa
    # --------------------------
    print(" Generando figura completa en sección transversal (cargas, isolíneas y flechas)...")
    figCs, ax = plt.subplots(figsize=(12, 6))

    # Formato de definición de línea para PlotCrossSection
    line_dict = {'line': csGeometry}

    # Generación de sección transversal con mapa de color, isolíneas y vectores
    crossSectionFlowVectorGenerator(
        gwf,
        cbc_file=cbcFile,
        head_file=hdsFile,
        line=line_dict,
        ax=ax,
        kstpkper=(0, 0),
        istep=1,
        jstep=1,
        plotArray=True,       # Renderiza mapa de color de cargas hidráulicas
        plotContour=True,     # Dibuja isolíneas equipotenciales
        contourLevels=10,     # Número de niveles de isolíneas
        normalize=False,
        color='black',
        alpha=0.8
    )

    #ax.grid(True, lw=0.3, ls='--', alpha=0.5)
    ax.set_title(f"Sección Transversal de Cargas e Isolíneas con Flujo - {meshName}")
    ax.set_xlabel("Distancia a lo largo del corte (m)")
    ax.set_ylabel("Elevación (m.s.n.m.)")

    csImgPath = os.path.join(imgDir, "flow_vectors_cross_section.png")
    figCs.savefig(csImgPath, dpi=300, bbox_inches='tight')
    plt.close(figCs)

    print(f" Figuras guardadas correctamente en:\n   -> {imgDir}")

print("\nProcesamiento de todos los casos completado.")