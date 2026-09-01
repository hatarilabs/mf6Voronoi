import flopy ## Org
from mf6Voronoi.tools.vtkGen import Mf6VtkGenerator ## Org
from mf6Voronoi.utils import initiateOutputFolder ## Org
import matplotlib.pyplot as plt
import geopandas as gpd
import os, json

with open('vtkGenerationData.json') as jsonFile:
    vtkGenerationDict = json.load(jsonFile)
    
for modelName, modelDict in vtkGenerationDict.items():

    print("\n Working for model: %s"%modelName)

    # load simulation
    simName = 'mf6Sim' ## Org

    modelWs = os.path.abspath(os.path.join('../../../examples/mf6models',modelName))

    sim = flopy.mf6.MFSimulation.load(sim_name=simName, version='mf6', ## Org
                                exe_name='../../../../bin/mf6', ## Org
                                sim_ws=modelWs) ## Org
    
    vtkDir = '../../checkFiles/postProcessing/vtkGeneration/'+modelName ## Org
    initiateOutputFolder(vtkDir) ## Org

    mf6Vtk = Mf6VtkGenerator(sim, vtkDir) ## Org

    mf6Vtk.loadModel(modelDict['gwfname']) ## Org

    #show output data
    headObj = mf6Vtk.gwf.output.head() ## Org
    headObj.get_kstpkper() ## Org


    gwf = sim.get_model(modelDict['gwfname']) ## Org

    #generate model geometry as vtk and parameter array
    mf6Vtk.generateGeometryArrays() ## Org

    #generate parameter vtk
    mf6Vtk.generateParamVtk() ## Org

    #generate bc and obs vtk
    mf6Vtk.generateBcObsVtk(nper=0) ## Org

    mf6Vtk.generateHeadVtk(nper=0, crop=True) ## Org

    mf6Vtk.generateWaterTableVtk(nper=0) ## Org

