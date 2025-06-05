import flopy ## Org
from mf6Voronoi.tools.graphs2d import FlowVectorGenerator
import matplotlib.pyplot as plt
import geopandas as gpd
import os, json

with open('flowDirectionVectorsData.json') as jsonFile:
    flowDirectionVectorsDict = json.load(jsonFile)
    
for modelName, modelDict in flowDirectionVectorsDict.items():

    print("\n Working for model: %s"%modelName)

    # load simulation
    simName = 'mf6Sim' ## Org

    modelWs = os.path.abspath(os.path.join('../../../examples/mf6models',modelName))

    sim = flopy.mf6.MFSimulation.load(sim_name=simName, version='mf6', ## Org
                                exe_name='../../../../bin/mf6', ## Org
                                sim_ws=modelWs) ## Org
    
    gwf = sim.get_model(modelDict['gwfname']) ## Org

    if modelDict['backgroundImage']:
        backgroundPath = os.path.join('../../../examples/datasets',modelName,'Png')
        backgroundImageDict = {
            'fig':os.path.join(backgroundPath,'backgroundImage.png'),
            'wrl':os.path.join(backgroundPath,'backgroundImage.pgw')
        }
        fig = FlowVectorGenerator(gwf, backgroundImageDict=backgroundImageDict,
                                  layer=modelDict['layer'], plotGrid=True)
    else:
        fig = FlowVectorGenerator(gwf, layer=modelDict['layer'], plotGrid=True)

    figPath='../../checkFiles/postProcessing/flowDirectionVectors/'+modelName+'.png'
    fig.savefig(figPath)