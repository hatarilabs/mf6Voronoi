import os, re, time
import flopy
import sys
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from scipy.interpolate import griddata
from mf6Voronoi.utils import isRunningInJupyter, printBannerHtml, printBannerText

def FlowVectorGenerator(gwf, backgroundImageDict=None, 
                        xul=None, yul=None, rotation=None,
                        kstpkper=(0,0),
                        plotGrid = True,
                        plotContour = True,
                        contourLevels = 10,
                        layer = 0,
                        istep = 4,
                        jstep = 4,
                        scale = 10,
                        normalize=True
                        ):

    if backgroundImageDict != None:
        try: 
            figPath = backgroundImageDict['fig'] 
            wrlPath = backgroundImageDict['wrl']
        except KeyError:
            print('ERROR: Something went wrong with the background image')
        #open figure and world file:
        bcgImg = mpimg.imread(figPath)
        wrlFile = open(wrlPath).read().split('\n')
        cellsize = float(wrlFile[0])
        left = float(wrlFile[4])
        right = float(wrlFile[4]) + bcgImg.shape[1]*cellsize
        bottom = float(wrlFile[5]) - bcgImg.shape[0]*cellsize
        top = float(wrlFile[5])      

    # Build geometry vtk
    if gwf.modelgrid.grid_type == 'structured' or gwf.modelgrid.grid_type == 'vertex':
        dis = gwf.get_package('dis')
        try:
            xul = dis.xorigin.data
            yul = dis.yorigin.data
            rotation = dis.angrot.data

            head = gwf.output.head().get_data(kstpkper)
            bud = gwf.output.budget()
            spdis = bud.get_data(text='DATA-SPDIS')[kstpkper[1]]

            hactive = head[head>gwf.hdry]
            hactive = head[head<gwf.hnoflo]
            levels = np.linspace(hactive.min(),hactive.max(),contourLevels)

            fig, ax = plt.subplots(figsize=(10, 10)) ## Org

            qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

            arrowColor = 'steelblue'

            #initialize graph
            pmv = flopy.plot.PlotMapView(gwf, layer=layer)
            #add layers
            pmv.plot_grid(colors='turquoise', lw=0.3, alpha=0.5,ax=ax, zorder=1)
            
            if plotGrid:
                pmv.plot_array(head, masked_values=[1e+30], cmap='YlGnBu',
                               ax=ax, zorder=2, alpha=0.5)
                arrowColor = 'azure'
            elif backgroundImageDict:
                arrowColor = 'azure'

            if plotContour:
                mf6Contour = pmv.contour_array(head, levels=levels, linewidths=3.,
                                               cmap='YlGnBu',
                                               ax=ax, zorder=3)
                texts = plt.clabel(mf6Contour, inline=True,
                                   fontsize=8, fmt="%.1f")

                for text in texts:
                    text.set_color('deepskyblue')
                    text.set_path_effects([
                        PathEffects.withStroke(linewidth=2, foreground='white')
                    ])
            pmv.plot_vector(qx, qy, normalize=normalize, color=arrowColor, 
                            scale_units='width',scale=scale, 
                            istep=istep,
                            jstep=jstep,
                            ax=ax, zorder=4)
            if backgroundImageDict != None:
                image = ax.imshow(bcgImg, extent=(left, right, bottom, top),
                                  zorder=0, alpha=0.7)


        except AttributeError:
            print('ERROR: Something went wrong with model data, check the dis file')

        
    # elif gwf.modelgrid.grid_type == 'vertex':
    #     pass
    elif gwf.modelgrid.grid_type == 'unstructured':
        print("DISU: This dicretization type is not supported")
    else:
        print("No Dis file was found")

    budget = gwf.output.budget()
    budgetData = budget.get_data(kstpkper=kstpkper)