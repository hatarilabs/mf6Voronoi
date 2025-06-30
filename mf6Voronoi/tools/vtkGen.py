import os, re, time
import flopy
import sys
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.interpolate import griddata
from mf6Voronoi.utils import isRunningInJupyter, printBannerHtml, printBannerText

class Mf6VtkGenerator:
    def __init__(self, sim, vtkDir):
        self.sim = sim
        self.vtkDir = vtkDir

        #insert banner
        if isRunningInJupyter():
            printBannerHtml()
        else:
            printBannerText()

        print('\n/---------------------------------------/')
        print('\nThe Vtk generator engine has been started')
        print('\n/---------------------------------------/')

    def listModels(self):
        print("\n Models in simulation: %s"%self.sim.model_names)
        
    def loadModel(self, modelName):
        self.gwf = self.sim.get_model(modelName)
        self.packageList = self.gwf.get_package_list()
        print("Package list: %s"%self.packageList)
    
    def saveArray(self,outputDir, array, arrayName):
        if array is not None:
            np.save(os.path.join(outputDir,'param_'+arrayName),array)    
            
    def timing_decorator(func):
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            execution_time = end_time - start_time
            print(f"Vtk file took {execution_time:.4f} seconds to be generated.")
            return result
        return wrapper
    
    @timing_decorator
    def exportBc(self,bcon,nper):
        #open package and filter names
        bcObj = self.gwf.get_package(bcon)
        if bcObj.has_stress_period_data: 
            bcObjSpdNames = bcObj.stress_period_data.dtype.names[1:]
            print('Working for %s package, creating the datasets: %s'%(bcon,bcObjSpdNames))
            #create a flat index
            flatIndexList = []
            # flatIndexTupleList = []

            #get the first stress period that has the active bc
            bcKeys = list(bcObj.stress_period_data.data.keys())
            
            try:
                tempVtkGeom = self.vtkGeom.copy()
                nCells = tempVtkGeom.n_cells
                tempVtkGeom.cell_data['cell_id'] = np.arange(nCells)
                for name in bcObjSpdNames:
                    tempVtkGeom.cell_data[name] = np.zeros(nCells)

                try: 
                    for row in bcObj.stress_period_data.data[nper]:
                        #print(row.cellid)
                        flatIndex = np.ravel_multi_index(row.cellid,self.gwf.modelgrid.shape) #works for dis and disv
                        flatIndexList.append(flatIndex)
                    
                    #working with EVT that has no dataframe option
                    try:
                        bcSpdDf = bcObj.stress_period_data.dataframe[nper]
                    except AttributeError:
                        auxRec = bcObj.stress_period_data.get_data(nper)
                        bcSpdDf = pd.DataFrame.from_records(auxRec)

                    for index, row in bcSpdDf.iterrows():
                        if 'cellid_row' in bcSpdDf.columns:
                            cellid = (int(row.cellid_layer), int(row.cellid_row), int(row.cellid_column))
                        elif 'cellid_cell' in bcSpdDf.columns:
                            cellid = (int(row.cellid_layer), int(row.cellid_cell))
                        elif 'cellid' in bcSpdDf.columns:
                            cellid = row.cellid
                        else:
                            print('[Error] Something went wrong with cell indexing')

                        flatIndex = np.ravel_multi_index(cellid,self.gwf.modelgrid.shape) 
                        for name in bcObjSpdNames:
                            if name != '':
                                if not isinstance(row[name],str):
                                    tempVtkGeom[name][flatIndex] = row[name]
                                else:
                                    pass
                                    #print('[WARNING] The following dataset %s is a string and will be filled by zeros'%name)

                    tempVtkFilter = tempVtkGeom.extract_cells(flatIndexList) 
                    tempVtkFilter.save(os.path.join(self.vtkDir,'%s_kper_%s.vtk'%(bcon,nper)))      
                except KeyError:
                    print('[WARNING] There is no data for the required stress period')
                

            except IndexError:
                print('[WARNING] There is no data for the required stress period')
        elif 'TS' in bcon:
            print('[WARNING] This boundary condition %s is a time series and wont be considered'%bcon)
        else:
            bcObjBlkNames = bcObj.blocks['period'].datasets.keys()
            print('Working for %s package, creating the datasets: %s'%(bcon,bcObjBlkNames))
            tempVtk = self.vtkGeom.extract_cells(range(self.gwf.modelgrid.ncpl))
            try:
                for name in bcObjBlkNames:
                    dataSet = bcObj.blocks['period'].datasets[name]
                    if dataSet.has_data():        
                        tempVtk.cell_data[name] = dataSet.get_data(nper).flatten()
                tempVtk.save(os.path.join(self.vtkDir,'%s_kper_%s.vtk'%(bcon,nper)))
            except AttributeError:
                print('[WARNING] There is no data for the required stress period')
        

        
    @timing_decorator
    def exportObs(self,obs):
        #open package 
        obsObj = self.gwf.get_package(obs)
        obsKey = list(obsObj.continuous.data.keys())[0]
        obsObjNames = obsObj.continuous.data[obsKey].dtype.names[:2]
        bcObjArray = obsObj.continuous.data[obsKey].id
        
        #create a flat index
        flatIndexList = []
        for cell in bcObjArray:
            flatIndex = np.ravel_multi_index(cell,self.gwf.modelgrid.shape)
            flatIndexList.append(flatIndex)
        #empty object
        cropVtk = pv.UnstructuredGrid()

        for index, cell in enumerate(flatIndexList):
            tempVtk = self.vtkGeom.extract_cells(cell)
            for name in obsObjNames:
                tempVtk.cell_data[name] = np.array([obsObj.continuous.data[obsKey][name][index]])
            cropVtk += tempVtk

        #save and return
        cropVtk.save('../Vtk/%s.vtk'%obs)
        #return cropVtk
        
    def generateGeometryArrays(self):
        dis = self.gwf.get_package('DIS')
        dis.export(self.vtkDir, fmt='vtk', smooth=True, binary=False)
        
        if self.gwf.get_package('DIS') is not None:
            dis = self.gwf.get_package('DIS')
            idomainArray = dis.idomain.array
            self.saveArray(self.vtkDir, idomainArray, 'idomain')

        if self.gwf.get_package('IC') is not None:
            ic = self.gwf.get_package('IC')
            strtArray = ic.strt.array
            self.saveArray(self.vtkDir, strtArray, 'strt')

        if self.gwf.get_package('STO') is not None:
            sto = self.gwf.get_package('STO')
            iconvertArray = sto.iconvert.array
            syArray = sto.sy.array
            ssArray = sto.ss.array
            self.saveArray(self.vtkDir, iconvertArray, 'iconvert')
            self.saveArray(self.vtkDir, syArray, 'sy')
            self.saveArray(self.vtkDir, ssArray, 'ss')

        if self.gwf.get_package('NPF') is not None:
            npf = self.gwf.get_package('NPF')
            icelltypeArray = npf.icelltype.array
            kArray = npf.k.array
            k22Array = npf.k22.array
            k33Array = npf.k33.array
            wetdryArray = npf.wetdry.array
            self.saveArray(self.vtkDir, icelltypeArray, 'icelltype')
            self.saveArray(self.vtkDir, kArray, 'k')
            self.saveArray(self.vtkDir, k22Array, 'k22')
            self.saveArray(self.vtkDir, k33Array, 'k33')
            self.saveArray(self.vtkDir, wetdryArray, 'wetdry')
        
        # Build geometry vtk
        if self.gwf.modelgrid.grid_type == 'structured':
            vtkFile = os.path.join(self.vtkDir,"dis.vtk")
        elif self.gwf.modelgrid.grid_type == 'vertex':
            vtkFile = os.path.join(self.vtkDir,"disv.vtk")
        elif self.gwf.modelgrid.grid_type == 'unstructured':
            vtkFile = os.path.join(self.vtkDir,"disu.vtk")
        else:
            print("No Dis file was found")
        gwfGrid = pv.read(vtkFile)
        gwfGrid.clear_cell_data()
        #self.geomPath = ('../Model/Vtk/modelGeometry.vtk')
        self.geomPath = os.path.join(self.vtkDir,'modelGeometry.vtk')
        gwfGrid.save(self.geomPath)
        os.remove(vtkFile)
        self.vtkGeom = gwfGrid
        
    def generateParamVtk(self):
        vtkParam = pv.read(self.geomPath)
        paramList = [file for file in os.listdir(self.vtkDir) if file.startswith('param')]
        
        for param in paramList:
            paramName = re.split('[_;.]',param)[1]
            paramValues = np.load(os.path.join(self.vtkDir,param),allow_pickle=True)
            vtkParam.cell_data[paramName] = paramValues.flatten()
        paramPath = os.path.join(self.vtkDir,'modelParameters.vtk')
        vtkParam.save(paramPath)
        print("Parameter Vtk Generated")
        
    def generateBcObsVtk(self, nper):
        vtkGrid = pv.read(self.geomPath)
        bcList = [x for x in self.packageList if re.search(r'\d',x) and not re.search('obs',x,re.IGNORECASE)]
        obsList = [x for x in self.packageList if re.fullmatch('obs',x,re.IGNORECASE)]
        #print(bcList)
        for bc in bcList:
            print('\n/--------%s vtk generation-------/'%bc)
            self.exportBc(bc, nper)
            print('/--------%s vtk generated-------/\n'%bc)
        for obs in obsList:
            self.exportObs(obs)
            print("%s btk generated"%obs)

    def generateHeadVtk(self, nper, nstp=0, crop=False):
        headObj = self.gwf.output.head()
        heads = headObj.get_data(kstpkper=(nstp,nper))
        waterTable = flopy.utils.postprocessing.get_water_table(heads).flatten()

        geomVtk = pv.read(self.geomPath)
        geomVtk.clear_cell_data()
        geomVtk.cell_data['heads'] = heads.flatten()

        if crop:
            botmArray = self.gwf.modelgrid.botm.flatten()
            #wtArray = np.repeat(waterTable,self.gwf.modelgrid.nlay).filled()
            wtArray = np.hstack([waterTable for i in range(self.gwf.modelgrid.nlay)]) 
            activeArray = np.where(wtArray > botmArray, 1, 0)
            geomVtk.cell_data['active'] = activeArray

        # geomVtk = geomVtk.threshold(value=(1,1), scalars="heads")
        # geomVtk = geomVtk.cell_data_to_point_data()
        geomVtk.save(os.path.join(self.vtkDir,'waterHeads_kper_%s.vtk'%nper))

    def generateWaterTableVtk(self, nper, nstp=0):
        headObj = self.gwf.output.head()
        heads = headObj.get_data(kstpkper=(nstp,nper))
        waterTable = flopy.utils.postprocessing.get_water_table(heads).flatten()

        grid = self.gwf.modelgrid

        xCell = grid.xcellcenters.flatten()
        yCell = grid.ycellcenters.flatten()
        zCell = waterTable.data

        xyVerts = grid.verts

        # Interpolate z-values at new xy locations
        zVert = griddata((xCell, yCell), zCell, xyVerts, method='nearest')

        # Working with grid type
        if type(self.gwf.modelgrid) == flopy.discretization.vertexgrid.VertexGrid:            
            faces = np.hstack([cell[3:] for cell in grid.cell2d])
        elif type(self.gwf.modelgrid) == flopy.discretization.structuredgrid.StructuredGrid:
            nrowvert = grid.nrow + 1
            ncolvert = grid.ncol + 1
            npoints = nrowvert * ncolvert
            iverts = []
            for i in range(grid.nrow):
                for j in range(grid.ncol):
                    iv1 = i * ncolvert + j  # upper left point number
                    iv2 = iv1 + 1
                    iv4 = (i + 1) * ncolvert + j
                    iv3 = iv4 + 1
                    iverts.append([4, iv1, iv2, iv3, iv4])  
            faces = np.array(iverts).flatten() 
        else:
            print('Your grid type is not supported')
        
        points = list(zip(grid.verts[:,0],grid.verts[:,1],zVert))
        mesh = pv.PolyData(points, faces)

        mesh.cell_data['waterTable'] = waterTable

        mesh = mesh.cell_data_to_point_data()
        mesh.save(os.path.join(self.vtkDir,'waterTable_kper_%d.vtk'%nper))
