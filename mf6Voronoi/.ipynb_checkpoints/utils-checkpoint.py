import geopandas as gpd
import os
import folium
import io
import fiona
import json
import pathlib
import numpy as np
import rasterio
import rasterio.mask
import tempfile

from django.http import HttpResponseRedirect, HttpResponse, Http404, JsonResponse
from django.views.generic import View, FormView
from django.contrib.auth.mixins import LoginRequiredMixin
from mf6.home.models import Project
import shutil
from django.core.exceptions import ValidationError


def read_shp_from_zip(file):
    zipshp = io.BytesIO(open(file, 'rb').read())
    with fiona.BytesCollection(zipshp.read()) as src:
        crs = src.crs
        gdf = gpd.GeoDataFrame.from_features(src, crs=crs)
    return gdf


def write_shp_as_zip(zipLoc,zipDest,baseName):
    shutil.make_archive(base_dir=zipLoc,
        root_dir=zipDest,
        format='zip',
        base_name=baseName)

def shpFromZipAsFiona(file):
    zipshp = io.BytesIO(open(file, 'rb').read())
    fionaObj = fiona.BytesCollection(zipshp.read())
    return fionaObj


def remove_files_and_folder(path_to_file, folder=True):
    folder=os.path.dirname(path_to_file)

    if os.path.isfile(path_to_file):
        os.remove(path_to_file)
        print("File has been deleted")
    else:
        print("File does not exist")
    #for filename in os.listdir(folder):
    #    file_path=os.path.join(folder,filename)
        """
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    if folder:
        os.rmdir(folder)
    """

class LoguedContextedView(LoginRequiredMixin,View):
    #get tokens and filter projects of current user
    def getToken(self,request):
        current_path = request.get_full_path()
        token = current_path.split('/')[2]
        #project = Project.objects.get(token=token)
        if request.user == Project.objects.get(token=token).user:
            request.session['token'] = token
            #request.session['project'] = project
            return token
        else:
            raise Http404

    def getProject(self,request):
        token = self.getToken(self.request)
        return Project.objects.get(token=token)

    def dispatch(self, request, *args, **kwargs):
        self.getToken(request)
        #self.get_context_data(**kwargs)
        return super(LoguedContextedView, self).dispatch(request, *args, **kwargs)

    #for the namespaces in the sidebar
    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        ctx.update({
            'token':self.getToken(self.request)
            })
        return ctx

def file_size(value): # add this to some file where you can import it from
    limit = 4 * 1024 * 1024
    if value.size > limit:
        raise ValidationError('File too large, size should not exceed 2 MB.')
    else:
        return value

def add_to_map(geometry,mapa):
    if geometry['geometry'].geom_type == 'Point':
        lat = geometry['geometry'].y
        lon = geometry['geometry'].x
        folium.Marker(location=[lat,lon],popup=geometry['id']).add_to(mapa)
    else:
        sim_geo = gpd.GeoSeries(geometry['geometry'])#.simplify(tolerance=0.001)
        geo_j = sim_geo.to_json()
        style = {'fillColor': '#5AAACD',"weight": 1,"opacity": 0.65,"fillOpacity": 0}
        geo_j = folium.GeoJson(data=geo_j,
                           style_function=lambda x: style)#{'fillColor': 'orange'})
        try:
            folium.Popup(str(geometry['id']),permanent=True,parse_html=True).add_to(geo_j)
        except:
            pass
        geo_j.add_to(mapa)

    return 0

def save_json(dbmodel,project,bc_name,bc_DataDict):
    #function to save json for boundaries conditions
    #dbmodel = WellsBcJson
    #project, bc_name = json name for output file
    #bc_DataDict is the stress period for the bc, example well_spd
    try:
        dbmodel.objects.get(project=project).delete()
    except dbmodel.DoesNotExist:
        pass

    modelDict = dbmodel()
    modelDict.project = project

    outName = bc_name
    outAbsPath = os.path.join(modelDict.getAbsDir,outName)
    if not os.path.isdir(modelDict.getAbsDir):
        os.makedirs(modelDict.getAbsDir,exist_ok=True)
    #bc_DataFrame.to_json(outAbsPath)
    with open(outAbsPath,"w") as write_file:
        json.dump(bc_DataDict,write_file)

    modelDict.json_File.name = os.path.join(modelDict.getRelDir,outName)
    modelDict.save()

    return 0

def save_flopy(file_path,code_lines_list):
    file = open(file_path,'a')
    for code in code_lines_list:
        file.write(code)
    file.close()
    return 0

def getXY(src):
    data = src.read(1)
    height,width = data.shape
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rasterio.transform.xy(src.transform, rows, cols)
    lons = [item.tolist() for item in xs]
    lats = [item.tolist() for item in ys]
    return lons,lats

def positive(array):
    pos_count,neg_count = 0,0
    for row in array:
        for num in row:
            if num>=0:
                pos_count += 1
            else:
                neg_count += 1

    return pos_count*100/array.size


def clip_raster(shp_path,rst_path):
    # Open the shapefile containing the clipping polygon
    clip = gpd.read_file(shp_path)

    # Open the raster file to be clipped
    with rasterio.open(rst_path) as src:
        out_image, out_transform = rasterio.mask.mask(src, clip.geometry, crop=True)
        out_meta = src.meta.copy()

        with tempfile.NamedTemporaryFile(suffix='.tif', delete=False) as tmpfile:
            with rasterio.open(tmpfile.name, 'w', **src.profile) as dst:
                dst.write(out_image)

        tmpfile_path = tmpfile.name

    # Update the metadata of the output file
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    return out_image,tmpfile_path,out_transform,out_meta