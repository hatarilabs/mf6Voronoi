from setuptools import setup, find_packages

setup(
    name='mf6Voronoi',
    version='0.1.0',
    description='A package to generate Voronoi grids for MODFLOW 6 DISV.',
    author='Hatarilabs',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'geopandas',
        'shapely',
        'rasterio',
        'matplotlib',
        'mapclassify', # for explore()
        'folium'       # for explore()
    ],
    extras_require={
        'test': ['pytest'],
    },
)
