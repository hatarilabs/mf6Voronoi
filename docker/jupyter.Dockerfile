FROM jupyter/base-notebook:python-3.11
LABEL maintainer="mf6Voronoi Developers"
LABEL repo="https://github.com/hatarilabs/mf6Voronoi"

COPY docker/requirements.txt /build-context/requirements.txt
WORKDIR /build-context/

RUN pip install -r requirements.txt

WORKDIR $HOME

