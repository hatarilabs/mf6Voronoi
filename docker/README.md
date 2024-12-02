### Build mf6Voronoi Docker Image Locally

Clone PyVista and run the following at the top level of the project:

```bash
docker build -t mf6voronoi-jupyterlab -f docker/jupyter.Dockerfile .
docker run -p 8888:8888 --name mf6voronoi mf6voronoi-jupyterlab:latest
```

To remove the container and the image

```bash
docker rm -f mf6voronoi
docker image rm mf6voronoi-jupyterlab:latest
```

If you want to work with Docker compose with mounted files

```bash
docker build -t mf6voronoi-jupyterlab -f docker/jupyter.Dockerfile .
docker compose -f docker/dockercompose.yml up 
```