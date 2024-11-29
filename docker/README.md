### Build mf6Voronoi Docker Image Locally

Clone PyVista and run the following at the top level of the project:

```bash
docker build -t mf6voronoi-jupyterlab -f docker/jupyter.Dockerfile .
```