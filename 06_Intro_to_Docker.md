# Introduction to Docker and Docker Containers

## What is Docker?

Alternative to Virtual Machines, which contain an entire OS (Linux, Windows, or MacOS) in order to run software.

Containers contain just enough linux utilities and associated code to run that code on any machine. Some of them contain command line utilities, some of them can contain web servers/IDEs, some have databases. These containers are managed by *Docker Daemon*.

They contain a limited file system that is isolated from your system, but you can connect other parts of your file system to them with *Volumes*.

## What problems is Docker meant to Solve?

> The key benefit of Docker is that it allows users to **package an application with all of its dependencies into a standardized unit** for software development. - https://docker-curriculum.com/#what-is-docker-

- **Software administration** - containers are independent of each other (very good for system administrators)
    - Gives you granular control over the software
    - Can take down one service without affecting the others
- **Containers are OS independent** - can run a docker container on Linux, Mac OS, Windows and it will work
- **Reproducibility** - Can precisely specify which versions of packages and software to use, so you can *replicate analysis in a reproducible manner*.

## What software is containerized?

Almost all Bioinformatics Software and Web Stacks!

- RStudio Connect - very fast to spin up an instance of AWS (15 minutes)
- rocker project - contains RStudio and R Dependencies for specific kinds of analyses: https://www.rocker-project.org/
- Bioconductor: https://www.bioconductor.org/help/docker/
- GATK - https://gatk.broadinstitute.org/hc/en-us/articles/360035889991--How-to-Run-GATK-in-a-Docker-container
- Web Stacks (web server/scripting software/database system) - https://github.com/hanafiah/docker-webstack
- Databases - https://docs.docker.com/engine/examples/postgresql_service/

## Terminology

![](image/docker.png)

- **Images** - The blueprints of our application which form the basis of containers. In the demo above, we used the docker pull command to download the busybox image.
- **Containers** - Created from Docker images and run the actual application. We create a container using docker run which we did using the busybox image that we downloaded. A list of running containers can be seen using the docker ps command.
- **Docker Daemon** - The background service running on the host that manages building, running and distributing Docker containers. The daemon is the process that runs in the operating system which clients talk to.
- **Docker Client** - The command line tool that allows the user to interact with the daemon. More generally, there can be other forms of clients too - such as Kitematic which provide a GUI to the users.
- **Docker Hub** - A registry of Docker images. You can think of the registry as a directory of all available Docker images. If required, one can host their own Docker registries and can use them for pulling images.
- **Volume** - A File Directory that persists beyond a container
- **Docker Compose** - utility that lets you connect multiple containers with a volume and with each other. 
- **Kubernetes** - utility that lets you run Docker (and other container systems) on a HPC cluster.


## Drawbacks to Docker

Docker requires **system-level access** to your machine (the equivalent of root access), so you can't install it by yourself on a remote machine where you don't have access privileges. You'll need to ask adminstration to install Docker and add you to the `docker` group. 

**Everything inside a container disappears when you finish running the container**. That means that files that you add to a container will disappear after you run them. You fix this by mounting **Volumes** to your container. 

## Docker Images

Docker containers can be pulled from a number of sources, but DockerHub is the most common.

On `state`, you should have access to the `gatk` container already. Check it exists by running

```
docker images
```

confirm that you can see `broadinstitute/gatk` 

## `docker ps`

Check what is running using 

```
docker ps
```

## Broad Institute GATK Docker Image

https://hub.docker.com/r/broadinstitute/gatk

```
docker pull broadinstitute/gatk
```

## Docker Volumes

Volumes let you use files and folders outside of the built in container file system.


```
-v ~/BMI535:/gatk/my_data
```

Everything to the left of the colon (`~/BMI535`) is a file on your system that you want to map to the container. This needs to be an *absolute* path.

Everything to the right of the colon (`/gatk/my_data`) is the *Volume* in the container's filesystem that we want to map to.

## Using the GATK Docker Image

The main command you use to run docker is `docker run`. Here we're running 

```{r}
docker run -v ~/BMI535:/gatk/my_data -it broadinstitute/gatk
```

The `-it` flag opens an interactive terminal into the container. This will open a limited version of Unix.
The `-v` flag mounts your `~/BMI535/` folder here in a location called 

```
/gatk/my_data/
```

You'll be in the `/gatk/` folder when you open. Confirm that you can see the `my_data` folder:

```
ls
cd my_data
```

Which emans you can run `gatk` as usual in the container. Any files you write to `my_data` will persist

```{r}
gatk MarkDuplicates 
```

When you're done with the container, use

```
exit 
```

To get out of it.

## Making your own Docker containers: Dockerfiles

Any file in a folder called `Dockerfile` will be used to build a container. This contains instructions for installing all of the software and its dependencies.

The good news is that you don't have to build your Dockerfiles from scratch.

You can build on a previous Docker container by using a `FROM` command.

```
FROM rocker/r-ver:3.4.4    ### Build on top of the rocker/r-ver image

RUN R -e "install.packages('stringr')"  ##Run the system command R -e "install.packages"
```

We won't be building our own Dockerfiles, but we will leverage the built in system in `mybinder.org` to build our Docker containers.

Once you have a Dockerfile, you can build the image using

```
docker build
```

If you're in the folder.

## mybinder.org: a place for sharing Reproducible Research

[Slides introducing mybinder](https://docs.google.com/presentation/d/1y2HrtsmERC9hfliJkWpJVQ28mzJNmoWNW44XizXaHok/edit?usp=sharing)

We'll be using `mybinder.org` to test out a reproducible notebook. It uses a utility called `repo2docker` that converts a github repository to a Docker Container.

You can share your analyses in multiple ways: a Jupyter Notebook, an RStudio Project, or even a Shiny application. 

`mybinder.org` uses donated compute time. You are limited to 1 GB memory and 40 Gb of disk space.

# Acknowledgements

https://docker-curriculum.com/
