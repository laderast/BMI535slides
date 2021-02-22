# Introduction to Docker and Docker Containers

## What is Docker?

Alternative to Virtual Machines, which contain an entire OS (Linux, Windows, or MacOS) in order to run software in them.

Containers contain just enough linux utilities and associated code to run that code on any machine. Some of them contain command line utilities, some of them can contain web servers/IDEs.

They contain a limited file system that is isolated from your system, but you can connect other parts of

## What problems is Docker meant to Solve?

> The key benefit of Docker is that it allows users to **package an application with all of its dependencies into a standardized unit** for software development. - https://docker-curriculum.com/#what-is-docker-

- **Software administration** - containers are independent of each other
    - Gives you granular control over the software.
- **Containers are OS independent** - can run a docker container on Linux, Mac OS, Windows and it will work
- **Reproducibility** - Can precisely specify which versions of packages and software to use, so you can *replicate analysis in a reproducible manner*


## What software is containerized?

Almost all Bioinformatics Software and Web Stacks

- rocker project - contains RStudio and R Dependencies for specific kinds of analyses: https://www.rocker-project.org/
- Bioconductor: https://www.bioconductor.org/help/docker/
- GATK - https://gatk.broadinstitute.org/hc/en-us/articles/360035889991--How-to-Run-GATK-in-a-Docker-container
- Web Stacks (web server/scripting software/database system) - https://github.com/hanafiah/docker-webstack

## Terminology

- **Container** - contains application software, dependencies, and linux.  "By leveraging the low-level mechanics of the host operating system, containers provide most of the isolation of virtual machines at a fraction of the computing power.""
- **Volume** - A File Directory that persists beyond a container
- **Docker Compose** - utility that lets you connect multiple containers with a volume and with each other. 
- **Kubernetes** - utility that lets you run Docker (and other container systems) on a HPC cluster.

## Drawbacks to Docker

Docker requires **system-level access** to your machine (the equivalent of root access), so you can't install it by yourself on a remote machine where you don't have access privileges. You'll need to ask adminstration to install Docker and add you to the `docker` group. 

**Everything inside a container disappears when you finish running the container**. That means that files that you add to a container will disappear after you run them. You fix this by mounting **Volumes** to your container. 

## Docker Images

Docker containers can be pulled from a number of sources, but DockerHub is the most common.

On `state`, you should have access to the `gatk` container already.

## Broad Institute GATK Docker Image

https://hub.docker.com/r/broadinstitute/gatk

```
docker pull broadinstitute/gatk
```

## Using the GATK Docker Image

```{r}
docker run -v ~/BMI535:/gatk/my_data -it broadinstitute/gatk
```

This will open a limited version of Unix and mounts your `~/BMI535/` folder here in a location called 

```
/gatk/my_data/
```

Which emans you can run `gatk` as usual.


## Making your own Docker containers: Dockerfiles

