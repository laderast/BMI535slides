# What Data Analysts Need to Know About Docker

![](image/containers-vs-virtual-machines.jpg)
https://www.weave.works/blog/a-practical-guide-to-choosing-between-docker-containers-and-vms

## What is Docker?

Docker uses *containers*. Containers contain just enough linux utilities and associated code to run that code on any machine. Some of them contain command line utilities, some of them can contain web servers/IDEs, some have databases. These containers are managed by *Docker Daemon*, a service that runs on a machine.

They contain a limited file system that is isolated from your system, but you can connect other parts of your file system to them with *Volumes*.

> Docker lets you package software and its dependencies in a way that when you run it, the same output can be reproduced across any machine that runs Docker, regardless of OS/processor architecture.

## What problems are Docker meant to Solve?

> The key benefit of Docker is that it allows users to **package an application with all of its dependencies into a standardized unit** for software development. - https://docker-curriculum.com/#what-is-docker-

- **Reproducibility** - Can precisely specify which versions of packages and software to use, so you can *replicate analysis in a reproducible manner*.
- **Containers are OS independent** - can run a docker container on Linux, Mac OS, Windows and it will work
- **Containers are lightweight** - require much less system overhead and resources
- **Software administration** - containers are independent of each other (very good for system administrators)
    - Gives you granular control over the software
    - Can take down one service without affecting the others

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

- **Images** - The blueprints of our application which form the basis of containers.
- **Containers** - Created from Docker images and run the actual application. We create a container using `docker run` based on an image. A list of running containers can be seen using the `docker ps` command. When a container is stopped using `docker stop`, all files in it fail to persist.
- **Docker Daemon** - The background service running on the host that manages building, running and distributing Docker containers. The daemon is the process that runs in the operating system which clients talk to.
- **Docker Client** - The command line tool that allows the user to interact with the daemon. More generally, there can be other forms of clients too - such as Kitematic which provide a GUI to the users.
- **Docker Hub** - A registry of Docker images. You can think of the registry as a directory of all available Docker images. If required, one can host their own Docker registries and can use them for pulling images.
- **Volume** - A File Directory that persists beyond a container.
- **Docker Compose** - utility that lets you connect multiple containers with a volume and with each other. 
- **Kubernetes** - utility that lets you run Docker (and other container systems) on a HPC cluster. Handles distributing an application over an allocation.

Some of these definitions were adapted from: https://docker-curriculum.com/

## Containers versus images

A **Docker Image** is a collection of software and its dependencies that can be pulled from DockerHub or another Registry. It is *read only* - you cannot modify the contents.

A **Docker Container** is an instance of the image running on a machine. 

Multiple containers can be run from a single docker image, and these containers will be independent of each other.

More information here: https://phoenixnap.com/kb/docker-image-vs-container

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

For more information about the Docker container:

https://hub.docker.com/r/broadinstitute/gatk

To download the image, I ran:

```
docker pull broadinstitute/gatk
```

## `docker run`

We will be using `docker run` to use our container. It works like a unix shell - we will open up an interactive terminal and then run `gatk` code inside the container.

The way you use a container varies - you'll need to check the documentation for each container. 

Some of them run as `services`, such as the `rocker` containers, which runs in the background on your machine and you'll access them via a your browser.

## Docker Volumes


![Docker Volumes](image/docker-volumes.jpg)
https://codeblog.dotsandbrackets.com/persistent-data-docker-volumes/


Volumes let you use files and folders outside of the built in container file system. In our `docker run` statement, we'll make the following volume.

```
-v ~/BMI535:/gatk/my_data
```

Everything to the left of the colon (`~/BMI535`) is a file on your system that you want to map to the container. This needs to be an *absolute* path.

Everything to the right of the colon (`/gatk/my_data`) is the *Volume* in the container's filesystem that we want to map to.

You can specify multiple volumes using multiple `-v` flags:

```
-v ~/BMI535:/gatk/my_data
-v ~/var_data/:/gatk/var_data
```

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

Which emans you can run `gatk` as usual in the container. Any files you write to `my_data` will persist.

Make sure that you've converted your `sam` file to a `bam` file, and sorted it before running `MarkDuplicates`

```{r}
gatk MarkDuplicates -I sortedSRR702072.bam -O SRR702072marked-dupes.bam -M metrics.txt
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
FROM rocker/r-ver:3.4.4    ### Build on top of the rocker/r-ver image with tag "3.4.4"

RUN R -e "install.packages('stringr')"  ##Run the system command R -e "install.packages"
```

We won't be building our own Dockerfiles, but we will leverage the built in system in `mybinder.org` to build our Docker containers.

Once you have a Dockerfile, you can build the image using

```
docker build
```

If you're in the folder.

## Dockerfile Example

```
FROM rocker/binder:latest

## Declares build arguments
ARG NB_USER
ARG NB_UID

## Copies your repo files into the Docker Container
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}

## Become normal user again
USER ${NB_USER}

## Run an install.R script, if it exists.
RUN if [ -f install.R ]; then R --quiet -f install.R; fi 
```

## What is a Binder?

![](image/binder_example.jpg)

A binder is a unit of reproducible research. If I give you a binder link that contains my code, you will be able to reproduce all of the analyses that it contains.

A binder usually lives in a repository, usually GitHub. 

## mybinder.org: a place for sharing Reproducible Research

![](image/binder.jpg)


We'll be using `mybinder.org` to test out a reproducible notebook. It uses a utility called `repo2docker` that converts a github repository to a Docker Container.

You can share your analyses in multiple ways: a Jupyter Notebook, an RStudio Project, or even a Shiny application. 

`mybinder.org` uses donated compute time. You are limited to 1 GB memory and 40 Gb of disk space.

## Some Lessons We learned in Starting a Binder

https://github.com/biodev/HNSCC_Notebook

[A Practical Guide to Reproducible Papers](https://biodata-club.github.io/talks/repro_paper.pdf)


- Requires a lot of distributed expertise (tracking efforts/who did what)
    - Code review should be done within your group
- Code will still be very messy! (that’s ok)
- Making a reproducible container involves a lot of crying and frustration
- Binder can break (can be difficult to maintain)
- Try to version software using version tags


# Acknowledgements

https://docker-curriculum.com/
https://biodata-club.github.io/talks/repro_paper.pdf
