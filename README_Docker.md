Containers
==========
With this release, I am including infrastructure to help those using
Docker containers.  This may help people who have difficulty building
the software from scratch - although you will need to install and
configure docker instead
[docker](https://docs.docker.com/get-docker/).

In the distro, there is a Dockerfile which I have used to generate a Docker
image that has been pushed to Docker Hub. Retrieve to your local Docker
environment via

```
docker pull ianwatson/lilly_medchem_rules:v1.2
```

The image is built on a gcc image, then goes to some trouble to
install Ruby, via `rvm`.  Perhaps a better approach might be to start
with a Ruby image and add gcc to that.

The image has a CMD that will consume a smiles file via stdin and
output surviving molecules to stdout.
```
docker run -i -a stdout -a stdin ianwatson/lilly_medchem_rules:v1.2 < file.smi
```

There is a script in the container `run_from_docker_image.sh` that will
mount a local file directory containing a smiles file, into the image
and run the rules on that input file. The bad*.smi files and the
log files will be written to that same directory, while the surviving
molecules will be written to stdout. Easy would be to copy that file
to somewhere outside the image and use it there.

I would welcome any feedback on this initial container implementation,
since I am no expert. But I do want to provide an alternative for people
encountering build difficulties.
