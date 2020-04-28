#!/bin/bash

# Runs the Lilly Medchem Rules using the executables within the
# official container.
# Run with a single argument, which should be a smiles file to
# process.  Output will be to stdout.  The bad*.smi files and log
# files will be written to the same directory as the input file. Yes,
# the container writes to the outside file system, via the mount below.

container='ianwatson/lilly_medchem_rules:v1.2'

if [[ -z "$@" ]] ; then
  echo "Must specify a file to run" >&2
  exit 1
fi

path=$(readlink -e $@)
if [[ -z "${path}" ]] ; then
  echo "Missing or invalid file ${@}" >&2
  exit 1
fi

outside_directory=$(dirname $path)
outside_fname=$(basename $path)

# We will mount $outside_directory at $docker_dir inside the container.
docker_dir=/mutable/outside/world

docker run --rm --mount type=bind,source=$outside_directory,destination=${docker_dir} \
--entrypoint "bash" ${container} \
-c ". /etc/profile && rvm use 2.7.1 > /dev/null && cd ${docker_dir}\
   && /Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb ${outside_fname}"
