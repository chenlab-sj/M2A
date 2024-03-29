#!/usr/bin/env bash
# Builds resources and then invokes dx build to build the DNAnexus app.
#
# This invokes docker and dx, so be sure they are in your path.
#
# Usage: ./build.sh [--no-docker-build] VERS [DX_BUILD_ARGS...]
#
# Parameters:
# --no-docker-build  if specified, does not perform docker build
# VERS:              the version you are trying to build
# DX_BUILD_ARGS      these arguments are passed to dx build
#
# The build performs the following operations:
# 1. Checks the dxapp.json to confirm the version matches VERS (fails if not)
# 2. Performs docker build, unless --no-docker-build is specified
# 3. Saves the docker image, gzipped, to the resources directory
# 4. Performs dx build

# Print usage instructions if no arguments were given
if [ "$#" == 0 ]; then sed -q '2,/^$/p' $0; exit 0; fi

# Get arguments
if [ "$1" == "--no-docker-build" ]
then no_docker_build=1; shift
fi
VERS=$1
shift
# Remaining arguments are passed to dx build

# Validate version
dxapp_version=`cat dxapp.json | tr -d '\t "' | sed 's/,$//' | awk -F : '$1 == "version" { print $2; exit }'`
if [ "$VERS" != "$dxapp_version" ]
then
  echo "Version from dxapp.json is $dxapp_version, not $VERS" >&2
  exit 1
fi

set -ex

if [ "$no_docker_build" == "" ]
then docker build -t stjude/m2a:$VERS ..
fi

mkdir -p resources/stjude
#docker save stjude/m2a:$VERS | gzip -c > resources/stjude/m2a-docker.tar.gz

if [ -d "resources/stjude/cwl" ]
then
  rm -r resources/stjude/cwl
fi 
mkdir -p resources/stjude/cwl
cp -r ../cwl/* resources/stjude/cwl/

dx build -a "$@"
