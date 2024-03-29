#!/bin/bash
# M2A 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {

    echo "Value of curatedSites: '$curatedSites'"
    echo "Value of chipBigwig: '$chipBigwig'"
    echo "Value of inputBigwig: '$inputBigwig'"
    echo "Value of model: '$model'"
    echo "Value of ref_name: '$ref_name'"
    echo "Value of promoterDefinitions: '$promoterDefinitions'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".
    echo ""
    echo "=== Setup ==="
    echo "  [*] Downloading input files ..."
    dx download "$curatedSites" -o $curatedSites_name

    dx download "$chipBigwig" -o $chipBigwig_name

    dx download "$inputBigwig" -o $inputBigwig_name

    dx download "$model" -o $model_name


    # We allow the user to choose common  genomes from a dropdown
    # (ref_name), with one option being Custom, which needs to be accompanied
    # by a promoter defintions file being specified as an input
    # (promoterDefinitions).  Here, we pull the reference from the appropriate
    # place and determine the file basename.  If neither or both of the two
    # options are specified, we error out, as the reference cannot be
    # unambiguously determined in those cases.
    echo "  [*] Downloading reference files ..."
    
    promoter=
    if [ "$ref_name" != "" -a ${ref_name:0:6} != "Custom" ]
    then
      if [ "$promoterDefinitions" != "" ]
      then
        echo "Could not determine which genome to use." >&2
        echo "A custom promoter definitions was provided but Genome was not set to custom." >&2
        exit 1
      fi
      if [ "$ref_name" == "GRCh37-lite" ]
      then
        promoter="promoterDefinitions.GRCh37-lite.txt"
        dx download -o $promoter project-F5444K89PZxXjBqVJ3Pp79B4:file-Fq1f2X09gg7FKPGk3gfkvGfK
      elif [ "$ref_name" == "GRCh38" ]
      then
        promoter="promoterDefinitions.GRCh38.txt"
        dx download -o $promoter -r project-F5444K89PZxXjBqVJ3Pp79B4:file-Fq1f2VQ9gg7BjVv60vf535fV 
      fi
    else
      if [ "$promoterDefinitions" == "" ]
      then
        echo "Could not determine which genome to use." >&2
        echo "A custom promoter definitions was not provided but Genome was set to custom." >&2
        exit 1
      fi
      promoter=$promoterDefinitions_name
      dx download -o $promoter $promoterDefinitions
    fi

    echo "  [*] Loading container image ..."
    image_tarfile_path=/stjude/m2a-docker.tar
    if [ -e $image_tarfile_path.gz ]
    then gunzip $image_tarfile_path.gz
    fi
    docker load -i $image_tarfile_path

    echo "  [*] Printing python information ..." 
    which python3
    python3 --version
    #echo $PYTHONPATH
    origin_python_paths=$PYTHONPATH
    #python3 -m site
    #python3 -c "import site; print(site.getsitepackages())"

    echo "  [*] Resetting PYTHONPATH ..."
    PYTHONPATH=''
    path_lst=$(python3 -c "import site; print(site.getsitepackages())")
    for item in $path_lst
    do
        path=$(echo $item | sed "s/\[//g; s/'//g; s/,//g; s/\]//g")
        export PYTHONPATH=$path:$PYTHONPATH
    done
    #echo $PYTHONPATH

    echo "  [*] Installing dependencies ..."
    # Install python3, pip3 and cwltool
    apt-get update
    apt-get -y upgrade
    apt-get -y install python3 python3-pip
    apt-get -y install python-dev libxml2-dev libxslt1-dev zlib1g-dev virtualenv
    #virtualenv -p python3 venv
    #source venv/bin/activate
    #pip3 install --upgrade pip3
    pip3 install cwlref-runner
    cwltool --version

    #echo " [*] Installing miniconda ..."
    #wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" -O miniconda.sh 
    #/bin/bash miniconda.sh -b -p /opt/conda/  
    #/opt/conda/bin/conda create -n python3 python=3.7 cwlref-runner  
    #/opt/conda/bin/conda activate python3 
    
    echo "=== Execution ==="
    
    # Don't make assumptions about the tag that was used when the image was
    # built, other than that it should be "xenocp".  Since this is run in a
    # clean environment, and we only did a single docker load, the method
    # below should be a safe way to determine the image.
    echo "  [*] Determine image ID ..."
    image_id=`docker images -q stjude/m2a | head -n 1`
    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # Generate CWL input
    cat > inputs.yml <<EOF
chipBigwig:
  class: File
  path: $chipBigwig_name
inputBigwig:
  class: File
  path: $inputBigwig_name
curated: 
  class: File
  path: $curatedSites_name
promoterDefinitions:
  class: File
  path: $promoter
model: 
  class: File
  path: $model_name
EOF
    # Run M2A
    mkdir -p out
    cwl-runner --parallel --outdir out /stjude/cwl/m2a.cwl inputs.yml 


    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.
    mkdir -p out/predictions
    mv out/Predictions* out/predictions/
    #predictions=$(dx upload predictions --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.
    mkdir -p out/features
    mv out/Features*.txt out/features

    mkdir -p out/features_h5
    mv out/Features*.h5 out/features_h5

    #for i in "${!chromosomePromoters[@]}"; do
    #    dx-jobutil-add-output chromosomePromoters "${chromosomePromoters[$i]}" --class=array:file
    #done
    #dx-jobutil-add-output predictions "$predictions" --class=file
    export PYTHONPATH=$origin_python_paths
    dx-upload-all-outputs --parallel
    echo "All done!"

}
