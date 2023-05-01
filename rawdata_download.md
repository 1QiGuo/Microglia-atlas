# Step1: Download FastQ from SRA database using nf-core_fetchngs

## Create a conda environment on OSC

```
conda create rawdatadownload
conda activate rawdatadownload
cd /fs/ess/PAS1475/guoqi/nf_core
```

## Install nf-core

### Install Nextflow on OSC

```
#download the executable package 
#(what is the bash mean?)
wget -qO- https://get.nextflow.io | bash
#make the binary executable on system
chmod +x nextflow
#Optionally, move the nextflow file to a directory accessible by your 
#$PATH variable (this is only required to avoid remembering and typing the 
#full path to nextflow each time you need to run it).
```

## Install singularity on OSC
Singularity is an open source container platform designed to be simple, fast, and secure.
OSC can not use docker and singularity is an alternative to docker.
OSC has singularity so we don't need to install it.
```
ingularity --version
#output(how to show the output?)
apptainer version 1.1.6
```
