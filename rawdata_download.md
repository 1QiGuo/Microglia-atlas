# Step1: Download FastQ from SRA database using nf-core_fetchngs

## Create a conda environment on OSC

```
conda create rawdatadownload
conda activate rawdatadownload
cd /fs/ess/PAS1475/guoqi/nf_core
```

## Install nf-core

### Install Nextflow on OSC

Nextflow enables scalable and reproducible scientific workflows using software containers.

```
#change java version on OSC to download nextflow
module load java/12.0.2
#download the executable package 
#(what is the bash mean?)
wget -qO- https://get.nextflow.io | bash
#make the binary executable on system
chmod +x nextflow
#Optionally, move the nextflow file to a directory accessible by your 
#$PATH variable (this is only required to avoid remembering and typing the 
#full path to nextflow each time you need to run it).
```
### Set the path of Nextflow to your your PATH environment variable.

```
export PATH=$PATH:$PWD/sratoolkit.3.0.0-mac64/bin
#check path
echo $PATH
```
Alternatively, we moved nextflow binary file to `/.conda/envs/rawdatadownload/bin/` directly.

## Install singularity on OSC
Singularity is an open source container platform designed to be simple, fast, and secure.
OSC can not use docker and singularity is an alternative to docker and do not need permission (sudo).
OSC has singularity so we don't need to install it.
```
singularity --version
#output(how to show the output?)
apptainer version 1.1.6
```

## Install nf-core using conda

```
conda install nf-core
```

## Configuretion (optional)
To let nextflow knows how to fetch the required software.
```
nf-core download fetchngs 1.9
```
If you are using singularity, please use the nf-core download command to download images first, before running the pipeline. Setting the NXF_SINGULARITY_CACHEDIR or singularity.cacheDir Nextflow options enables you to store and re-use the images from a central location for future pipeline runs. Images is like a compressed folder including all dependences and package.

$NXF_SINGULARITY_CACHEDIR and images have been added to `/users/PAS1475/guoqi/.bashrc`.
![image](https://user-images.githubusercontent.com/96271990/235548782-e9b48417-3dad-4c31-a143-f95e69e28d33.png)


## Download the pipeline and test it on a minimal dataset with a single command
change your profile to singularity, because we want to use singularity as container.
```
#original code
#nextflow run nf-core/fetchngs -profile test,YOURPROFILE --outdir <OUTDIR>
#revised code
nextflow run nf-core/fetchngs -profile test,singularity --outdir ./out
```
You finished the configuration for nf-core! Now, let's start to use on our own data!

## Start to run on Faith's data.
```
nextflow run nf-core/fetchngs --input SRR_Acc_List_faith.csv --outdir ./fetch_faith
```
## Bug
1. vim nextflow.config delete description
2. vim sra_ids_to_runinfo.py delete 1) accession 2) parent_study 3)description



```
