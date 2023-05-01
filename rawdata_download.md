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
$PATH variable (this is only required to avoid remembering and typing the full path to nextflow each time you need to run it).
```
