# Preprocessing downloaded fastq files

## Download preprocessing pipeline

The pipeline 'scrnaseq' from nfcore will be applied to preprocess fastq files from 10X.

Pull the code of 'scrnaseq' from Github

```
cd personal directory
git clone https://github.com/nf-core/scrnaseq.git
```

## Debug the pipeline 'scrnaseq'

1. We found errors when scrnaseq pull images from the docker hub. In docker hub, cellranger is 'nfcore/cellranger' instead of 'nf-core/cellranger'.
Hence, we think there are typos in 'scrnaseq' pipeline. We revised code in three files named "main.nf".
<img width="1170" alt="image" src="https://github.com/1QiGuo/Microglia-atlas/assets/96271990/7550f9e9-41d8-4a54-ab8f-d7a58e974c49">

```
cd /fs/ess/PAS1475/guoqi/nf_core/scrnaseq/modules/nf-core/cellranger
#original code
container "nf-core/cellranger:7.1.0"
#revised code
container "nfcore/cellranger:7.1.0"
```

2. Rename fastq files
We found cellranger has its own compatible name for each fastq file. Hence, we change file names in a batch.
```
#check whether the right files are chosen
for f in S*; do echo $f; done
#
for f in SRX*; do srr_name="${f#*_}" srr_name="${srr_name%%_*}" new_name=; done
```
