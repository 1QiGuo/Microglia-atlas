# Construct a submit file

naming: "submit.sh"

```
#!/bin/bash
#SBATCH --job-name=jianying
#SBATCH --time=10:50:59
#SBATCH --output=jianying_cellranger
#SBATCH --account=PCON0022
#SBATCH --mem=50GB
#SBATCH --mail-type=BEGIN,END,FAIL

mkdir log
while read NAME 
do
   sbatch --job-name=$NAME.run --output=./log/$NAME.out --export=NAME=$NAME cellranger_count.sh #need a shell file including speific cellranger function
   sleep 0.1s
done < /fs/scratch/PCON0022/Megan/Qi_microglia/fastq/debug_qi_temp/sample_unix.txt #need a txt file including sample name
```

# Construct a cellranger script

naming: cellranger_count.sh (can be changed but need to be consistent with the above code)
```
#!/usr/bin/bash
#SBATCH --account PCON0022
#SBATCH --time=8:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks=16
#SBATCH --mem=64GB

wd=/fs/scratch/PCON0022/Megan/Qi_microglia/fastq/debug_qi_temp
CellRanger=/fs/ess/PCON0022/tools/cellranger-7.1.0/cellranger
FastqFolder=/fs/scratch/PCON0022/Megan/Qi_microglia/fastq/raw_data
Refer=/fs/project/PCON0022/tools/refdata-gex-mm10-2020-A

#########################

cd $wd
echo $NAME
${CellRanger} count --id=$NAME --transcriptome=${Refer} --fastqs=${FastqFolder}/ --sample=$NAME --localcores=16 --localmem=64
```

# A sample name file

naming: sample.txt
It includes the names of Fastq files.

```
Y11339_WenH_WT_V1G_1
Y11340_WenH_KO_V1G_1
```

remove "\r" if there is an error.

```
tr -d '\r' < sample.txt > sample_unix.txt
```

