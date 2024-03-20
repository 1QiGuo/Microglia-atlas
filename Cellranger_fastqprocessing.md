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
done < /fs/ess/PCON0022/Jianyingli/yard/sample.txt #need a txt file including sample name
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

wd=/fs/ess/PCON0022/Jianyingli/yard
CellRanger=/fs/ess/PCON0022/tools/cellranger-7.1.0/cellranger
FastqFolder=$wd
Refer=/fs/project/PCON0022/tools/refdata-gex-mm10-2020-A

#########################

cd $wd
echo $NAME
${CellRanger} count --id=$NAME --transcriptome=${Refer} --fastqs=${FastqFolder}/$NAME --sample=$NAME --localcores=16 --localmem=64
```

# A sample name file
naming: sample.txt
It include the names of fastq files.

```
Y11339_WenH_WT_V1G_1
Y11340_WenH_KO_V1G_1
```

# debug
```
tr -d '\r' < sample.txt > sample_unix.txt
```
