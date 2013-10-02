#!/bin/bash
#$ -cwd
#$ -j y
#$ -o /ifs/data/c2b2/dp_lab/shares/Data/Human/CMAP2/mas5.output
#$ -l mem=3G,time=144:0:0

source /ifs/home/c2b2/dp_lab/bc2252/.bashrc

/nfs/apps/R/2.15.1/bin/Rscript /ifs/home/c2b2/dp_lab/bc2252/projects/rscripts/MAS5.compress.R /ifs/data/c2b2/dp_lab/shares/Data/Human/CMAP2/CEL1/