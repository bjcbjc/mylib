# This script processes one folder at a time, reading all files
# in that directory named *.CEL, and processing them one at a time
# using MAS5 algorithm.
# UNIX systems are usually case sensistive, so you might need to rename
# your files, and/or replace the pattern defined below.

# It can be run on the C2B2 server using this command
# /nfs/apps/R/2.9.0/bin/Rscript MAS5.R

library(affy);
# The directory which has the CEL files, change to the correct one
dirstr <- '/LabData/Human/DataFromPapers/CCLE/CCLE_Expression.Arrays_2012-04-05.CEL';
setwd(dirstr);

fn <- dir('.',pattern='.CEL')
  if (length(fn)) {
    for (j in 1:length(fn)) {
      # Skip already processed files. If you want to force reprocessing
      # remove the if statement
      if (length(dir(pattern=sub('.CEL','.mas',fn[j])))) {
        print(paste('Skipping file', fn[j], 'since it was already processed.'));
      }
      else {
        data <- ReadAffy(filenames=fn[j])
        print(paste('MAS5ing...', fn[j]), sep='')
        mas <- mas5(data)
        write.exprs(mas, file=sub('.CEL','.mas',fn[j]));
        # MAS5 calls are the Absent/Marginal/Present detection.
        # They are not very useful, so they are ignored.
        # If you need them for your data, uncomment these lines
        #mascall <- mas5calls(data)
        #write.exprs(mascall, file=sub('.CEL','.call',fn[j]))
    }
  }
 setwd('..')
}

