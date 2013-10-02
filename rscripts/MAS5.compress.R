library(affy);
library(AnnotationDbi);

args <- commandArgs(trailingOnly=TRUE);

if (length(args) < 1) {
    print("Usage: Rscript ./MAS5_ENTREZG.R Directory1[,Directory2]\n");
    stop();
}

currentDir <- getwd();

folders <- c(unlist(strsplit(args[1], ",")));
for (i in 1:length(folders)) {
  setwd(folders[i]);
  fn <- dir('.',pattern='.CEL.gz$')
  if (length(fn)) {
    print(paste('..........',folders[i],'..........',sep=''))
    for (j in 1:length(fn)) {
      if (length(dir(pattern=sub('.CEL.gz','.mas',fn[j])))) {
        print(paste('Skipping file', fn[j], 'since it was already processed.'));
      }
      else {
        data <- ReadAffy(filenames=fn[j], compress=TRUE)
        print(paste('MAS5ing...', fn[j]), sep='')
        mas <- mas5(data)
        write.exprs(mas, file=sub('.CEL.gz','.mas',fn[j]));
        mascall <- mas5calls(data)
        write.exprs(mascall, file=sub('.CEL.gz','.call',fn[j]))
      }
    }
  }
  setwd(currentDir);
}
