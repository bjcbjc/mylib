library(Biobase, lib.loc="/ifs/home/c2b2/dp_lab/uda2001/.Rlibrary")
library(AnnotationDbi, lib.loc="/ifs/home/c2b2/dp_lab/uda2001/.Rlibrary")
library(affy);

args <- commandArgs(trailingOnly=TRUE);

if (length(args) < 2) {
    print("Usage: Rscript ./MAS5_ENTREZG.R Directory1[,Directory2] CDFName\n");
    stop();
}

currentDir <- getwd();

folders <- c(unlist(strsplit(args[1], ",")));
folders <- gsub("\"", "", folders);

#print(folders);
#print(args[2]);
#stop();

for (i in 1:length(folders)) {
  setwd(folders[i]);
  fn <- dir('.',pattern='.CEL')
  if (length(fn)) {
    print(paste('..........',folders[i],'..........',sep=''))
    for (j in 1:length(fn)) {
          if (length(dir(pattern=sub('.CEL','.ENTREZG.mas',fn[j])))) {
        print(paste('Skipping file', fn[j], 'since it was already processed.'));
      }
      else {
        data <- try(ReadAffy(filenames=fn[j], cdfname=args[2]))
        if(isTRUE(all.equal(class(data), "tryCatch"))) {
            writeLines("Model failed")
                              }
        print(.Last.value)
        print(paste('MAS5ing...', fn[j]), sep='')
        mas <- mas5(data)
        write.exprs(mas, file=sub('.CEL','.ENTREZG.mas',fn[j]));
        mascall <- mas5calls(data)
        write.exprs(mascall, file=sub('.CEL','.ENTREZG.call',fn[j]))
      }
    }
  }
  setwd(currentDir);
}
