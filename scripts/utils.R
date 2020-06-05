

table_writer_checker <- function(fn, data){
  towrite<-NULL
  if (file.exists(fn)){
    repeat {
      towrite <- readline(prompt=paste("File", fn, "exists. Do you want to overwrite? (y/n)"))
      if (towrite=="n"){
        print("File was not overwritten")
        break
      } else if (towrite=="y"){
        print(paste("File",fn, "rewritten"))
        write.table(data, fn, sep = "\t", quote = F, row.names = F)
        towrite<-"n"
        break
      } else {
        print("Input not recognized. Please enter either y or n")
      }
    }
  } else{
    write.table(data, fn, sep = "\t", quote = F, row.names = F)
  }
}

get_preprocessing_fn <- function(genome, gene, track) {
  paste(paste("annotated", get_genome_nm(genome), gene, track, sep="_"), ".txt", sep="")
}

get_genome_nm <- function(genome) {
  paste(strsplit(genome, "\\.")[[1]][3], "_", strsplit(genome, "\\.")[[1]][4], sep="")
}

move_download <-function(gene, animal, genome, downloads_path, fn, dir_nm=NULL) {
  fp<-paste(downloads_path, fn, sep="/") #file path
  if (dir_nm==NULL) {
    dir_nm<-paste(gene, animal, get_genome_nm(genome), "input_data", sep="/") #directory path
  }
  if (file.exists(paste(dir_nm, fn, sep="/"))){
    print(paste("File already exists in", dir_nm, "repository. No files were moved or deleted"))
  }
  else{
    if (file.copy(fp, dir_nm)){
      print(paste(fp, "file moved to ", dir_nm, "!", sep=""))
      hold<-file.remove(fp) #assigning hold prevents unnecessary output
    }
    else{
      print(paste("Error: Either file", fp, "or directory", dir_nm, "does not exist. Please check for both before continuing."))
    }
  }
}