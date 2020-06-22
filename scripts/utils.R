table_writer_checker <- function(fp, data){
  #Helper function to check if there already exists a file at fp - asks if user wants to overwrite if table does exist
  #fp <- string with the file path name
  #data <- data to be written
  towrite<-NULL
  if (file.exists(fp)){
    repeat {
      towrite <- readline(prompt=paste("File", fp, "exists. Do you want to overwrite? (y/n)"))
      if (towrite=="n"){
        print("File was not overwritten")
        break
      } else if (towrite=="y"){
        print(paste("File",fp, "rewritten"))
        write.table(data, fp, sep = "\t", quote = F, row.names = F)
        towrite<-"n"
        break
      } else {
        print("Input not recognized. Please enter either y or n")
        print(towrite)
      }
    }
  } else{
    write.table(data, fp, sep = "\t", quote = F, row.names = F)
  }
}

read_table_helper <- function(fp, header=F, sep="") {
  #Helper function to check if table exists - produces readable error message if table isn't found
  #fp <- string with the file path name
  #header, sep: same as read.table input
  tryCatch({read.table(fp, header = header, stringsAsFactors = F, sep=sep)
  }, error = function(e) {
    stop(paste("Expected ", fp, " but did not find it. Is the file named correctly?", sep=""))
  })
}

get_preprocessing_fn <- function(genome, gene, track) {
  #Combines genome, gene, and track names together to produce file name for annotations to be saved to from UCSC table browser
  #genome <- name of the genome package
  #gene <- name of the gene
  #track <- name of the track used to acquire the features
  paste(paste("annotated", get_genome_nm(genome), gene, track, sep="_"), ".txt", sep="")
}

get_genome_nm <- function(genome) {
  #Helper function that gets the genome name from package name by removing details about the package and species.
  #genome <- name of the genome package
  paste(strsplit(genome, "\\.")[[1]][3], "_", strsplit(genome, "\\.")[[1]][4], sep="")
}

move_download <-function(gene, species, genome, downloads_path, fn, dir_nm=NULL) {
  #Helper function that moves file donwloaded with name according to get_preprocessing_fn from the downloads repository to correct directory
  #downloads_path <- full path to folder where the file was downloaded
  #dir_nm <- (optional) if path doesn't match convention, specifies where the file should be moved to
  fp<-paste(downloads_path, fn, sep="/") #file path
  if (is.null(dir_nm)) {
    dir_nm<-paste("data", gene, species, get_genome_nm(genome), "input_data", sep="/") #directory path
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