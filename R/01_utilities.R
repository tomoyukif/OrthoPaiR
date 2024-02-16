#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5delete h5write
.h5overwrite <- function(obj, file, name){
    file <- H5Fopen(file)
    on.exit(H5Fclose(file))
    if(H5Lexists(h5loc = file, name = name)){
        h5delete(file = file, name = name)
    }
    h5write(obj = obj, file = file, name = name)
}

#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5createGroup
.h5creategroup <- function(file, name){
    file <- H5Fopen(file)
    on.exit(H5Fclose(file))
    if(!H5Lexists(h5loc = file, name = name)){
        h5createGroup(file, name)
    }
}

#' @importFrom rtracklayer import.gff
.importAllGFF <- function(fn){
    for(i in seq_along(fn)){
        if(i == 1){
            out <- import.gff(fn[i])
        } else {
            out <- c(out, import.gff(fn[i]))
        }
    }
    return(out)
}

#' @importFrom rhdf5 h5createFile
.makeHDF5 <- function(hdf5_path){
    if(file.exists(hdf5_path)){
        message(hdf5_path, " already exists.")
    } else {
        h5createFile(file = hdf5_path)
    }
    return(hdf5_path)
}
