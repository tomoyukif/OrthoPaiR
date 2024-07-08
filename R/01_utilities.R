#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5delete h5write
# Function to overwrite an HDF5 dataset
.h5overwrite <- function(obj, file, name){
    # Open the HDF5 file
    file <- H5Fopen(file)

    # Ensure the file is closed when the function exits
    on.exit(H5Fclose(file))

    # Check if the specified name exists in the file
    if(H5Lexists(h5loc = file, name = name)){
        # Delete the existing dataset if it exists
        h5delete(file = file, name = name)
    }

    # Write the new object to the HDF5 file
    h5write(obj = obj, file = file, name = name)
}

#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists h5createGroup
# Function to create a group in an HDF5 file
.h5creategroup <- function(file, name){
    # Open the HDF5 file
    file <- H5Fopen(file)

    # Ensure the file is closed when the function exits
    on.exit(H5Fclose(file))

    # Check if the specified group name does not exist in the file
    if(!H5Lexists(h5loc = file, name = name)){
        # Create a new group in the HDF5 file
        h5createGroup(file, name)
    }
}

#' @importFrom rtracklayer import.gff
# Function to import multiple GFF files
.importAllGFF <- function(fn){
    for(i in seq_along(fn)){
        if(i == 1){
            # Import the first GFF file
            out <- import.gff(fn[i])

        } else {
            # Concatenate subsequent GFF files
            out <- c(out, import.gff(fn[i]))
        }
    }
    # Return the concatenated GFF data
    return(out)
}

#' @importFrom rhdf5 h5createFile
# Function to create an HDF5 file
.makeHDF5 <- function(hdf5_path, overwrite){
    if(file.exists(hdf5_path)){
        # Notify if the file already exists
        message(hdf5_path, " already exists.")

        if(overwrite){
            message("Since user specified 'overwrite = TRUE', the exsiting HDF5 will be overwritten.")

            unlink(x = hdf5_path)

            # Create a new HDF5 file
            h5createFile(file = hdf5_path)
        }

    } else {
        # Create a new HDF5 file
        h5createFile(file = hdf5_path)
    }

    # Return the HDF5 file path
    return(hdf5_path)
}
