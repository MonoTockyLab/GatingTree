# Copyright 2024 Masahiro Ono
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create FlowObject by importing flow cytometric data
#'
#' This function constructs a new FlowObject using either data imported from CSV files or data provided directly as arguments.
#'
#' @param path Path to CSV files.
#' @param select Whether to select files via an interactive window.
#' @param sample_file When select is FALSE, a vector of sample file names.
#' @param Data A data frame containing flow cytometric expression data. If provided, file reading is skipped.
#' @param sampledef An object of class SampleDef or a list containing sample definitions. If provided, sampledef population is skipped.
#' @param experiment_name Name of the experiment.
#' @return A FlowObject
#' @export
#' @examples
#' \dontrun{
#' x  <-  CreateFlowObject()
#'}
#' @import utils graphics grDevices methods stats
#' @family Initialization
CreateFlowObject <- function(path = '.', select = TRUE, sample_file = NULL, Data = NULL, sampledef = NULL, experiment_name = NULL) {
      if (is.null(experiment_name) || experiment_name == "") {
        experiment_name <- readline("Enter Experiment Name: ")
        if (is.null(experiment_name) || experiment_name == "") {
            experiment_name <- "Flow Experiment"
        }
    }
      
      check_variables_consistency <- function(files) {
        if (length(files) == 0) {
          stop("No files provided.")
        }
        
        # Initialize variables
        common_vars <- NULL
        inconsistent_files <- character(0)
        all_vars_list <- list()
        
        # Read column names from each file
        for (i in seq_along(files)) {
          vars <- tryCatch({
            colnames(read.table(files[i], nrows = 2, sep = ',', header = TRUE))
          }, error = function(e) {
            stop(paste("Error reading file:", files[i], "\n", e$message))
          })
          
          all_vars_list[[files[i]]] <- vars
          
          if (is.null(common_vars)) {
            common_vars <- vars
          } else {
            common_vars <- intersect(common_vars, vars)
          }
        }
        
        # Identify files with inconsistent columns
        for (file in files) {
          vars <- all_vars_list[[file]]
          if (!setequal(vars, common_vars)) {
            inconsistent_files <- c(inconsistent_files, basename(file))
          }
        }
        
        # Return list with common variables and inconsistent files
        return(list(common_vars = common_vars, inconsistent_files = inconsistent_files))
      }

      
      
      # Initialize variables
      metadata <- list(
        experiment = experiment_name,
        analysis_date = Sys.Date(),
        R_version = R.version$version.string
      )

      if (!is.null(Data)) {
         # Use provided Data
         if (!is.data.frame(Data)) {
           stop("Data must be a data frame.")
         }
         # Ensure 'file' column exists in Data
         if (!'file' %in% colnames(Data)) {
           stop("Data must contain a 'file' column indicating sample file names.")
         }
         commonvariables <- colnames(Data)
         sample_files <- unique(Data$file)
         prep <- list(samplefile = sample_files, variables = commonvariables)
       } else {
         # Read data from files
         if (select) {
           list_files <- list.files(path, pattern = 'csv', full.names = TRUE)
           samplefile <- select.list(list_files, graphics = TRUE, title = "Choose sample CSV files to be analyzed", multiple = TRUE)
           if (length(samplefile) == 0) {
             stop("No files selected.")
           }
         } else {
           if (is.null(sample_file) || length(sample_file) == 0) {
             stop("No sample files provided.")
           }
           samplefile <- sample_file
         }
         
         # Check variables consistency
         result <- check_variables_consistency(samplefile)
         commonvariables <- result$common_vars
         inconsistent_files <- result$inconsistent_files
         
         if (length(commonvariables) == 0) {
           stop("No common variables found across files.")
         }
         
         if (length(inconsistent_files) > 0) {
           warning("Some sample files do not have consistent column names:\n", paste(inconsistent_files, collapse = ", "))
         }
         
         # Read and combine data
         Data <- data.frame()
         cellcount.total <- numeric(length(samplefile))
         names(cellcount.total) <- basename(samplefile)
         
         for (i in seq_along(samplefile)) {
           tmpdata <- tryCatch({
             read.table(samplefile[i], sep = ',', header = TRUE)
           }, error = function(e) {
             stop(paste("Error reading file:", samplefile[i], "\n", e$message))
           })
           
           # Retain only common variables
           tmpdata <- tmpdata[, commonvariables, drop = FALSE]
           
           cellcount.total[i] <- nrow(tmpdata)
           tmpdata$file <- basename(samplefile[i])
           Data <- rbind(Data, tmpdata)
           cat(".")
         }
         cat("\n")
         
         # Create prep list
         prep <- list(samplefile = samplefile, variables = commonvariables)
       }

    
    Clustering <- list(original_cell_id = seq_len(nrow(Data)))
    Gating <- list(ancestor = 'top')
    Transformation <- list(logdata_parameters = list())
    QCdata <- list(cellcount.total = as.data.frame(table(Data$file)))


    # Handle sampledef
     if (is.null(sampledef)) {
       # Create sampledef directory and file
       output_sampledef <- file.path(path, "sampledef")
       if (!dir.exists(output_sampledef)) {
         dir.create(output_sampledef)
       }
       sample_def_file <- file.path(output_sampledef, "sample.csv")
       write.table(
         data.frame(file = unique(Data$file), group = ""),
         file = sample_def_file,
         quote = FALSE,
         row.names = FALSE,
         sep = ','
       )
       message("Please edit the sample definitions in ", sample_def_file)
       
       # Initialize sampledef as empty list; user will populate it later
       sampledef <- list()
     } else {
       # Use provided sampledef
       if (!inherits(sampledef, "SampleDef") && !is.list(sampledef)) {
         stop("sampledef must be an object of class SampleDef or a list.")
       }
       sampledef <- sampledef
       sampledef[['sampledef']] <- sampledef
     }

 out <- new("FlowObject",
              Data = Data,
              metadata = metadata,
              sampledef = sampledef,
              prep = prep,
              Clustering = Clustering,
              Gating = Gating,
              QCdata = QCdata,
              Transformation = Transformation
   )
   
    return(out)
}

