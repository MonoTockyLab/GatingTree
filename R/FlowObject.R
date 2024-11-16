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

#' @title SampleDef
#'
#' @description A class to represent a sample definition object
#'
#' @name SampleDef
#' @keywords classes
#' @docType class
setClass("SampleDef")


#' Create or update a SampleDef object in a FlowObject
#'
#' This function updates a FlowObject with a SampleDef object created from a CSV file or a provided data frame.
#'
#' @param x A FlowObject.
#' @param sample_df Optional data frame containing the sample definitions with columns 'file' and 'group'.
#' @param path A character string specifying the path to the directory containing the CSV file. Default is './sampledef'.
#' @param file A character string specifying the name of the CSV file. Default is 'sample.csv'.
#' @param group_column A character string specifying the column name for sample grouping in the CSV file.
#' @param confirm If TRUE, prompts confirmation of editing the CSV file.
#' @return A modified FlowObject containing a SampleDef object.
#' @rdname SampleDef
#' @aliases SampleDef
#' @export
#' @examples
#' \dontrun{
#' x <- SampleDef(x)
#' }
#' @family Initialization

SampleDef <- function(x, sample_df = NULL, path = 'sampledef', file = 'sample.csv', group_column = 'group', confirm = TRUE) {
    if (!inherits(x, "FlowObject")) {
        stop("Use a FlowObject for x.")
    }

    if (is.null(x@sampledef)) {
        x@sampledef <- list()
    }

    if (!is.null(sample_df)) {
        if (!all(c('file', 'group') %in% colnames(sample_df))) {
            stop("Data frame must contain 'file' and 'group' columns.")
        }
        x@sampledef$sampledef <- sample_df
    } else {
        if (confirm) {
            answer <- askYesNo("Have you edited sampledef/sample.csv?", default = FALSE)
            if (!answer) {
                stop("Please edit the sampledef/sample.csv file before proceeding.")
            }
        }
        scanlines <- scan(file = file.path(path, file), what = "character", sep = "\n", quiet = TRUE)
        tpsn <- strsplit(scanlines, split = ',')
        sn <- do.call(rbind, lapply(tpsn[-1], function(row) data.frame(file = row[1], group = row[2], stringsAsFactors = FALSE)))
        colnames(sn) <- c('file', 'group')
        
        if(!group_column %in% colnames(sn)) {
            stop(paste("Group column", group_column, "not found in the data."))
        }
        x@sampledef$sampledef <- sn
    }

    cat("Groups are: \n")
    cat(paste(levels(as.factor(x@sampledef$sampledef$group)), collapse = "\n"))
    
    return(x)
}


#' Show SampleDef
#' @param x A FlowObject.
#' @export
#' @examples
#' \dontrun{
#' showSampleDef(x)
#' }
showSampleDef <- function(x){
    if(!is(x, "FlowObject")){
        stop("Use a FlowObject for x. \n")
    }

    
    if(is.null(x@sampledef$sampledef)){
        cat("Apply SampleDef to your FlowObject. \n")
    }
    
    print(x@sampledef$sampledef)

}


#' A class representing a FlowObject object
#'
#' This class provides a representation of a FlowObject object.
#'
#' @keywords classes
#' @slot Data A data frame containing flow cytometric expression data
#' @slot sampledef A SampleDef object containing sample file definitions.
#' @slot metadata A list containing annotation data
#' @slot prep A list containing sample file definitions.
#' @slot Clustering A list containing cell identities and clustering data
#' @slot Gating A list containing gating information.
#' @slot QCdata A list containing quality control information.
#' @slot Transformation A list containing settings for data transformation, including the log data settings.
#' @importFrom methods setClass setMethod
#' @export

setClass("FlowObject", slots = list(
        Data = "data.frame",
        sampledef = "list",
        metadata = "list",
        prep = "list",
        Clustering = "list",
        Gating = "list",
        QCdata = "list",
        Transformation = "list"
    ),
validity = function(object){
           return(TRUE)
})





#' Initialize a FlowObject
#'
#' @param .Object A FlowObject
#' @param Data A data frame containing the raw flow cytometry data.
#' @param sampledef A SampleDef object containing sample file definitions.
#' @param metadata A list containing annotation data
#' @param prep A list of character string vectors defining the samples and controls.
#' @param Clustering A list containing cell identities and clustering data
#' @param Gating A list containing outputs of Gating.
#' @param QCdata A list containing the quality control information for the flow cytometry data.
#' @param Transformation A list containing settings for data transformation, including the log data settings.
#'
#' @return A FlowObject
#'
#' @export
#' @keywords internal

setMethod("initialize", "FlowObject",
  function(.Object,
           Data = data.frame(),
           sampledef = new("SampleDef"),
           QCdata = list(),
           Transformation = list(),
           prep = list(),
           Gating = list(),
           metadata = list(),
           Clustering = list()) {
    .Object@Data <- Data
    .Object@sampledef <- sampledef
    .Object@metadata <- metadata
    .Object@prep <- prep
    .Object@Clustering <- Clustering
    .Object@Gating <- Gating
    .Object@QCdata <- QCdata
    .Object@Transformation <- Transformation

    return(.Object)
  }
)




#' Show method for the FlowObject class
#'
#' Displays summary information for various slots of the FlowObject object
#'
#' @param object An object of the FlowObject class
#' @export
#' @method show FlowObject
#' @importFrom methods setClass setMethod
#' @keywords methods internal

show.FlowObject <- function(object) {
     cat("A FlowObject Class \n")
    cat(paste("Experiment: ", object@metadata$experiment, "\n"))
    cat(paste("Assay: ", nrow(object@Data), "cells and", ncol(object@Data), "variables \n"))
    
    if(length(object@sampledef) > 0 && !is.null(object@sampledef$sampledef)){
        l <- length(levels(as.factor(object@sampledef$sampledef$group)))
        n <- nrow(object@sampledef$sampledef)
        cat(paste("Experiment: ",l, "groups and ", n, "samples \n"))
    }
    cat("Available data: ")
    if(length(object@Data) > 0){
        cat(paste(colnames(object@Data)[1:5], collapse = ", "))
        cat('...\n')
    }

    if(length(object@Transformation$logdata_parameters) > 0){
        cat("LogData ")
    }

    if(length(object@QCdata$negative_gate_def) != 0){
        cat("DefineNegative ")
    }
    
    if(length(object@Gating$GatingTreeObject) != 0){
        cat("GatingTreeObject ")
    }
    
    cat("\n")

  }

setMethod("show", "FlowObject", show.FlowObject)



#' Create a new FlowObject
#'
#' This function constructs a new FlowObject using various inputs related to flow cytometry analysis,
#' including data, sample definitions, metadata, and data processing configurations.
#'
#' @param Data A data frame containing flow cytometric expression data.
#' @param sampledef An object of class SampleDef, specifying sample definitions.
#' @param metadata A list containing annotation data for the flow cytometry experiment.
#' @param prep A list of character string vectors defining the samples and controls.
#' @param Clustering A list containing cell identities and clustering data.
#' @param Gating A list containing outputs of gating, including GatingTree objects.
#' @param QCdata A list containing the quality control information for the flow cytometry data.
#' @param Transformation A list containing settings for data transformation, including the log data settings.
#' @return An object of class FlowObject, which encapsulates all provided data and settings for flow cytometry analysis.
#' @export
flowObject <- function(Data = data.frame(),
                        sampledef = new("SampleDef"),
                        metadata = list(),
                        prep = list(),
                        Clustering = list(),
                        Gating = list(),
                        QCdata = list(),
                        Transformation = list()
                ){
  

  if(!is.data.frame(Data)) {
    stop("Data must be a data frame")
  }
  if(!inherits(sampledef, "SampleDef")) {
    stop("sampledef must be an object of class SampleDef")
  }
  if(!is.list(metadata)) {
    stop("metadata must be a list")
  }

  if(!is.list(prep)) {
    stop("prep must be a list")
  }

  if(!is.list(Clustering)) {
    stop("Clustering must be a list")
  }

  if(!is.list(Gating)) {
    stop("Gating must be a list")
  }

  if(!is.list(QCdata)) {
    stop("QCdata must be a list")
  }

  if(!is.list(Transformation)) {
    stop("Transformation must be a list")
  }

  out <- new("FlowObject",
              Data = Data,
              metadata = metadata,
              sampledef = sampledef,
              prep = prep,
              Gating = Gating,
              Clustering = Clustering,
              QCdata = QCdata,
              Transformation = Transformation
            )
  
  return(out)
}
