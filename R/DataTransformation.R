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

#' Log fluorescence data
#' @param x A FlowObject containing non-logged, raw fluorescence data
#' @param graphics Whether to use a graphic device to choose variables for log transformation.
#' @param prompt Whether to invoke a prompt for asking optional questions.
#' @param variables A character vector for defining the variables to be log-transformed.
#' @return logged data. FlowObject will include the slot  logdata_parameters in the slot Transformation, which includes a list object of logged channel names and parameter s used for log transformation: x_log = log10 (x - s + 1)
#' @export
#' @examples
#' \dontrun{
#' x  <-  LogData(x) #x is a FlowObject
#'}
#' @family Data Transformation

LogData <- function(x, graphics = TRUE, variables = NULL, prompt = FALSE){
    if(inherits(x,'FlowObject')){
        
        data <- x@Data
    }
    
    if(is.data.frame(x)){
        data <- x
    }
    
    if(prompt){
        if(length(x@Transformation$logdata_parameters)!=0){
            ans <- askYesNo("Logged data found. Do you want to erase your previous logged data? \n")
            if(!ans){
                stop("Aborted. \n")
            }else{
                x@Transformation$logdata_parameters <- list()
                lg <- grepl(pattern= 'logdata', colnames(x@Data))
                x@Data <- x@Data[,!lg]
                data <- x@Data
            }
        }

        
    }else{
        x@Transformation$logdata_parameters <- list()
        lg <- grepl(pattern= 'logdata', colnames(x@Data))
        x@Data <- x@Data[,!lg]
        data <- x@Data
        
    }

    
    if(!inherits(x,'FlowObject')){stop("Input should be a FlowObject")
        }
        
        choices <- colnames(data)
        lg <- !(grepl(pattern = 'FSC',choices))& !(grepl(pattern = 'SSC',choices))&!(grepl(pattern = 'file',choices))&!(grepl(pattern = 'Angle',choices))&!(grepl(pattern = 'Intensity',choices))&!(grepl(pattern = 'Time',choices))&!(grepl(pattern = 'asinhdata',choices))
        
        choices <- choices[lg]
        
        if(is.null(variables)){
            selected <- select.list(choices, graphics = graphics, title = "Data to be logged", multiple =TRUE)
            
        }else{
            if(!all(variables %in% choices)){
                stop("Some variables were not found in the data. \n")
            }else{
                selected <- variables
            }
        }

        data_logged <- data[,selected]
        
        s <- rep(0, ncol(data_logged))
        data_logged <- as.matrix(data_logged)
        
        data_logged <- logTransformData(data_logged, 0.005, 0)
        
        if(prompt){
            ans <- askYesNo("Do you want to name your channel data? \n")
            if(ans){
                tmpname <- colnames(data_logged)
                outputname <- tmpname
                for(i in 1:length(tmpname)){
                    outputname[i] <- readline(paste("Name", tmpname[i], ": "))
                    cat("\n")
                }
                outputnamedf <- data.frame(original_variable = tmpname, custom_variable = outputname)
                x@Transformation$channel_names <- outputnamedf
            
            colnames(data_logged) <- paste(selected, 'logdata', sep='.')
            
            lg <- colnames(data_logged) %in% tmpname
            data_logged <- cbind(data[,!lg], data_logged)

            }else{
                colnames(data_logged) <- paste(selected, 'logdata', sep='.')
                lg <- colnames(data_logged) %in% selected
                data_logged <- cbind(data[,!lg], data_logged)

            }
        }else{
            colnames(data_logged) <- paste(selected, 'logdata', sep='.')
            lg <- colnames(data_logged) %in% selected
            data_logged <- cbind(data[,!lg], data_logged)
            
        }


           out <- x
           out@Data <- data_logged
           out@Transformation$logdata_parameters <- list(
                logged.channel = selected,
                s = s)
                
           return(out)

}





#' Moderate Extreme Negative Outliers with Random Values within Noise
#'
#' This function moderates extreme negative outliers in the specified variables
#' of a FlowObject by replacing these values with random numbers. These random
#' numbers are drawn from a normal distribution determined by the data lying
#' within the `negative_gate_def` thresholds in the FlowObject. It is typically
#' used to handle outliers in flow cytometry data after defining negative gates
#' using the DefineNegative function.
#'
#' @param x A FlowObject that has already been processed using DefineNegative.
#' @param var A character vector specifying the variables (markers) in the FlowObject
#'   for which the moderation of extreme negative values should be performed. If NULL,
#'   the user is prompted to select variables interactively.
#' @param output The output directory name for output files
#' @param plot Logical, whether to produce diagnotic plots.

#'
#' @return A modified FlowObject in which extreme negative values in the specified
#'   variables are replaced with random numbers based on the distribution of values
#'   within the defined negative gates. This modification is intended to reduce the
#'   impact of extreme outliers on subsequent analyses.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'x' is a valid FlowObject with required preprocessing:
#'   x <- NormAF(x, var = c("marker1", "marker2"))
#' }
#'
#' @export
#' @family Data Transformation


NormAF <- function(x, var = NULL, plot = FALSE) {

    if(!inherits(x, "FlowObject")){
        stop("Use a FlowObject for x.")
    }

    if (is.null(x@QCdata$negative_gate_def)) {
        stop("Perform DefineNegatives. \n")
    }
    
    neg_gate_def <- x@QCdata$negative_gate_def
    neg_gate_def <- as.data.frame(neg_gate_def)
    choices <- neg_gate_def$variable[neg_gate_def$negative.gate !=0]
    data <- x@Data[,choices]
    num_vars <- ncol(data)
    num_pairs <- num_vars / 2
    var_name <- sub(colnames(data), pattern = '.logdata', replacement = '')
    
    n <- round(nrow(neg_gate_def)^0.5) + 1
    
    if(plot){
        par(mfrow = c(n, n))
        if(num_vars >= 2){
            if(nrow(data) > 2000){
                sampling <- sample(1:nrow(data), 2000)
                data <- data[sampling,]
                
            }
            
        }else{
            if(length(data) > 2000){
                sampling <- sample(1:length(data), 2000)
                data <- data[sampling]
                
            }
        }
        
        for (i in seq(1, num_vars-1, by=2)) {
          x_var <- data[, i]
          y_var <- data[, i+1]
          
          plot(x_var, y_var, xlab = var_name[i], ylab = var_name[i+1], main = paste(var_name[i], "vs", var_name[i+1]), pch = '.')
        }

    }

    
    
    if(is.null(var)){
        choices <- neg_gate_def$variable[neg_gate_def$negative.gate !=0]
        var <- select.list(choices, graphics = TRUE, title = "Select variables for which autofluorescence should be corrected before UMAP", multiple =TRUE)

    }
    
    if(length(x@QCdata$rawData)==0){
        x@QCdata$rawData <- x@Data[,var]
    }else{
        tmp_cn <- setdiff(var, colnames(x@QCdata$rawData))
        x@QCdata$rawDat <- cbind(x@QCdata$rawDat, x@Data[, tmp_cn])
        
    }
    
    
    neg_gate_def <- neg_gate_def[neg_gate_def$variable %in% var,]

    
    if(plot){
        pdf(paste0(output, '/ModerateExtremeNegatives_thresholds.pdf'))
        par(mfrow = c(n,n))
        for (i in 1:nrow(neg_gate_def)) {
            channel_name <- neg_gate_def$variable[i]
            neg_threshold <- neg_gate_def$negative.gate[i]
            texpression <- x@Data[[channel_name]]
            out <- moderate_negative(texpression,neg_threshold, main = channel_name)
            x@Data[[channel_name]] <- out
        }
        data <- x@Data[,x@QCdata$negative_gate_def$variable]
        n <- round(ncol(data)^0.5)
        
        par(mfrow = c(n, n))
        if(length(var)>= 2){
            if(nrow(data) > 2000){
                sampling <- sample(1:nrow(data), 2000)
                data <- data[sampling,]
            }
        }else{
            if(length(data) > 2000){
                sampling <- sample(1:length(data), 2000)
                data <- data[sampling]
                
            }
        }

        num_vars <- ncol(data)
        num_pairs <- num_vars / 2
        var_name <- sub(colnames(data), pattern = '.logdata', replacement = '')
        for (i in seq(1, num_vars-1, by=2)) {
          x_var <- data[, i]
          y_var <- data[, i+1]
          
          plot(x_var, y_var, xlab = var_name[i], ylab = var_name[i+1], main=paste(var_name[i], "vs", var_name[i+1]), pch = '.')
          abline(v = neg_gate_def[var_name[i],2], h = neg_gate_def[var_name[i+1],2], col = 2)
        }
        dev.off()
    }else{
        for (i in 1:nrow(neg_gate_def)) {
            channel_name <- neg_gate_def$variable[i]
            neg_threshold <- neg_gate_def$negative.gate[i]
            texpression <- x@Data[[channel_name]]
            out <- moderate_negative(texpression,neg_threshold, main = channel_name)
            x@Data[[channel_name]] <- out
        }
    }


    return(x)
}


#' Interactively determine gates to define negative (and positive) regions for each marker expression
#'
#' This function allows interactive setting of thresholds on plots to define negative and
#' positive gates for markers. Users can choose to plot data as 2D scatter plots or as density plots,
#' and may utilize pseudocolor to enhance visualization.
#'
#' @param x A FlowObject containing the data.
#' @param select Logical, whether to interactively select markers to be processed.
#'               If FALSE, all logged expression variables will be used. Default is TRUE.
#' @param max_cells_displayed Maximum number of cells displayed in plots for determining a negative threshold.
#' @param y_axis_var The variable to use for the y-axis in the interactive plot, or 'Density' to use a density plot.
#'                   If NULL (default), the user will be prompted to select a variable from x@Data.
#' @param pseudocolour Logical, whether to use pseudocolor based on density in scatter plots.
#'                     Default is TRUE. The option FALSE will use monochrome.
#' @return A modified FlowObject containing updated threshold values for autofluorescence
#'         in the selected variables in the slot QCdata. You can repeat DefineNegatives to renew
#'         or adjust some or all of the negative thresholds for variables.
#' @export
#' @examples
#' \dontrun{
#' x <- DefineNegatives(x)
#' }
#' @family Data Transformation


DefineNegatives <- function(x, select = TRUE, max_cells_displayed = 30000, y_axis_var = NULL, pseudocolour = TRUE) {
    if (!inherits(x, "FlowObject")) {
        stop("Input must be a FlowObject.")
    }

    # Extract data and verify required components
    data <- x@Data
    if (is.null(x@Transformation$logdata_parameters$logged.channel)) {
        stop("No logged channels found in x@Transformation$logdata_parameters$logged.channel.")
    }
    logged_channels <- x@Transformation$logdata_parameters$logged.channel
    logged_channel_fullnames <- paste(logged_channels, 'logdata', sep = '.')

    # Initialize negative_gate_def if necessary
    if (is.null(x@QCdata$negative_gate_def)) {
        x@QCdata$negative_gate_def <- data.frame(
            variable = logged_channel_fullnames,
            negative.gate = apply(data[, logged_channel_fullnames, drop = FALSE], 2, min)
        )
    }

    # Marker selection
    if (!select) {
        var <- logged_channel_fullnames
    } else {
        selected_channels <- select.list(
            logged_channels,
            graphics = TRUE,
            title = "Select markers for defining negative gates",
            multiple = TRUE
        )
        if (length(selected_channels) == 0) {
            stop("No markers selected.")
        }
        var <- paste(selected_channels, 'logdata', sep = '.')
    }

    if (is.null(y_axis_var)) {
        choices <- colnames(data)[!grepl(colnames(data), pattern = 'file')]
        y_axis_var <- select.list(
            c(choices, 'Density'),
            graphics = TRUE,
            title = "Select variable for y-axis or choose 'Density' for density plot",
            multiple = FALSE
        )
        if (length(y_axis_var) == 0 || y_axis_var == "") {
            stop("No y-axis variable selected.")
        }
    } else {
        if (y_axis_var != 'Density' && !y_axis_var %in% colnames(data)) {
            stop(paste("y_axis_var", y_axis_var, "not found in data."))
        }
    }
    
    cat("Using y-axis variable:", y_axis_var, "\n")
    y.data <- data[[y_axis_var]]

    # Set sample size
    sample_size <- min(max_cells_displayed, nrow(data))

    selected_gates <- list()

    # Interactive gating for each marker
    for (channel in var) {
        if (!channel %in% colnames(data)) {
            warning("Channel", channel, "not found in data. Skipping.")
            next
        }

        # Sample data if necessary
        if (nrow(data) <= sample_size) {
            tmpdata <- data[[channel]]
            ty <- y.data
        } else {
            sample_indices <- sample(seq_len(nrow(data)), sample_size)
            tmpdata <- data[sample_indices, channel]
            ty <- y.data[sample_indices]
        }
        channel_original <- channel
        channel <- sub(channel, pattern = '.logdata', replacement = ' (log)')

        # Interactive gating loop
        repeat {
            if (y_axis_var == 'Density') {
                plot(density(tmpdata), main = paste("Density of", channel), xlab = channel, ylab = "Density")
            } else {
                ty <- if (nrow(data) <= sample_size) data[[y_axis_var]] else data[sample(1:nrow(data), sample_size), y_axis_var]
                
                if(!pseudocolour){
                    plot(x = tmpdata, y = ty, main = channel, pch = '.', col = 4, xlab = channel, ylab = y_axis_var)
                }else{
                    smoothScatter(x = tmpdata, y = ty, main = channel, xlab = channel, ylab = y_axis_var, colramp = colorRampPalette(c("white", rev(rainbow(10, end = 4/6)))))

                }
                


            }

            cat("Click on the plot to set the negative gate for", channel, "\n")
            scalegate <- locator(type = 'p', col = 2)

            if (is.null(scalegate)) {
                if (askYesNo("No point selected. Do you want to retry?")) {
                    next
                } else {
                    break
                }
            } else {
                abline(v = scalegate$x[length(scalegate$x)], col = 2)
                if (askYesNo("Happy with your gating?")) {
                    selected_gates[[channel_original]] <- scalegate$x[length(scalegate$x)]
                    break
                }
            }
        }

    }


    # After gates have been selected interactively...
    if (length(selected_gates) > 0) {
        # Create a data frame from selected gates
        final_gates <- data.frame(
            variable = names(selected_gates),
            negative.gate = unlist(selected_gates),
            stringsAsFactors = FALSE
        )

        # Update existing entries or add new ones
        existing_vars <- x@QCdata$negative_gate_def$variable
        for (i in seq_along(final_gates$variable)) {
            var_name <- final_gates$variable[i]
            if (var_name %in% existing_vars) {
                # Update the existing entry
                x@QCdata$negative_gate_def$negative.gate[existing_vars == var_name] <- final_gates$negative.gate[i]
            } else {
                # Add new entry
                x@QCdata$negative_gate_def <- rbind(x@QCdata$negative_gate_def, final_gates[i, ])
            }
        }
        # Reset the row names after updating
        rownames(x@QCdata$negative_gate_def) <- NULL
    } else {
        warning("No gates were defined.")
    }

    return(x)
}


#' Import Negative Gate Definitions into FlowObject
#'
#' This function imports a specified negative gate definition data frame into a FlowObject.
#' The data frame must contain exactly two columns: 'variable' and 'negative.gate'.
#'
#' @param x A FlowObject.
#' @param negative_gate_def A data frame containing the negative gate definitions.
#' @return The modified FlowObject with updated negative gate definitions.
#' @export
#' @examples
#' \dontrun{
#' x <- import_negative_gate_def(x, negative_gate_def = my_gate_definitions)
#' }
import_negative_gate_def <- function(x, negative_gate_def) {
    if (!inherits(x, "FlowObject")) {
        stop("Use a FlowObject for x.")
    }
    
    if (!is.data.frame(negative_gate_def) || length(colnames(negative_gate_def)) != 2 ||
        !all(c("variable", "negative.gate") %in% colnames(negative_gate_def))) {
        stop("negative_gate_def must be a data frame with columns 'variable' and 'negative.gate'.")
    }
    
    x@QCdata$negative_gate_def <- negative_gate_def
    return(x)
}

#' Export Negative Gate Definitions from FlowObject
#'
#' This function exports the negative gate definitions from a FlowObject to a CSV file.
#' It will stop if the negative gate definitions are not found within the FlowObject.
#'
#' @param x A FlowObject containing negative gate definitions.
#' @param filename The name of the file to which the negative gate definitions will be written.
#' @return The FlowObject, unchanged, with the side effect of writing to a file.
#' @export
#' @examples
#' \dontrun{
#' export_negative_gate_def(x, filename = "negative_gate_def.csv")
#' }
export_negative_gate_def <- function(x, filename = 'negative_gate_def.csv') {
    if (!inherits(x, "FlowObject")) {
        stop("Use a FlowObject for x.")
    }
    
    if (is.null(x@QCdata$negative_gate_def)) {
        stop("Negative gate definitions are not defined. Perform DefineNegatives first.")
    }
    
    df <- x@QCdata$negative_gate_def
    write.csv(df, filename, row.names = FALSE)
    return(x)
}

#' Plot the threshold for for each marker expression as per determined by DefineNegatives
#' @param x A FlowObject.
#' @param y_axis_var The variable to use for the y-axis in the generated plots, or 'Density' to use a density plot.
#' @param output If TRUE, plots are generated as a file. The default is FALSE and shows plots in X window.
#' @param panel The number of panels to be included. If NULL, all panels are included in the output plot.
#' @param outputFile When output is TRUE, this defines the filename of the output file.
#' @return The same FlowObject is returned for safety.
#' @export
#' @examples
#' \dontrun{
#' PlotDefineNegatives(x)
#'}
#' @family Data Transformation

PlotDefineNegatives <- function(x, y_axis_var = NULL, output = FALSE, outputFile = "DefineNegativePlot.pdf", panel = NULL){
    if(!inherits(x, "FlowObject")){
        stop("Use a FlowObject for x. \n")
    }
    
    if(length(x@QCdata$negative_gate_def) == 0){
        cat("Perform DefineNegatives before using this function. ")
    }

    
    data <- x@QCdata$negative_gate_def
    
    if (is.null(y_axis_var)) {
        choices <- colnames(x@Data)[!grepl(colnames(x@Data), pattern = 'file')]
        y_axis_var <- select.list(
            c(choices, 'Density'),
            graphics = TRUE,
            title = "Select variable for y-axis or choose 'Density' for density plot",
            multiple = FALSE
        )
        if (length(y_axis_var) == 0 || y_axis_var == "") {
            stop("No y-axis variable selected.")
        }
    } else {
        if (y_axis_var != 'Density' && !y_axis_var %in% colnames(x@Data)) {
            stop(paste("y_axis_var", y_axis_var, "not found in data."))
        }
    }
    
    
    if(output){
        if(grepl(".pdf", outputFile, ignore.case = TRUE)){
            pdf(file = outputFile)
        } else {
            jpeg(filename = outputFile, width = 480*3, height = 480*3)
        }
    }
    
    if(y_axis_var == 'Density'){
        if(is.null(panel)){
            par(mfrow = c(4, 4))
            for(i in 1:nrow(data)){
                variable_name <- sub(data$variable[i], pattern = '.logdata', replacement = ' (log)')
                tmpdata <- x@Data[, data$variable[i]]
                plot(density(tmpdata), main = variable_name, xlab = 'Expression', ylab = 'Density')
                abline(v = data$negative.gate[i], col = "red")
            }
        } else{
            par(mfrow = c(1, panel))
            for(i in 1:panel){
                variable_name <- sub(data$variable[i], pattern = '.logdata', replacement = ' (log)')
                tmpdata <- x@Data[, data$variable[i]]
                plot(density(tmpdata), main = variable_name, xlab = 'Expression', ylab = 'Density')
                abline(v = data$negative.gate[i], col = "red")
            }
        }

    } else {
        
        if(is.null(panel)){
            n <- ceiling(sqrt(nrow(data)))
            par(mfrow = c(n, n))
            for(i in 1:nrow(data)){
                variable_name <- sub(data$variable[i], pattern = '.logdata', replacement = ' (log)')
                tmpdata <- x@Data[, data$variable[i]]
                
                
                y_data <- x@Data[[y_axis_var]]

                smoothScatter(tmpdata, y_data, main = variable_name, xlab = variable_name, ylab = y_axis_var, colramp = colorRampPalette(c("white", rev(rainbow(10, end = 4/6)))))
                
                abline(v = data$negative.gate[i], col = "red")
            }
        }else{
            n <- ceiling(sqrt(panel))
            par(mfrow = c(1, panel))
            for(i in 1:panel){
                variable_name <- sub(data$variable[i], pattern = '.logdata', replacement = ' (log)')
                tmpdata <- x@Data[, data$variable[i]]
                y_data <- x@Data[[y_axis_var]]

                smoothScatter(tmpdata, y_data, main = variable_name, xlab = variable_name, ylab = y_axis_var, colramp = colorRampPalette(c("white", rev(rainbow(10, end = 4/6)))))
                
                abline(v = data$negative.gate[i], col = "red")
            }
        }
        

    }

    if(output){
        dev.off()
    }

    return(x)
    
}
    


#' To moderate extreme negative values in an expression vector
#' @param expression A numeric vector
#' @param neg A negative threshold value
#' @param main title of plot
#' @return A moderated expression data
#' @examples
#' \dontrun{
#' expression  <-  moderate_negative(x)
#'}
#'
#' @keywords internal

moderate_negative <- function(expression, neg, main){
    xlim = c(min(expression), max(expression))
    plot(density(expression), xlim = xlim, col = 8, main = main)
    abline(v = neg, col = 2)
    median_neg <- median(expression[expression < neg])
    abline(v = median_neg, col = 3)
    mean <- min(median_neg, neg)
    higher_neg_values <- expression[which(expression > median_neg & expression <= neg)]
    lower_neg_values <- expression[which(expression <= median_neg)]
    extreme_neg_indices <- which(expression < median_neg)
    
    sd_value <- sd(higher_neg_values)
    n_extreme_neg <- length(extreme_neg_indices)
    replacement_values <- rnorm(n_extreme_neg, mean, sd = sd_value)
    sorted_indices <- extreme_neg_indices[order(expression[extreme_neg_indices])]
    sorted_replacement_values <- sort(replacement_values)
    
    expression[sorted_indices] <- sorted_replacement_values
    
    par(new = T)
    plot(density(expression), xlim = xlim, col = 4, lty = 2, ann = FALSE, xaxt = 'n', yaxt = 'n')
    return(expression)
}









#' Plot Effects of Autofluorescnce Modelling
#'
#' This function moderates extreme negative outliers in the specified variables
#' of a FlowObject by replacing these values with random numbers. These random
#' numbers are drawn from a normal distribution determined by the data lying
#' within the `negative_gate_def` thresholds in the FlowObject. It is typically
#' used to handle outliers in flow cytometry data after defining negative gates
#' using the DefineNegative function.
#'
#' @param x A FlowObject that has already been processed using DefineNegative and NormAF.
#' @param graphics Logical. The default is FALSE, showing variable options in console.
#' @param output A character for the name of output directory.
#' @param max_cells_displayed The number of cells displayed in plots for determining a negative therehold.
#'
#' @return The same FlowObject is returned for safety.
#'
#' @examples
#' \dontrun{
#'   x <- PlotNormAF(x)
#' }
#'
#' @export
#' @family Data Transformation


PlotNormAF <- function(x, graphics = FALSE, max_cells_displayed = 30000, output = 'output') {

    if(!inherits(x, "FlowObject")){
        stop("Use a FlowObject for x.")
    }

    if (is.null(x@QCdata$negative_gate_def)) {
        stop("Perform DefineNegatives. \n")
    }
    
    if(!file.exists(output)){
        dir.create(output)
    }
    
    neg_gate_def <- x@QCdata$negative_gate_def
    neg_gate_def <- as.data.frame(neg_gate_def)
    choices_logdata <- neg_gate_def$variable[neg_gate_def$negative.gate !=0]
    choices <- sub(choices_logdata, pattern = '.logdata', replacement = '')
    
    variable1 <- select.list(choices, graphics = graphics, title = "The first variable to be plotted:", multiple = FALSE)
    variable2 <- select.list(choices, graphics = graphics, title = "The second variable to be plotted:", multiple = FALSE)
    
    variable_chosen <- c(variable1, variable2)
    
    variable_logged <- paste(variable_chosen, '.logdata', sep= '')

    rawdata <- as.matrix(x@QCdata$rawData[,variable_logged])
    raw_logged <- logTransformData(rawdata, 0.005, 0)
    af_logged <- x@Data[,variable_logged]
    
    
    if(nrow(raw_logged) > max_cells_displayed){
            sampling <- sample(1:nrow(raw_logged), max_cells_displayed)
            raw_logged_spl <- raw_logged[sampling,]
            af_logged_spl <- af_logged[sampling,]
            
        }else{
            raw_logged_spl <- raw_logged
            af_logged_spl <- af_logged
        }

    par(mfrow = c(2, 2))
    par(mar = c(4,3,1,1))
    par(cex.axis = 0.6)
    par(mgp = c(1.2, 0.2, 0))
    par(tck = -0.025)
    cex.lab = 1
    
    x_var <- raw_logged_spl[,1]
    y_var <- raw_logged_spl[,2]
    title <- paste('Log-Transformation:', variable_chosen[1], "vs", variable_chosen[2])
    plot(x_var, y_var, xlab = variable_chosen[1], ylab = variable_chosen[2], main = title, pch = '.', col = rgb(0,0,0,alpha = 0.2))
    abline(v = neg_gate_def[variable_chosen[1],2], h = neg_gate_def[variable_chosen[2],2], col = 2)

    x_var <- af_logged_spl[,1]
    y_var <- af_logged_spl[,2]
    title <- paste('AF Modelling:', variable_chosen[1], "vs", variable_chosen[2])
    plot(x_var, y_var, xlab = variable_chosen[1], ylab = variable_chosen[2], main = title, pch = '.', col = rgb(0,0,0,alpha = 0.2))
    abline(v = neg_gate_def[variable_chosen[1],2], h = neg_gate_def[variable_chosen[2],2], col = 2)

    neg_gate_def_chosen <- neg_gate_def[variable_logged,]

    for (i in 1:2) {
        channel_name <- neg_gate_def_chosen$variable[i]
        neg_threshold <- neg_gate_def_chosen$negative.gate[i]
        texpression <- raw_logged[,i]
        out <- moderate_negative(texpression,neg_threshold, main = channel_name)
    }

    return(x)
}
