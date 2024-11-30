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

#' Two-Dimensional Plot with Gating and Quadrant Statistics
#'
#' This function plots a two-dimensional flow cytometry data highlighting the gated regions
#' and calculates quadrant statistics. It handles data preprocessing to manage outliers by
#' replacing extreme negative values with values from a normal distribution within specified
#' gates. The function is ideal for analyzing and visualizing flow cytometry data, providing
#' insights into the distribution of cell populations across specified markers.
#'
#' @param x A FlowObject that has already been processed using DefineNegative and NormAF.
#' @param graphics Logical, if TRUE, enables graphical selection of markers and gating thresholds via a GUI; otherwise, selections must be input manually. Defaults to FALSE.
#' @param max_cells_displayed The maximum number of cells to display in the plots, which can help manage performance and clarity in visualizing dense datasets.
#' @param output The directory path where the output plots and data summaries will be saved.
#' @param markers Optionally, a vector of marker names to be used for gating; if NULL, the function prompts for selection.
#' @param states A vector indicating the gating state ('positive' or 'negative') for each marker; if NULL, the function prompts for selection.
#' @param gating Logical, if TRUE, applies gating based on the markers and states provided; defaults to TRUE.
#' @param split_group Logical, if TRUE, splits the data by 'group' variable within the dataset for separate analysis and plotting; defaults to TRUE.
#'
#' @return Returns the same FlowObject for safety.
#'
#' @examples
#' \dontrun{
#'   plotFlow2D(x)
#' }
#'
#' @export

plotFlow2D <- function(x, graphics = FALSE, output = 'output', markers = NULL, states = NULL, gating = TRUE, split_group = TRUE, max_cells_displayed = 30000) {
    
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
    rownames(neg_gate_def) <- x@QCdata$negative_gate_def$variable
    neg_gate_def <- as.data.frame(neg_gate_def)
    choices_logdata <- neg_gate_def$variable[neg_gate_def$negative.gate !=0]
    choices <- sub(choices_logdata, pattern = '.logdata', replacement = '')
    X <- x@Data
    X <- X[, c(choices_logdata, 'file')]
    colnames(X) <- sub(colnames(X), pattern = '.logdata', replacement = '')
    
    X <- merge(X, x@sampledef$sampledef, by = 'file')
    
    if(gating) {
        if(is.null(markers) || length(markers) == 0) {
            markers <- select.list(choices, graphics = graphics, title = "Markers to be used for Gating:\n", multiple = TRUE)
            if(is.null(markers)) { # In case user cancels the selection
                stop("No markers selected. Exiting.")
            }
        }
        
        if(!all(markers %in% choices)) {
            stop("Selected markers must be within provided choices.")
        }
        if(is.null(states) || length(states) != length(markers)) {
            states <- vector("character", length(markers))
            for(i in 1:length(markers)) {
                ttitle <- paste("Choose whether positive or negative cells are gated for", markers[i], "\n")
                states[i] <- select.list(c("positive", "negative"), graphics = graphics, title = ttitle, multiple = FALSE)
            }
        } else if(length(states) != length(markers)) {
            stop("The length of states and markers must be the same.")
        } else  if(!all(states %in% c("positive", "negative"))) {
            if(all(states %in% c("Positive", "Negative", "positive", "negative"))){
                states <- tolower(states)
            } else {
                stop("States are either Positive or Negative.")
            }
            
        }
    }
    
    
    if(gating){
        logic <- rep(TRUE, nrow(X))
        for(i in 1:length(states)){
            if(states[i]=='positive'){
                logic <- logic & X[,markers[i]] > neg_gate_def[markers[i], "negative.gate"]
                
            }else{
                logic<- logic & X[,markers[i]] <= neg_gate_def[markers[i], "negative.gate"]
            }
        }
        X <- X[logic,]
    }
    show(summary(logic))
    show(dim(X))
    
    variable1 <- select.list(choices, graphics = graphics, title = "The first variable to be plotted (x-axis):", multiple = FALSE)
    variable2 <- select.list(choices, graphics = graphics, title = "The second variable to be plotted (y-axis):", multiple = FALSE)
    
    variable_chosen <- c(variable1, variable2)
    variable_logged <- paste(variable_chosen, '.logdata', sep= '')
    
    par(mfrow = c(2, 2))
    par(mar = c(4,3,1,1))
    par(cex.axis = 0.6)
    par(mgp = c(1.2, 0.2, 0))
    par(tck = -0.025)
    cex.lab = 1
    
    xrange <- range(X[,variable1])
    yrange <- range(X[,variable2])
    stateSymbols <- ifelse(states == "positive", "+", "-")
    
    if(split_group){
        X$group <- as.factor(X$group)
        
        X_list <- split(X, f = X$group)
        
        for(i in 1:length(X_list)){
            tmpX <- X_list[[i]]
            show(dim(tmpX))
            
            x_var <- tmpX[,variable1]
            y_var <- tmpX[,variable2]
            
            vline <- neg_gate_def[variable_logged[1], 2]
            hline <- neg_gate_def[variable_logged[2], 2]
            
            quadrant <- ifelse(y_var > hline, ifelse(x_var > vline, "UR", "UL"),
            ifelse(x_var > vline, "LR", "LL"))
            
            tdf <- data.frame(x_var = x_var, y_var = y_var, file = tmpX$file, quadrant = quadrant)
            quad_stats <- aggregate(cbind(count = file) ~ file + quadrant, data = tdf, FUN = length)
            quad_stats$total <- ave(quad_stats$count, quad_stats$file, FUN = sum)
            quad_stats$percent <- (quad_stats$count / quad_stats$total) * 100
            list_by_quadrant <- split(quad_stats$percent, quad_stats$quadrant)
            summary_stats <- lapply(list_by_quadrant, function(x) {
              c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
            })

            quad_summary <- do.call(rbind, summary_stats)
            quad_summary <- as.data.frame(quad_summary)
            quad_summary$quadrant <- rownames(quad_summary)
            rownames(quad_summary) <- NULL
            
            cat(paste(names(X_list)[i], ':\n'))

            quad_table <- table(quadrant)
            
            if(gating){
                title <- paste('Gate:', paste(paste0(markers, stateSymbols), collapse = ' '), '; Group:', names(X_list)[i])
            }else{
                title <- paste('Group:', names(X_list)[i])
            }
            if(!is.null(max_cells_displayed)){
                if(length(x_var)> max_cells_displayed){
                    sampling <- sample(1:length(x_var), max_cells_displayed)
                }else{
                    sampling <- 1:length(x_var)
                }
            }else{
                sampling <- 1:length(x_var)
            }
            plot(x_var[sampling], y_var[sampling], xlab = variable_chosen[1], ylab = variable_chosen[2], main = title, pch = '.', col = rgb(0,0,0,alpha = 0.2),
            xlim = xrange, ylim = yrange)
            abline(v = vline, h = hline, col = 2)
            
            x_right_pos <- 0.6*(max(x_var) - min(x_var))+min(x_var)
            x_left_pos <-min(x_var)
            y_low_pos <- 0.3*(hline - min(y_var))+min(y_var)

              
            text(x = c(x_left_pos, x_right_pos, x_left_pos, x_right_pos),
            y = c(y_low_pos, y_low_pos, max(y_var) * 0.95, max(y_var) * 0.95),
            labels = sprintf("%.2f%%\n(sd:%.2f%%)", quad_summary$mean, quad_summary$sd),
            
            pos = c(4, 4, 4, 4), cex = 1, col = 'darkblue')
        }
        
    }else{
        x_var <- X[,variable1]
        y_var <- X[,variable2]
        
        # Thresholds for gates
        vline <- neg_gate_def[variable_logged[1], 2]
        hline <- neg_gate_def[variable_logged[2], 2]
        
        # Classification of cells into quadrants
        quadrant <- ifelse(y_var > hline, ifelse(x_var > vline, "UR", "UL"),
        ifelse(x_var > vline, "LR", "LL"))
        
        tdf <- data.frame(x_var = x_var, y_var = y_var, file = X$file, quadrant = quadrant)
        quad_stats <- aggregate(cbind(count = file) ~ file + quadrant, data = tdf, FUN = length)
        quad_stats$total <- ave(quad_stats$count, quad_stats$file, FUN = sum)
        quad_stats$percent <- (quad_stats$count / quad_stats$total) * 100
        list_by_quadrant <- split(quad_stats$percent, quad_stats$quadrant)
        summary_stats <- lapply(list_by_quadrant, function(x) {
          c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
        })

        quad_summary <- do.call(rbind, summary_stats)
        quad_summary <- as.data.frame(quad_summary)
        quad_summary$quadrant <- rownames(quad_summary)
        rownames(quad_summary) <- NULL

        cat('Quadrant Percentages: \n')
        show(quad_stats)
        if(gating){
            title <- paste('Gate:', paste(paste0(markers, stateSymbols), collapse = ' '))
        }else{
            title <- paste('All Cells')
        }
        if(!is.null(max_cells_displayed)){
            if(length(x_var)> max_cells_displayed){
                sampling <- sample(1:length(x_var), max_cells_displayed)
            }else{
                sampling <-1:length(x_var)
            }
        }else{
            sampling <-1:length(x_var)
        }
        plot(x_var[sampling], y_var[sampling], xlab = variable_chosen[1], ylab = variable_chosen[2], main = title, pch = '.', col = rgb(0,0,0,alpha = 0.2))
        abline(v = vline, h = hline, col = 2)
        
        x_right_pos <- 0.7*(max(x_var) - min(x_var))+min(x_var)
        x_left_pos <- 0.3*(max(x_var) - min(x_var))+min(x_var)
        y_low_pos <- 0.3*(hline - min(y_var))+min(y_var)
        
        text(x = c(x_left_pos, x_right_pos, x_left_pos, x_right_pos),
        y = c(max(y_var) * 0.95, max(y_var) * 0.95, y_low_pos, y_low_pos),
        labels = sprintf("%.2f%%\n(sd:%.2f%%)", quad_summary$mean, quad_summary$sd),
        
        pos = c(2, 4, 2, 4), cex = 1, col = 'darkblue')
    }

    
    return(invisible(x))
}






