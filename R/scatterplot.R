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

#' Plot Node Scatter Plot
#'
#' This function creates a scatter plot visualization of node percentages
#' based on the given path through a gating tree. It computes and displays
#' statistical summaries including means and standard deviations by group,
#' adds error bars to the plot, and marks significant differences with asterisks.
#'
#' @param x A TockyObject or FlowObject with GatingTreeObject
#' @param path A character vector specifying the path through the gating hierarchy.
#'
#' @return A ggplot object representing the node percentage scatter plot with
#'         error bars and optionally asterisks denoting statistical significance.
#'
#' @examples
#' \dontrun{
#' PlotNodeScatterPlot(x, c("path", "to", "node"))
#' }
#' @importFrom ggplot2 ggplot geom_jitter geom_errorbar labs theme_bw theme element_text
#' @importFrom stats wilcox.test
#' @importFrom gridExtra grid.arrange
#'
#' @export


PlotNodeScatterPlot <- function(x, path){
    
    node  <-  x@Gating$GatingTreeObject
    unique_maxima_df <- x@Gating$GatingTreeDF
    sampledef <- x@sampledef$sampledef
    
    rootdata <-  node$RootData
    total_cell_per_file <- table(rootdata$file)
    total_cell_per_file <- data.frame(file = names(total_cell_per_file), total = as.vector(total_cell_per_file))
    
    average_proportion <- p_value_vec <- rep(0, nrow(unique_maxima_df))
    percent_list <- as.list(p_value_vec)
    names(percent_list) <- unique_maxima_df$combination
    tmp <- findNodeByPath(node, path)
    node_data <- rootdata[rootdata$indices %in% tmp$CurrentNodeIndices,]
    total_cell_num <- table(node_data$file)
    total_cell_num <- data.frame(file = names(total_cell_num), node_total = as.vector(total_cell_num))
    total_cell_num <- merge(total_cell_num, data.frame( total_cell_per_file), by = 'file')
    total_cell_num <- merge(total_cell_num, sampledef, by = 'file')
    total_cell_num$node_percentage <- 100*total_cell_num$node_total/total_cell_num$total
    p_value = wilcox.test(data = total_cell_num, node_percentage~group)$p.value
    cat(show(p_value))
    
    df <- total_cell_num
    colnames(df)[4:5] <- c("Group", "Percentage")
    mean_df <- aggregate(Percentage ~ Group, data = df, FUN = function(x) mean(x, na.rm = TRUE))
    sd_df <- aggregate(Percentage ~ Group, data = df, FUN = function(x) sd(x, na.rm = TRUE))
    names(mean_df) <- c("Group", "Mean")
    names(sd_df) <-  c("Group", "SD")
    summary_df <- merge(mean_df, sd_df, by = "Group")
    ylab <- "Node Percentage"
    path <- paste(path, collapse = '_')
    title <- gsub("(.logdata)|[,]", "", gsub("[.](pos)", "+", path))
    title <- gsub("(.logdata)|[,]", "", gsub("[.](neg)", "-", title))
    title <- gsub("_", "", title)
    ylim_range <- c(0, NA)  # Adjust based on your data range
    
    gg <- ggplot() +
    geom_jitter(data = df, aes(x = Group, y = Percentage, colour = Group),
    position = position_jitter(width = 0.02)) +  # Corrected line
    geom_errorbar(data = summary_df, aes(x = Group, ymin = Mean - SD, ymax = Mean + SD, group = Group), width = 0.15, color = "grey") +
    geom_errorbar(data = summary_df, aes(x = Group, ymin = Mean, ymax = Mean, group = Group), width = 0.2, color = "grey") +
    labs(y = ylab, title = title) +
    ylim(ylim_range) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
    
    if (p_value < 0.05) {
        asterisk = ifelse(p_value < 0.001, "***", ifelse(p_value < 0.001, "**", "*"))
        y_position_asterisk <- max(summary_df$Mean + 0.5*summary_df$SD, na.rm = TRUE)
        
        gg <- gg + geom_text(aes(x = 1.5, y = y_position_asterisk, label = asterisk),
        vjust = 0,
        size = 6)
    }
    
    return(gg)
}

