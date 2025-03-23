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

#' Analyze Cell Overlaps Between Nodes And Generate Heatmap
#'
#' Computes and visualizes the percentage overlap between nodes in a `FlowObject`, producing a heatmap of overlap percentages. Optionally, it can also display bar plots of specified parameters such as Enrichment, Entropy, or Average Proportion for the selected nodes.
#'
#' @param x A `FlowObject` after PruneGatingTree.
#' @param select_nodes Logical value indicating whether to select nodes interactively. Default is `FALSE`, which selects the top `n` nodes based on significance.
#' @param n Integer specifying the number of top nodes to select if `select_nodes` is `FALSE`. Default is `25`.
#' @param graphics Logical indicating whether to use graphical selection when `select_nodes` is `TRUE`. Default is `TRUE`.
#' @param margin Numeric value specifying the margins for the heatmap plot. Default is `8`.
#' @param parameter Character string specifying the parameter to plot alongside the heatmap. Can be `"Enrichment"`, `"Entropy"`, or `"Average Proportion"`. Default is `NULL`, which does not plot additional parameters.
#' @return A matrix containing the percentage overlap between the selected nodes.
#' @export
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' AnalyzeNodeOvealaps(x)
#'
#' # Select nodes interactively
#' AnalyzeNodeOvealaps(x, select_nodes = TRUE)
#' }
#' @importFrom gplots heatmap.2
#' @importFrom grDevices dev.new topo.colors
#' @importFrom graphics par barplot
#' @importFrom utils select.list
#' @family Unsupervised GatingTree Analysis

AnalyzeNodeOvealaps <- function(x, select_nodes = FALSE, n = NULL, graphics = TRUE, margin = 8, parameter = NULL){

    if(length(x@Gating$PrunedGatingTreeDF)==0){
        stop("Apply GatingTree analysis and Prune your tree before using this function. \n")
    }
    
    
    node <- x@Gating$GatingTreeObject
    sig_df <- x@Gating$PrunedGatingTreeDF
    sampledef_df <- x@sampledef$sampledef
    
    root_indices <- x@Gating$GatingTree$RootData$indices
    data <- x@Data[root_indices,]

    if(!is.null(n)){
        num <- min(n, nrow(sig_df))
    }else{
        num <- nrow(sig_df)
    }
    
    
    choices <- sig_df$markers_up_to_max
    entropy_sd <- scale(1/(10^sig_df$entropy))
    enrichment_sd <- scale(sig_df$max_enrichment)
    average_prop_sd <- scale(sig_df$average_proportion)
    sig_df$composite_score <- entropy_sd + enrichment_sd + average_prop_sd
    sig_df <- sig_df[order(sig_df$composite_score, decreasing = TRUE),]
    tmp <- gsub(".logdata", "", choices)
    tmp <- gsub(".neg", "-", tmp)
    tmp <- gsub(".pos", "+", tmp)
    choice_names <- gsub("_", "", tmp)
    
    tmp_char_df <- data.frame(Choices = choices, Choice_names = choice_names)


    if(select_nodes){
        selected <- select.list(choice_names, graphics = graphics, title = "Select Nodes to Plot", multiple =TRUE)
        logic <- tmp_char_df$Choice_names %in% selected
        selected_nodes <- tmp_char_df$Choices[logic]
        
        tmp <- gsub(".logdata", "", selected_nodes)
        tmp <- gsub(".neg", "-", tmp)
        tmp <- gsub(".pos", "+", tmp)
        node_names <- gsub("_", "", tmp)
        
        nodecells_list <- as.list(selected_nodes)
        names(nodecells_list) <- node_names

        nodecells_colour <- nodecells_list

        
        for(sel in selected_nodes){
            cellcounttable <- c()
            df <- c()
            df_summary <- c()
            path <- sel
            path <- unlist(strsplit(path, split='_'))
            tmp <- findNodeByPath(node, path)
            nodecells_list[[sel]]<- tmp$CurrentNodeIndices
            
            if(!is.null(parameter)){
                if(parameter == 'Enrichment'){
                    nodecells_colour[[sel]] <- tmp$CurrentEnrichment
                    
                }else{
                    if(parameter == 'Entropy'){
                        nodecells_colour[[sel]] <- tmp$CurrentEntropy
                        
                    } else {
                        nodecells_colour[[sel]] <- tmp$AverageProportion
                        
                    }
                    
                    
                }
            }
        }
        
        
    }else{
        selected_nodes <- tmp_char_df$Choices[1:num]
        
        tmp <- gsub(".logdata", "", selected_nodes)
        tmp <- gsub(".neg", "-", tmp)
        tmp <- gsub(".pos", "+", tmp)
        node_names <- gsub("_", "", tmp)
        
        nodecells_list <- as.list(1:num)
        names(nodecells_list) <- sig_df$markers_up_to_max[1:num]
        nodecells_colour <- nodecells_list
        for(i in 1:num){
            cellcounttable <- c()
            df <- c()
            df_summary <- c()
            path <- sig_df$markers_up_to_max[i]
            path <- unlist(strsplit(path, split='_'))
            tmp <- findNodeByPath(node, path)
            nodecells_list[[i]]<- tmp$CurrentNodeIndices
            if(!is.null(parameter)){
                if(parameter == 'Enrichment'){
                    nodecells_colour[[i]] <- tmp$CurrentEnrichment
                    
                }else{
                    if(parameter == 'Entropy'){
                        nodecells_colour[[i]] <- tmp$CurrentEntropy
                        
                    } else {
                        nodecells_colour[[i]] <- tmp$AverageProportion
                        
                    }
                    
                    
                }
            }
        }
    }
    
    calculate_overlap <- function(vector1, vector2) {
         intersection <- length(intersect(vector1, vector2))
        union <- length(union(vector1, vector2))
        
         if (union == 0) {
            return(0)
        } else {
            return((intersection / union) * 100)
        }
    }
    
    nodecells_names <- tmp_char_df$Choice_names[match(names(nodecells_list), tmp_char_df$Choices)]

    if(!is.null(parameter)){
        nodecells_colour <- unlist(nodecells_colour)
        names(nodecells_colour) <- nodecells_names
    }
    overlap_matrix <- matrix(nrow = length(nodecells_list), ncol = length(nodecells_list),
    dimnames = list(nodecells_names, nodecells_names))
    

    
    for (i in seq_along(nodecells_list)) {
        for (j in seq_along(nodecells_list)) {
            
            overlap_matrix[i, j] <- calculate_overlap(nodecells_list[[i]], nodecells_list[[j]])
        }
    }
    
    colours <- topo.colors(75)
    breaks <- seq(0, 100, length.out = length(colours) + 1)

    hm <- heatmap.2(overlap_matrix,
              col = colours,
              breaks = breaks,
              key = TRUE,
              symkey = FALSE,
              density.info = "none",
              trace = "none",
              cexRow = 0.8,
              cexCol = 0.8,
              margins = c(margin,margin))
              

              if(!is.null(parameter)){
                  dev.new()
                  par(mar = c(5, margin, 4, 3))
                  
                  nodecells_colour <- nodecells_colour[hm$rowInd]

                  if(parameter == 'Enrichment'){
                      title <-  'Enrichment Score'
                      barplot(nodecells_colour, las = 1, horiz = TRUE, main = title, xlab = title, cex.names = 0.8)
                  }else{
                      if(parameter == 'Entropy'){
                          title <-  'Gating Entropy'
                          barplot(nodecells_colour, las = 1, horiz = TRUE, main = title, xlab = title, cex.names = 0.8)
                      }else{
                          title <-  'Average Proportion'
                          nodecells_colour <- nodecells_colour*100
                          barplot(nodecells_colour, las = 1, horiz = TRUE, main = title, xlab = 'Average Percentage (%)', cex.names = 0.8)
                      }
                  }
        }
        
        out <- list(heatmap = hm, overlap_matrix = overlap_matrix, selected_nodes = choices)
        
        x@Gating$NodeOverlaps <- out
        return(invisible(x))
    
}
 

#' Extract Top Nodes from Node Overlaps in a FlowObject
#'
#' This function analyzes node overlaps stored within a `FlowObject` and extracts top nodes based on clustering and composite scores derived from entropy, enrichment, and average proportion.
#' It performs hierarchical clustering on the dendrogram of node overlaps, selects the best node from each cluster based on a composite score, and updates the gating tree with these nodes.
#'
#' @param x A `FlowObject` containing the gating tree and node overlaps data.
#' @param cluster_num An integer specifying the number of clusters to cut the dendrogram into. Defaults to 6.
#'
#' @return The modified `FlowObject` with updated `ExtractPrunedGatingTree` and `ExtractPrunedGatingTreePaths` slots indicating the pruned gating tree and the paths of selected top nodes.
#'
#' @examples
#' \dontrun{
#' x <- ExtractTopNodes(x, cluster_num = 6)
#' }
#'
#' @export
#'
#' @importFrom stats cutree as.hclust
#' @family Unsupervised GatingTree Analysis
ExtractTopNodes <- function(x, cluster_num = 6){
    output_node_overlaps<- x@Gating$NodeOverlaps
    row_clusters <- cutree(as.hclust(output_node_overlaps[[1]]$rowDendrogram), k = cluster_num)
    
    new_row_id <- row_clusters[output_node_overlaps[[1]]$rowInd]

    prn_df <- x@Gating$PrunedGatingTreeDF
    entropy_sd <- scale(1/(10^prn_df$entropy))
    enrichment_sd <- scale(prn_df$max_enrichment)
    average_prop_sd <- scale(prn_df$average_proportion)
    prn_df$composite_score <- entropy_sd + enrichment_sd + average_prop_sd
    
    tmp_best_path <- c()
    for(i in 1:cluster_num){
        
        tmp_paths <- names(new_row_id[new_row_id==i])
        tmp_paths <- gsub("\\+", ".logdata.pos_", tmp_paths)
        tmp_paths <- gsub("\\-", ".logdata.neg", tmp_paths)
        
        tmp_df <- prn_df[prn_df$markers_up_to_max %in% tmp_paths,]
        max_index <- which.max(tmp_df$composite_score)
        tmp_best_path <- c(tmp_best_path, tmp_df$markers_up_to_max[max_index])

    }
    
    node <- x@Gating$GatingTreeObject
    pruned_node <- update_nodes_by_paths(node = node, paths = tmp_best_path)
    pruned_node <- prune_tree(pruned_node)
    
    x@Gating$ExtractPrunedGatingTree <- pruned_node
    x@Gating$ExtractPrunedGatingTreePaths <- tmp_best_path
    
    return(x)
}

#' Convert Gating Tree to Graphical Representation
#'
#' This function converts a gating tree from a `FlowObject` into a graph object
#' which can then be visualized or exported. It supports different modes of operation,
#' allowing users to work with either all nodes, pruned nodes, or extracted top nodes
#' based on previous analyses.
#'
#' @param x A `FlowObject` containing the gating tree data.
#' @param mode A character string specifying which part of the gating tree to use:
#'   - 'All': Uses all nodes from the gating tree.
#'   - 'Pruned': Uses nodes pruned by `PruneGatingTree`.
#'   - 'Extracted': Uses top nodes extracted by `ExtractTopNodes`.
#'   The default is 'Pruned'.
#' @param size_factor A numeric value used to adjust the size of the nodes in the graph.
#'   Default is 1.
#' @param fontsize An integer specifying the font size used for node labels in the graph.
#'   Default is 20.
#' @param all_labels A logical indicating whether to label all nodes or not.
#'   Default is TRUE.
#'
#' @return A graph object that can be further processed with graph rendering or exporting functions.
#'
#' @examples
#' \dontrun{
#'   data(FlowObjectExample)
#'   graph <- convertGatingTreeToGraph(FlowObjectExample)
#'   render_graph(graph, width = 600, height = 600)
#'   export_graph(graph, file_name= "PrunedGatingTree.pdf")
#' }
#'
#' @export
#' @family GatingTree Visualization
convertGatingTreeToGraph <- function(x, mode = 'Pruned', size_factor=1, fontsize=20, all_labels = TRUE){
    
    if(mode=='Pruned'){
        
        node <- x@Gating$PrunedGatingTreeObject
    }
    if(mode=='All'){
        
        node <- x@Gating$GatingTreeObject
    }
    if(mode=='Extracted'){
        
        node <- x@Gating$ExtractPrunedGatingTree
    }

    datatree2 <- convertToDataTree(node)
    graph <- convert_to_diagrammer(datatree2, size_factor=1, fontsize=20, all_labels = all_labels)
    return(graph)
    
}

