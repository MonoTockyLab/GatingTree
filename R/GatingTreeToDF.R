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

#' Convert Gating Tree Object to Data Frames
#'
#' This function processes a Gating Tree object from a flow cytometry analysis package,
#' extracting and calculating various statistics such as enrichment, entropy, and delta
#' values for each gate. It returns the object with additional data frames attached that
#' contain detailed analysis results.
#'
#' @param x FlowObject.
#' @return FlowObject with a GatingTreeDF.
#' @details The function traverses the gating tree recursively to collect enrichment and
#'   entropy values from each node. It calculates deltas of enrichment values as differences
#'   between consecutive gates. Additionally, it computes a maximum enrichment statistic
#'   for each path through the tree, sorting these values to identify the most significant
#'   gates. Interaction effects (IE) and a detailed breakdown of changes in enrichment values
#'   are also calculated and stored in the returned object.
#' @examples
#' \dontrun{
#' x <- GatingTreeToDF(x)
#' }
#' @export
#' @importFrom stats setNames
#' @family GatingTree

GatingTreeToDF <- function(x){
    if(inherits(x,'FlowObject')){
        node <- x@Gating$GatingTreeObject
    }else{
        node <- x
    }
    enrichments <- collect_all_enrichment(node)
    entropies <- collect_all_entropy(node)
    d_values <- sapply(enrichments, function(x) x[1], simplify = FALSE)
    d_vector <- unlist(d_values)
    d_vector <- d_vector[!duplicated(d_vector)]
    
    delta_E <- list()
    for (i in 1:length(enrichments)) {
        current_enrichments <- unlist(enrichments[i])
        l = length(current_enrichments)
        tmp = c()
        for(j in 2:l){
            tmp <- c(tmp, current_enrichments[j] - current_enrichments[j-1])
            names(tmp)[j-1] <- names(current_enrichments[j])
        }
        
        delta_E[[i]] <- tmp
        names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
    }
    
    calculated_deltas <- list()
    for (i in 1:length(enrichments)) {
        current_enrichments <- enrichments[[i]]
        l <- length(current_enrichments)
        if(l <=1){
            if(!is.na(current_enrichments)){
                delta_E[[i]] <- current_enrichments
            }else{
                delta_E[[i]] <- 0
            }
            
            names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
            next
        }
        tmp <- vector("numeric", l - 1)
        names_tmp <- vector("character", l - 1)
        
        k <- 1
        for (j in 2:l) {
            marker_pair <- paste(names(current_enrichments)[1:j], collapse = "_")
            
            if (!marker_pair %in% names(calculated_deltas)) {
                tmp[k] <- current_enrichments[j] - current_enrichments[j - 1]
                names_tmp[k] <- names(current_enrichments)[j]
                calculated_deltas[[marker_pair]] <- TRUE  # Mark this pair as calculated
                k <- k + 1
            }
        }
        
        if (k > 1) {
            tmp <- tmp[1:(k - 1)]  # trim the unused portion
            names_tmp <- names_tmp[1:(k - 1)]
            names(tmp) <- names_tmp
            delta_E[[i]] <- tmp
            names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
        }
    }

    
    marker_deltas <- list()
    for (i in 1:length(delta_E)) {
        for (marker in names(delta_E[[i]])) {
            if (!is.null(marker_deltas[[marker]])) {
                marker_deltas[[marker]] <- c(marker_deltas[[marker]], delta_E[[i]][[marker]])
            } else {
                marker_deltas[[marker]] <- delta_E[[i]][[marker]]
            }
        }
    }
    
    delta_df <- do.call(rbind, lapply(names(marker_deltas), function(marker) {
        data.frame(marker = marker, delta = marker_deltas[[marker]])
    }))
    
   # rownames(delta_df) <- NULL

    
    process_enrichment <- function(enrichment, index, entropies) {
        max_value <- max(enrichment)
        max_position <- which.max(enrichment)
        markers_up_to_max <- names(enrichment)[1:max_position]
        markers_sorted <- sort(markers_up_to_max)
        combination_name <- paste(markers_sorted, collapse = "_")
        markers_up_to_max = paste(markers_up_to_max, collapse= '_')
        
        corresponding_entropy <- entropies[[index]][max_position]
        
        list(
        max_value = max_value,
        combination = combination_name,
        markers_up_to_max = markers_up_to_max,
        entropy = corresponding_entropy
        )
    }
    
    maxima_results <- mapply(process_enrichment, enrichments, seq_along(enrichments), MoreArgs = list(entropies = entropies), SIMPLIFY = FALSE)
    maxima_df <- do.call(rbind, lapply(maxima_results, function(x) {
        data.frame(max_enrichment = x$max_value,entropy = x$entropy,  combination = x$combination, markers_up_to_max = x$markers_up_to_max, stringsAsFactors = FALSE)
    }))
    
    rownames(maxima_df) <- NULL
    maxima_df <- maxima_df[!duplicated(maxima_df$combination),]
    sl <- order(maxima_df$max_enrichment, decreasing = TRUE)
    maxima_df <- maxima_df[sl,]

    marker_deltas_df <- do.call(rbind, lapply(names(marker_deltas), function(marker) {
      data.frame(Marker = marker, DeltaE_value = marker_deltas[[marker]])
    }))
    
    if(inherits(x,'FlowObject')){
        x@Gating$GatingTree$delta_df <- marker_deltas_df
        x@Gating$GatingTreeDF <- maxima_df
        return(x)

    }else{
        return(maxima_df)
    }

}

#' Plot Delta Enrichment or Interaction Effects for Each Marker
#'
#' This function takes a complex object containing gating tree data and performs
#' statistical tests to evaluate the significance of delta enrichment or interaction
#' effects across markers. It plots the results as boxplots and returns the modified
#' object with additional annotations based on the analysis.
#'
#' @param x FlowObject.
#' @return Returns the object `x` with statistical values
#' @details The function conducts Kruskal-Wallis tests to determine the significance of
#'   differences across markers, followed by Dunn's test for post-hoc analysis with Bonferroni
#'   correction if significant.
#' @examples
#' \dontrun{
#'   x <- PlotDeltaEnrichment(x)
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme
#' @importFrom dunn.test dunn.test
#' @importFrom rlang sym
#' @family GatingTree

PlotDeltaEnrichment <- function(x){
    marker_deltas_df <- x@Gating$GatingTree$delta_df
    kruskal_test_delta <- kruskal.test(marker_deltas_df[['DeltaE_value']] ~ marker_deltas_df[['Marker']])
    
    if (kruskal_test_delta$p.value < 0.05) {
        posthoc_results_delta <- dunn.test(marker_deltas_df$DeltaE_value, marker_deltas_df$Marker, method = "bonferroni", table = FALSE)
    }
    
    mean_delta_values <- aggregate(DeltaE_value ~ Marker, data = marker_deltas_df, FUN = mean)
    colnames(mean_delta_values) <- c("Marker", "MeanDeltaE_value")

    p_vec <- posthoc_results_delta$P.adjusted
    names(p_vec) <- posthoc_results_delta$comparisons
    markers <- names(p_vec)
    markers <- strsplit(markers, " - ")
    marker_counts <- table(unlist(markers))
    summary_table <- data.frame(
    Marker = names(marker_counts),
    Count_of_Significant_Differences = as.integer(marker_counts)
    )
    summary_table <- summary_table[order(-summary_table$Count_of_Significant_Differences), ]
    summary_table <- merge(summary_table, mean_delta_values, by = "Marker")
    summary_table <- summary_table[order(-summary_table$MeanDeltaE_value), ]
    
    markers <- summary_table$Marker
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](pos)", "+", markers))
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](neg)", "-", markers))
    summary_table$Marker <- markers
    
    markers <- marker_deltas_df[,1]
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](pos)", "+", markers))
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](neg)", "-", markers))
    marker_deltas_df$Marker <- factor(markers, levels = summary_table$Marker[order(-summary_table$MeanDeltaE_value)])

    p <- ggplot(marker_deltas_df, aes(x = !!sym("Marker"), y = !!sym("DeltaE_value"), fill = !!sym("Marker"))) +
    geom_boxplot()+
    labs(title = "Delta(Enrichment) Values for Each Marker",
    x = "Marker",
    y = "Delta Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plot(p)
    x@Gating$GatingTree$DeltaPlot <- p
    x@Gating$GatingTree$DeltaStats <- kruskal_test_delta
    return(x)
}

#' Plot Delta Enrichment or Interaction Effects for Each Marker Using Pruned GatingTree
#'
#' @param x FlowObject after applying PruneGatingTree.
#' @param q A numeric value between 0 and 1 as a quantile threshold for extracting top markers.
#' @return Returns the object `x` with statistical values
#' @details The function conducts Kruskal-Wallis tests to determine the significance of
#'   differences across markers, followed by Dunn's test for post-hoc analysis with Bonferroni
#'   correction if significant.
#' @examples
#' \dontrun{
#'   x <- PlotDeltaEnrichmentPrunedTree(x)
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme
#' @importFrom dunn.test dunn.test
#' @importFrom rlang sym
#' @family GatingTree

PlotDeltaEnrichmentPrunedTree <- function(x, q = 0){
    if(length(x@Gating$PrunedGatingTreeObject)==0){
        stop('Apply PruneGatingTree before using this function. \n')
        
    }
    
    node <- x@Gating$PrunedGatingTreeObject
  
    enrichments <- collect_all_enrichment(node)
    entropies <- collect_all_entropy(node)
    d_values <- sapply(enrichments, function(x) x[1], simplify = FALSE)
    d_vector <- unlist(d_values)
    d_vector <- d_vector[!duplicated(d_vector)]
    
    delta_E <- list()
    for (i in 1:length(enrichments)) {
        current_enrichments <- unlist(enrichments[i])
        l = length(current_enrichments)
        tmp = c()
        for(j in 2:l){
            tmp <- c(tmp, current_enrichments[j] - current_enrichments[j-1])
            names(tmp)[j-1] <- names(current_enrichments[j])
        }
        
        delta_E[[i]] <- tmp
        names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
    }
    
    calculated_deltas <- list()
    for (i in 1:length(enrichments)) {
        current_enrichments <- enrichments[[i]]
        l <- length(current_enrichments)
        if(l <=1){
            if(!is.na(current_enrichments)){
                delta_E[[i]] <- current_enrichments
            }else{
                delta_E[[i]] <- 0
            }
            
            names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
            next
        }
        tmp <- vector("numeric", l - 1)
        names_tmp <- vector("character", l - 1)
        
        k <- 1
        for (j in 2:l) {
            marker_pair <- paste(names(current_enrichments)[1:j], collapse = "_")
            
            if (!marker_pair %in% names(calculated_deltas)) {
                tmp[k] <- current_enrichments[j] - current_enrichments[j - 1]
                names_tmp[k] <- names(current_enrichments)[j]
                calculated_deltas[[marker_pair]] <- TRUE  # Mark this pair as calculated
                k <- k + 1
            }
        }
        
        if (k > 1) {
            tmp <- tmp[1:(k - 1)]  # trim the unused portion
            names_tmp <- names_tmp[1:(k - 1)]
            names(tmp) <- names_tmp
            delta_E[[i]] <- tmp
            names(delta_E)[i] <- paste(names(current_enrichments), collapse = '_')
        }
    }
    
    marker_deltas <- list()
    for (i in 1:length(delta_E)) {
        for (marker in names(delta_E[[i]])) {
            if (!is.null(marker_deltas[[marker]])) {
                marker_deltas[[marker]] <- c(marker_deltas[[marker]], delta_E[[i]][[marker]])
            } else {
                marker_deltas[[marker]] <- delta_E[[i]][[marker]]
            }
        }
    }
    
    delta_df <- do.call(rbind, lapply(names(marker_deltas), function(marker) {
        data.frame(marker = marker, delta = marker_deltas[[marker]])
    }))
    
    process_enrichment <- function(enrichment, index, entropies) {
        max_value <- max(enrichment)
        max_position <- which.max(enrichment)
        markers_up_to_max <- names(enrichment)[1:max_position]
        markers_sorted <- sort(markers_up_to_max)
        combination_name <- paste(markers_sorted, collapse = "_")
        markers_up_to_max = paste(markers_up_to_max, collapse= '_')
        
        corresponding_entropy <- entropies[[index]][max_position]
        
        list(
        max_value = max_value,
        combination = combination_name,
        markers_up_to_max = markers_up_to_max,
        entropy = corresponding_entropy
        )
    }
    
    maxima_results <- mapply(process_enrichment, enrichments, seq_along(enrichments), MoreArgs = list(entropies = entropies), SIMPLIFY = FALSE)
    maxima_df <- do.call(rbind, lapply(maxima_results, function(x) {
        data.frame(max_enrichment = x$max_value,entropy = x$entropy,  combination = x$combination, markers_up_to_max = x$markers_up_to_max, stringsAsFactors = FALSE)
    }))
    
    rownames(maxima_df) <- NULL
    maxima_df <- maxima_df[!duplicated(maxima_df$combination),]
    sl <- order(maxima_df$max_enrichment, decreasing = TRUE)
    maxima_df <- maxima_df[sl,]
    
    marker_deltas_df <- do.call(rbind, lapply(names(marker_deltas), function(marker) {
        data.frame(Marker = marker, DeltaE_value = marker_deltas[[marker]])
    }))
        
    kruskal_test_delta <- kruskal.test(marker_deltas_df[['DeltaE_value']] ~ marker_deltas_df[['Marker']])
    
    
    posthoc_results_delta <- dunn.test(marker_deltas_df$DeltaE_value, marker_deltas_df$Marker, method = "bonferroni", table = FALSE)

    mean_delta_values <- aggregate(DeltaE_value ~ Marker, data = marker_deltas_df, FUN = mean)
    colnames(mean_delta_values) <- c("Marker", "MeanDeltaE_value")
    p_vec <- posthoc_results_delta$P.adjusted
    names(p_vec) <- posthoc_results_delta$comparisons

    markers <- names(p_vec)
   
    
    markers <- strsplit(markers, " - ")
    marker_counts <- table(unlist(markers))
    summary_table <- data.frame(
    Marker = names(marker_counts),
    Count_of_Significant_Differences = as.integer(marker_counts)
    )
    summary_table <- summary_table[order(-summary_table$Count_of_Significant_Differences), ]
    summary_table <- merge(summary_table, mean_delta_values, by = "Marker")
    summary_table <- summary_table[order(-summary_table$MeanDeltaE_value), ]
    
    markers <- summary_table$Marker
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](pos)", "+", markers))
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](neg)", "-", markers))
    summary_table$Marker <- markers
    
    markers <- marker_deltas_df[,1]
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](pos)", "+", markers))
    markers <- gsub("(.logdata)|[,]", "", gsub("[.](neg)", "-", markers))
    marker_deltas_df$Marker <- factor(markers, levels = summary_table$Marker[order(-summary_table$MeanDeltaE_value)])
    
    
    if(!is.null(q) && is.numeric(q) && length(q)==1 && q >=0 & q <= 1 ){
        quantile_val <- tapply(marker_deltas_df$DeltaE_value, marker_deltas_df$Marker, mean)
        threshold_q <- quantile(quantile_val, q)
        sel_markers <- names(quantile_val[quantile_val > threshold_q])
        logic <- marker_deltas_df$Marker %in% sel_markers
        marker_deltas_df <- marker_deltas_df[logic, ]
        
        
    }
    
    
    p <- ggplot(marker_deltas_df, aes(x = !!sym("Marker"), y = !!sym("DeltaE_value"), fill = !!sym("Marker"))) +
    geom_boxplot()+
    labs(title = "Delta(Enrichment) Values for Each Marker",
    x = "Marker",
    y = "Delta Value") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for clarity
    plot(p)
    
    
    if(length(x@Gating$PrunedGatingTree)==0){
        x@Gating$PrunedGatingTree <- list()
    }
    x@Gating$PrunedGatingTree$delta_df <- marker_deltas_df
    x@Gating$PrunedGatingTree$DeltaPlot <- p
    x@Gating$PrunedGatingTree$DeltaStats <- kruskal_test_delta
    
    return(x)
}

#' Add Prune Markers to Nodes
#'
#' This function recursively adds prune markers to nodes of a tree based on specific
#' conditions related to enrichment, entropy, and average proportion.
#'
#' @param node The node to process.
#' @param min_average_proportion The minimum average proportion required to not prune a node.
#' @param theta The threshold added to the current enrichment for comparison.
#' @return The node with prune markers added.
#' @examples
#' \dontrun{
#'   node <- add_prune(node, 0.001, 0.02)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd

add_prune <- function(node,  min_average_proportion = 0.001, theta = 0) {
    if(!is.null(node$CurrentNodeName) && node$CurrentNodeName=='rootNode'){
        node[['prune']] <- FALSE
    }
    node$prune <- FALSE
    
    
    if (!is.null(node$Children) && length(node$Children) > 0) {
        for (child_name in names(node$Children)) {
            
            node$Children[[child_name]] <- add_prune(node$Children[[child_name]],  min_average_proportion, theta)
            
            logic <- (node$Children[[child_name]]$CurrentEnrichment >=  node$CurrentEnrichment + theta) &&
            (node$Children[[child_name]]$CurrentEntropy <= node$CurrentEntropy) && (node$Children[[child_name]]$AverageProportion > min_average_proportion)
            
            node$Children[[child_name]][['prune']] <- !logic
             

            
        }
    }
    
    return(node)
}

#' Prune Tree
#'
#' Recursively prunes a tree based on prune markers set by `add_prune`.
#'
#' @param node The root node of the tree.
#' @return The pruned tree or NULL if the root node is pruned.
#' @examples
#' \dontrun{
#'   pruned_tree <- prune_tree(root_node)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd

prune_tree <- function(node) {
    # Remove node if it's marked for pruning
    if (is.null(node$prune)) {
        return(NULL)
    }else{if(node$prune){return(NULL)}}
    
    # Process children if not pruning this node
    if (!is.null(node$Children)) {
        temp_children <- list()
        for (child_name in names(node$Children)) {
            pruned_child <- prune_tree(node$Children[[child_name]])
            if (!is.null(pruned_child)) {
                temp_children[[child_name]] <- pruned_child
            }
        }
        node$Children <- temp_children
    }
    
    return(node)
}

#' Update Nodes Based on External Data
#'
#' This function updates nodes to set prune flags based on a comparison with an external
#' data frame that specifies which node combinations should not be pruned.
#'
#' @param node The node to update.
#' @param path The current path of nodes, initially NULL.
#' @param res_df The data frame containing combinations that should not be pruned.
#' @return The node with updated prune flags.
#' @examples
#' \dontrun{
#'   updated_node <- find_and_update_nodes(node, NULL, res_df)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd

find_and_update_nodes <- function(node, path = NULL, res_df) {
    if (is.null(path)) path <- c()

    if (!is.null(node$CurrentPath)) {
        path <- node$CurrentPath
    }
    current_combination <- paste(path[-1], collapse="_")
    
    if (current_combination %in% res_df$combination) {
        node$prune <- FALSE
    }else{
        node$prune <- TRUE
    }
    if (!is.null(node$Children)) {
        for (child_name in names(node$Children)) {
            node$Children[[child_name]] <- find_and_update_nodes(node$Children[[child_name]], path, res_df)
            if (!node$Children[[child_name]]$prune) {
                node$prune <- FALSE
            }
        }
    }

    return(node)
}

#' Prune a Gating Tree Based on Statistical Criteria
#'
#' This function prunes a gating tree by applying various statistical thresholds to
#' the nodes based on entropy, enrichment, and average proportion metrics. Nodes that
#' do not meet the specified criteria are pruned from the tree. Additionally, p-values
#' are adjusted for multiple comparisons.
#' In the current implementation the output PrunedGatingTreeNodePercentages will have
#' nodes filtered only through a threshold for average proportion.
#'
#' @param x An object, expected to be of class 'FlowObject', containing
#'   gating tree data and metadata.
#' @param max_entropy Maximum allowable entropy for a node to remain in the gating tree.
#' @param min_enrichment Minimum enrichment required for a node to remain in the gating tree.
#' @param min_average_proportion Minimum average proportion of cells required for a node to
#'   remain in the gating tree.
#' @param p_adjust_method A character string indicating the method to be used for adjusting
#'   p-values for multiple comparisons. Defaults to 'BY' (Benjamini-Yekutieli).
#' @param theta A numeric threshold added to the enrichment values for each node.
#'   This threshold is used to enforce a criterion where only nodes showing a steady
#'   increase in enrichment, greater than this threshold, can be considered for retention
#'   in the pruned gating tree. The default is zero.
#' @return Returns the modified object 'x' with the gating tree pruned according to the
#'   specified parameters. The function also attaches the pruned gating tree and a data
#'   frame containing node statistics to the object.
#' @details The function first identifies nodes that meet specified entropy and enrichment
#'   criteria, computes statistical metrics for these nodes, and then prunes the gating
#'   tree based on these metrics and average proportion criteria. Adjusted p-values are
#'   calculated to account for multiple testing.
#' @examples
#' \dontrun{
#'   updated_object <- PruneGatingTree(x, min_enrichment = 0.5,max_entropy =0.5)
#' }
#' @importFrom dplyr count full_join mutate coalesce
#' @importFrom rlang sym
#' @export
#' @family GatingTree


PruneGatingTree <- function(x, max_entropy =0.9, min_enrichment = 0.1, min_average_proportion = 0.001, p_adjust_method = 'BY', theta = NULL){
    
    if(length(x@Gating$GatingTreeDF)==0){
        stop('Apply GatingTreeToDF before using this function. \n')
    }
    
    unique_maxima_df <- x@Gating$GatingTreeDF
    node <- x@Gating$GatingTreeObject
    sampledef <- x@sampledef$sampledef
    
    logic <-(unique_maxima_df$entropy <= max_entropy)&( unique_maxima_df$max_enrichment >= min_enrichment)
    unique_maxima_df <- unique_maxima_df[logic,]
    
    rootdata <-  node$RootData
    total_cell_per_file <- table(rootdata$file)
    total_cell_per_file <- data.frame(file = names(total_cell_per_file), total = as.vector(total_cell_per_file))
    average_proportion <- p_value_vec <- rep(0, nrow(unique_maxima_df))
    percent_list <- as.list(p_value_vec)
    names(percent_list) <- unique_maxima_df$markers_up_to_max#combination
    
    for(k in 1:nrow(unique_maxima_df)){
        path <- unique_maxima_df$markers_up_to_max[k]
        path <- unlist(strsplit(path, split='_'))
        
        tmp <- findNodeByPath(node, path)
        node_data <- rootdata[rootdata$indices %in% tmp$CurrentNodeIndices,]
        total_cell_num <- node_data %>%
          count(file, name = "node_total") %>%
          full_join(data.frame(total_cell_per_file), by = "file") %>%
          full_join(sampledef, by = "file") %>%
          mutate(node_total = coalesce(!!sym("node_total"), 0),  # Replace NA values with 0
                 node_percentage = !!sym("node_total") / !!sym("total"))
                 
        percent_list[[k]] <- total_cell_num
        p_value_vec[k] = p_value = wilcox.test(data = total_cell_num, node_percentage~group)$p.value
        total_cell_num$node_percentage[is.na(total_cell_num$node_percentage)] <- 0
        average_proportion[k] <- mean(total_cell_num$node_percentage)
        
    }

    res_df <- cbind(unique_maxima_df, data.frame(p_value = p_value_vec, average_proportion = average_proportion))
   # res_df <- res_df[!is.na(res_df$p_value),]
    res_df[['p_adjust']] <- p.adjust(res_df$p_value, method = p_adjust_method)
    
    if(!is.null(min_average_proportion) && !is.null(theta)){
        pruned_node <- add_prune(node, min_average_proportion = min_average_proportion, theta = theta)
        
    }else{
        if(!is.null(min_average_proportion) && is.null(theta)){
            pruned_node <- add_prune(node, min_average_proportion = min_average_proportion, theta = 0)
        }else{
            pruned_node <- add_prune(node, min_average_proportion = 0, theta = 0)
        }
    }
    
    
    tree_to_df <- function(node){
        enrichments <- collect_all_enrichment(node)
        entropies <- collect_all_entropy(node)
        
        process_enrichment <- function(enrichment, index, entropies) {
            max_value <- max(enrichment)
            max_position <- which.max(enrichment)
            markers_up_to_max <- names(enrichment)[1:max_position]
            markers_sorted <- sort(markers_up_to_max)
            combination_name <- paste(markers_sorted, collapse = "_")
            markers_up_to_max <- paste(markers_up_to_max, collapse= '_')
            
            corresponding_entropy <- entropies[[index]][max_position]
            
            list(
            max_value = max_value,
            combination = combination_name,
            markers_up_to_max = markers_up_to_max,
            entropy = corresponding_entropy
            )
        }
        
        maxima_results <- mapply(process_enrichment, enrichments, seq_along(enrichments), MoreArgs = list(entropies = entropies), SIMPLIFY = FALSE)
        maxima_df <- do.call(rbind, lapply(maxima_results, function(x) {
            data.frame(max_enrichment = x$max_value,entropy = x$entropy,  combination = x$combination, markers_up_to_max = x$markers_up_to_max, stringsAsFactors = FALSE)
        }))
        
        return(maxima_df)
    }
    
    df <- tree_to_df(pruned_node)
    res_df <- res_df[res_df$markers_up_to_max %in% df$markers_up_to_max,]#[res_df$combination %in% df$combination,]
    
    
    
    logic <- res_df$max_enrichment > min_enrichment & res_df$entropy < max_entropy & res_df$average_proportion >= min_average_proportion

    
    pruned_node <- find_and_update_nodes(pruned_node, res_df = res_df)
    pruned_node <- prune_tree(pruned_node)

    x@Gating$PrunedGatingTreeObject <- pruned_node
    x@Gating$PrunedGatingTreeDF <- res_df
    
    sampledef$file <- as.character(sampledef$file)
    extract_and_merge <- function(df, idx) {
        df$file <- as.character(df$file)
        merged_df <- merge(sampledef, df[, c("file", "node_percentage")], by = "file", all.x = TRUE)
        names(merged_df)[names(merged_df) == "node_percentage"] <-idx
        return(merged_df[, c("file", idx)])
    }
    
    merged_df <- sampledef
    for(i in 1:length(percent_list)){
        tdf <- percent_list[[i]]
        colnames(tdf)[colnames(tdf)=='node_percentage'] <- names(percent_list)[i]
        merged_df <- merge(merged_df, tdf[, c("file", names(percent_list)[i])], by = "file", all.x = TRUE)
    }
    
    x@Gating$PrunedGatingTreeNodePercentages <- merged_df
    
    return(x)
    
}

#' Extract All Nodes from GatingTree
#'
#' This function recursively extracts all nodes from the specified GatingTree object
#' within a FlowObject, either the full GatingTree or the pruned version.
#'
#' @param x A FlowObject containing a GatingTree.
#' @param pruned Logical indicating whether to extract nodes from the pruned GatingTree.
#'   Defaults to FALSE.
#' @return A data frame containing details of all nodes in the specified GatingTree.
#' @export
#' @examples
#' \dontrun{
#'   full_tree_nodes <- extractNodes(x)
#'   pruned_tree_nodes <- extractNodes(x, pruned = TRUE)
#' }
extractNodes <- function(x, pruned = FALSE) {
    if (!pruned) {
      if (is.null(x@Gating$GatingTreeObject) || length(x@Gating$GatingTreeObject) == 0) {
        stop("GatingTreeObject is empty.")
      }
      node <- x@Gating$GatingTreeObject
    } else {
      if (is.null(x@Gating$PrunedGatingTreeObject) || length(x@Gating$PrunedGatingTreeObject) == 0) {
        stop("PrunedGatingTreeObject is empty.")
      }
      node <- x@Gating$PrunedGatingTreeObject
    }
    
  return(node)
}


#' Convert Gating Tree to DiagrammeR Graph
#'
#' This function takes a gating tree object and constructs a graphical representation
#' using the DiagrammeR package. It visually represents the enrichment, entropy,
#' and optionally the average proportion of nodes in the gating tree.
#'
#' @param root The root node of the gating tree object, which must have properties like
#'   `isRoot`, `name`, `CurrentEnrichment`, `CurrentEntropy`, and optionally `AverageProportion`.
#' @param size_factor A multiplier for node sizes in the graph, allowing customization based on
#'   enrichment or average proportion values.
#' @param fontsize The size of font to be used in the graph.
#' @param average_proportion A logical flag indicating whether to use the average proportion
#'   of nodes to adjust node sizes and color gradient based on enrichment values.
#' @param all_labels Logical. If TRUE, all parameters are included in each node within the output graph.
#' @return Returns a DiagrammeR graph object representing the gating tree with nodes colored
#'   and sized according to specified metrics.
#' @details The function recursively traverses the gating tree, starting from the root,
#'   adding nodes to the graph with labels that include relevant metrics. Node size and color
#'   are determined by either the enrichment or entropy values, depending on the `average_proportion` flag.
#' @examples
#' \dontrun{
#'   # Assuming 'root' is your gating tree root node
#'   graph <- convert_to_diagrammer(root, size_factor = 1, average_proportion = TRUE)
#'   # To view the graph
#'   DiagrammeR::render_graph(graph)
#' }
#' @export
#' @importFrom DiagrammeR create_graph add_node add_edge render_graph
#' @family GatingTree Visualization

convert_to_diagrammer <- function(root, size_factor = 1, fontsize = 12, average_proportion = FALSE, all_labels = TRUE) {
    graph <- create_graph(directed = TRUE)
    node_index <- 1 # To keep track of nodes
    add_nodes_recursively <- function(node, parent_index = NULL) {
        if(node$isRoot) {
            label_text <- 'Root'
            node_size <- size_factor * 0.1
            tmp <- 0
            color_value <- rgb(tmp, 0, 1-tmp, alpha = tmp)  # Red gradient

            graph <<- add_node(graph, label = label_text,
                               node_aes = list(title = node$name, fillcolor=color_value, width = node_size, fontcolor = 'black', fontsize = fontsize))
        } else {
            formatted_entropy <- round(node$CurrentEntropy, 2)
            if (abs(formatted_entropy) < 0.001) formatted_entropy <- 0
            simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](pos)$", "+", node$name))
            simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](neg)$", "-", simplified_name))

            label_text <- if (all_labels) {
                sprintf("%s\nEnrichment: %.2f\nEntropy: %.2f\nAvg. %%: %.1f%%", simplified_name, node$CurrentEnrichment, formatted_entropy, node$AverageProportion * 100)
            } else {
                simplified_name
            }

            if(average_proportion) {
                node_size <- size_factor * node$AverageProportion
                tmp <- node$CurrentEnrichment
                tmp <- max(tmp, 0)
                tmp <- pmin(tmp, 1)
                color_value <- rgb(tmp, 0, 1-tmp, alpha = tmp)
            } else {
                node_size <- size_factor * node$CurrentEnrichment
                tmp <- node$CurrentEntropy
                tmp <- max(tmp, 0)
                tmp <- pmin(tmp, 1)
                color_value <- rgb(1-tmp, 0, tmp, alpha = (1-tmp) / 1.5)
            }

            graph <<- add_node(graph, label = label_text,
                               node_aes = list(title = node$name, fillcolor=color_value, width = node_size, fontcolor = 'black', fontsize = fontsize))

            if (!is.null(parent_index)) {
                graph <<- add_edge(graph, from = parent_index, to = node_index)
            }
        }
        current_index <- node_index
        node_index <<- node_index + 1
        if (!node$isLeaf) {
            for (child in node$children) {
                parent_index <<- current_index
                add_nodes_recursively(child, current_index)
            }
        }
    }
    add_nodes_recursively(root)
    return(graph)
}


#' Collect All Enrichment Values Recursively
#'
#' This function recursively collects enrichment values from a gating tree structure.
#' Starting from the provided node, it traverses all child nodes to aggregate their
#' enrichment values into a list.
#'
#' @param tree A tree node from which to start collecting enrichment values. The
#'   node structure is expected to have a `History` component containing `enrichment`
#'   values and possibly `Children` with further sub-nodes.
#' @return A list of enrichment values collected recursively from the node and its children.
#' @details Enrichment values are extracted from the `History` component of each node,
#'   if available. The function recursively explores all child nodes to collect their
#'   enrichment values.
#' @examples
#' \dontrun{
#'   enrichment_values <- collect_all_enrichment(root_node)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd

collect_all_enrichment <- function(tree) {
  enrichment_list <- list()

  if (!is.null(tree$History) && !is.null(tree$History$enrichment)) {
    enrichment_list <- c(enrichment_list, list(tree$History$enrichment))
  }

  if (length(tree$Children) > 0) {
    for (child_name in names(tree$Children)) {
      child_enrichments <- collect_all_enrichment(tree$Children[[child_name]])
      enrichment_list <- c(enrichment_list, child_enrichments)
    }
  }

  return(enrichment_list)
}



collect_peak_nodes <- function(tree) {
  node_names_with_peak <- list()

  if (!is.null(tree$Peak) && !is.null(tree$CurrentPath)) {
      tpath = paste(tree$CurrentPath, collapse = '_')
    node_names_with_peak <- c(node_names_with_peak, tpath)
  }

  if (length(tree$Children) > 0) {
    for (child_name in names(tree$Children)) {
      child_node_names <- collect_peak_nodes(tree$Children[[child_name]])
      node_names_with_peak <- c(node_names_with_peak, child_node_names)
    }
  }

  return(node_names_with_peak)
}


#' Collect All Entropy Values Recursively
#'
#' Similar to `collect_all_enrichment`, this function recursively collects entropy
#' values from a specified node and all its child nodes within a gating tree.
#'
#' @param tree A tree node from which to start collecting entropy values.
#'   Nodes are expected to have a `History` with `entropy` values and possibly
#'   `Children` with further sub-nodes.
#' @return A list of entropy values collected recursively from the node and its children.
#' @details Entropy values are gathered from the `History` component of each node,
#'   if present. The function recursively navigates through all child nodes to collect
#'   their entropy values.
#' @examples
#' \dontrun{
#'   entropy_values <- collect_all_entropy(root_node)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd

collect_all_entropy <- function(tree) {
  entropy_list <- list()

  if (!is.null(tree$History) && !is.null(tree$History$entropy)) {
    entropy_list <- c(entropy_list, list(tree$History$entropy))
  }

  if (length(tree$Children) > 0) {
    for (child_name in names(tree$Children)) {
      child_entropys <- collect_all_entropy(tree$Children[[child_name]])
      entropy_list <- c(entropy_list, child_entropys)
    }
  }

  return(entropy_list)
}




#' Update Gating Tree Nodes Based on Valid Paths
#'
#' This function traverses a gating tree and flags nodes for pruning if their
#' combination (derived from the node's current path) is not present in a given set
#' of valid paths. The function recursively updates each child node, and if any child node
#' is retained (i.e., not pruned), the parent node is also retained.
#'
#' @param node A node object from a gating tree. The node should contain a \code{CurrentPath}
#'   field (a character vector) and a \code{Children} list.
#' @param paths A character vector of valid node combination strings (paths) that should
#'   be retained. Nodes whose combination is not in this vector will be flagged with \code{prune = TRUE}.
#'
#' @return The updated node object with its \code{prune} flag set to \code{TRUE} for nodes
#'   not in the valid paths and \code{FALSE} for nodes that are valid or have valid descendants.
#'
#' @examples
#' \dontrun{
#' # Define valid paths (e.g., combinations of markers) to keep:
#' valid_paths <- c("CD4.logdata.pos_CD8.logdata.neg", "CD4.logdata.pos_CD19.logdata.pos")
#' # Update the gating tree starting from the root node:
#' updated_tree <- update_nodes_by_paths(rootNode, valid_paths)
#' }
#'
#' @keywords internal
#' @family GatingTree
#' @noRd
update_nodes_by_paths <- function(node, paths) {
    if (!is.null(node$CurrentPath)) {
        current_combination <- paste(node$CurrentPath[-1], collapse = "_")
    } else {
        current_combination <- ""
    }
    if (current_combination %in% paths) {
        node$prune <- FALSE
    } else {
        node$prune <- TRUE
    }
    if (!is.null(node$Children) && length(node$Children) > 0) {
        for (child_name in names(node$Children)) {
            node$Children[[child_name]] <- update_nodes_by_paths(node$Children[[child_name]], paths)
            if (!node$Children[[child_name]]$prune) {
                node$prune <- FALSE
            }
        }
    }
    
    return(node)
}



#' Extract and Finalize the Gating Tree Object
#'
#' This function performs a tree extraction by specifying node paths. The resulting tree will be stored in the \code{ExtractedGatingTreeObject} slot
#' of the FlowObject.
#'
#' @param x A FlowObject that has been previously processed with \code{GatingTreeToDF}
#'   and contains gating information in the \code{Gating} slot, including
#'   \code{GatingTreeDF} and \code{PrunedGatingTreeObject}.
#' @param node_paths A character vector of node paths to be used in the extraction
#'   process.
#'
#' @return The updated FlowObject with a new slot \code{Gating$ExtractedGatingTreeObject}
#'   containing the final, pruned gating tree.
#'
#' @details
#' Node paths should use the format in GatingTreeDF's \code{markers_up_to_max} slot.
#'
#' @examples
#' \dontrun{
#' # Assume 'x' is a FlowObject that has been processed with GatingTreeToDF and PruneGatingTree:
#' x <- ExtractGatingTree(x, node_paths = c("CD4.logdata.pos_CD8.logdata.neg",
#' "CD4.logdata.pos_CD19.logdata.pos"))
#' }
#'
#' @export
#' @family GatingTree
ExtractGatingTree <- function(x, node_paths){
    
    if(length(x@Gating$GatingTreeDF)==0){
        stop('Apply GatingTreeToDF before using this function. \n')
    }
    
    node <- x@Gating$PrunedGatingTreeObject
    pruned_node <- update_nodes_by_paths(node = node, paths = node_paths)
    pruned_node <- prune_tree(pruned_node)
    
    x@Gating$ExtractedGatingTreeObject <- pruned_node
    
    return(x)
    
}
