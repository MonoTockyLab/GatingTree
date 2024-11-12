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
    
    rownames(delta_df) <- NULL
    delta_df <- do.call(rbind, lapply(names(marker_deltas), function(marker) {
        data.frame(marker = marker, delta = marker_deltas[[marker]])
    }))
    
    rownames(delta_df) <- NULL
    
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
#' @param value A character string specifying which type of values to analyze and plot.
#'   Expected values are 'DeltaE' for delta enrichment or 'IE' for interaction effects.
#' @param significance A logical flag indicating whether to highlight significant markers
#'   based on the results of statistical tests.
#' @return Returns the original object `x` for safety.
#' @details The function conducts Kruskal-Wallis tests to determine the significance of
#'   differences across markers, followed by Dunn's test for post-hoc analysis with Bonferroni
#'   correction if significant. It prepares and displays a boxplot of the specified metrics
#'   (Delta Enrichment or Interaction Effects). The markers are ordered and displayed based on
#'   their significance and mean values.
#' @examples
#' \dontrun{
#'   x <- PlotDeltaEnrichment(x, value = "DeltaE", significance = TRUE)
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot labs theme_minimal theme
#' @importFrom dunn.test dunn.test
#' @importFrom rlang sym
#' @family GatingTree


PlotDeltaEnrichment <- function(x, significance = TRUE){
    marker_deltas_df <- x@Gating$GatingTree$delta_df
    kruskal_test_delta <- kruskal.test(marker_deltas_df[['DeltaE_value']] ~ marker_deltas_df[['Marker']])
    
    if (kruskal_test_delta$p.value < 0.05) {
        posthoc_results_delta <- dunn.test(marker_deltas_df$DeltaE_value, marker_deltas_df$Marker, method = "bonferroni", table = FALSE)
    }
    
    mean_delta_values <- aggregate(DeltaE_value ~ Marker, data = marker_deltas_df, FUN = mean)
    colnames(mean_delta_values) <- c("Marker", "MeanDeltaE_value")
    significant_threshold <- 0.05  # Placeholder, adjust based on your Bonferroni correction
    p_vec <- posthoc_results_delta$P.adjusted
    names(p_vec) <- posthoc_results_delta$comparisons
    
    if(significance){
        markers <- names(p_vec) [p_vec < significant_threshold]
    }else{
        markers <- names(p_vec)
    }
    
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for clarity
    plot(p)
    x@Gating$GatingTree$DeltaPlot <- p
    x@Gating$GatingTree$DeltaStats <- kruskal_test_delta
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


add_prune <- function(node,  min_average_proportion = 0.001, theta = 0.02) {
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


find_and_update_nodes <- function(node, path = NULL, res_df) {
    if (is.null(path)) path <- c()

    if (!is.null(node$CurrentPath)) {
        path <- node$CurrentPath  # Correcting path handling
    }
    current_combination <- paste(path[-1], collapse="_")  # Excluding 'rootNode'
    
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
#' @export
#' @family GatingTree


PruneGatingTree <- function(x, max_entropy =0.9, min_enrichment = 0.1, min_average_proportion = 0.001, p_adjust_method = 'BY', theta = 0){
    unique_maxima_df <- x@Gating$GatingTreeDF
    node <- x@Gating$GatingTreeObject
    sampledef <- x@sampledef$sampledef
    
    logic <-(unique_maxima_df$entropy < max_entropy)&( unique_maxima_df$max_enrichment > min_enrichment)
    unique_maxima_df <- unique_maxima_df[logic,]
    
    rootdata <-  node$RootData
    total_cell_per_file <- table(rootdata$file)
    total_cell_per_file <- data.frame(file = names(total_cell_per_file), total = as.vector(total_cell_per_file))
    average_proportion <- p_value_vec <- rep(0, nrow(unique_maxima_df))
    percent_list <- as.list(p_value_vec)
    names(percent_list) <- unique_maxima_df$combination
    
    for(k in 1:nrow(unique_maxima_df)){
        path <- unique_maxima_df$markers_up_to_max[k]
        path <- unlist(strsplit(path, split='_'))
        
        tmp <- findNodeByPath(node, path)
        node_data <- rootdata[rootdata$indices %in% tmp$CurrentNodeIndices,]
        total_cell_num <- table(node_data$file)
        if(nrow(total_cell_num) < nrow(sampledef)){
            p_value_vec = NA
            average_proportion = NA
        }else{
            total_cell_num <- data.frame(file = names(total_cell_num), node_total = as.vector(total_cell_num))
            total_cell_num <- merge(total_cell_num, data.frame( total_cell_per_file), by = 'file')
            total_cell_num <- merge(total_cell_num, sampledef, by = 'file')
            total_cell_num$node_percentage <- total_cell_num$node_total/total_cell_num$total
            percent_list[[k]] <- total_cell_num
            p_value_vec[k] = p_value = wilcox.test(data = total_cell_num, node_percentage~group)$p.value
            average_proportion[k] <- mean(total_cell_num$node_percentage)
        }
        
        
    }
    
    res_df <- cbind(unique_maxima_df, data.frame(p_value = p_value_vec, average_proportion = average_proportion))
    res_df <- res_df[!is.na(res_df$p_value),]
    res_df[['p_adjust']] <- p.adjust(res_df$p_value, method = p_adjust_method)
    
    pruned_node <- add_prune(node, min_average_proportion = min_average_proportion, theta)
    #   pruned_node <- prune_tree(pruned_node)
    
    pruned_node <- find_and_update_nodes(pruned_node, res_df = res_df)
    pruned_node <- prune_tree(pruned_node)
    
    tree_to_df <- function(node){
        enrichments <- collect_all_enrichment(node)
        entropies <- collect_all_entropy(node)
        
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
        
        return(maxima_df)
    }
    
    
    df <- tree_to_df(pruned_node)
    res_df <- res_df[res_df$combination %in% df$combination,]
    
    x@Gating$PrunedGatingTreeObject <- pruned_node
    x@Gating$PrunedGatingTreeDF <- res_df
    
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
#' @param average_proportion A logical flag indicating whether to use the average proportion
#'   of nodes to adjust node sizes and color gradient based on enrichment values.
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
#' @family GatingTree

convert_to_diagrammer <- function(root, size_factor = 1, average_proportion = FALSE, all_labels = TRUE) {
    graph <- create_graph(directed = TRUE)
    node_index <- 1 # To keep track of nodes
    add_nodes_recursively <- function(node) {
        if(node$isRoot){
            simplified_name <- "Root"
            label_text <- ''
            
            node_size <-  size_factor *0.1
            tmp <- 0

            color_value <- rgb(tmp, 0, 1-tmp, alpha = tmp)  # Red gradient

            graph <<- add_node(graph, label = label_text,
                               node_aes = list(title = node$name, fillcolor=color_value, width = node_size, fontcolor = 'black'))
        }
         if (!node$isRoot) {
             if(all_labels){
                formatted_entropy <- round(node$CurrentEntropy, 2)
                if (abs(formatted_entropy) < 0.001) formatted_entropy <- 0
                simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](pos)$", "+", node$name))
                simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](neg)$", "-", simplified_name))
                label_text <- sprintf("%s\nEnrichment: %.2f\nEntropy: %.2f\nAvg. %%: %.1f%%", simplified_name, node$CurrentEnrichment, formatted_entropy, node$AverageProportion * 100)
            }else{
                simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](pos)$", "+", node$name))
                simplified_name <- gsub("(.logdata)|[,]", "", gsub("[.](neg)$", "-", simplified_name))
                label_text <- simplified_name

            }
            if(average_proportion){
                node_size <-  size_factor*node$AverageProportion
                tmp <- node$CurrentEnrichment
                tmp <- max(tmp, 0)
                tmp <- pmin(tmp, 1)
                color_value <- rgb(tmp, 0, 1-tmp, alpha = tmp)  # Red gradient
            }else{
                node_size <-  size_factor*node$CurrentEnrichment
                tmp <- node$CurrentEntropy
                tmp <- max(tmp, 0)
                tmp <- pmin(tmp, 1)
                color_value <- rgb(1-tmp, 0, tmp, alpha = (1-tmp)/1.5)  # Red gradient
            }

            graph <<- add_node(graph, label = label_text,
                               node_aes = list(title = node$name, fillcolor=color_value, width = node_size, fontcolor = 'black'))

            if (!is.null(parent_index)) {
                graph <<- add_edge(graph, from = parent_index, to = node_index)
            }
        }
        current_index <- node_index
        node_index <<- node_index + 1
        if (!node$isLeaf) {
            for (child in node$children) {
                parent_index <<- current_index
                add_nodes_recursively(child)
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

collect_all_enrichment <- function(tree) {
  enrichment_list <- list()

  # Collect enrichment from the current node if it exists
  if (!is.null(tree$History) && !is.null(tree$History$enrichment)) {
    enrichment_list <- c(enrichment_list, list(tree$History$enrichment))
  }

  # Recursively collect enrichment from child nodes
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

  # Recursively collect $CurrentNodeName from child nodes that have a $Peak
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


collect_all_entropy <- function(tree) {
  entropy_list <- list()

  # Collect entropy from the current node if it exists
  if (!is.null(tree$History) && !is.null(tree$History$entropy)) {
    entropy_list <- c(entropy_list, list(tree$History$entropy))
  }

  # Recursively collect entropy from child nodes
  if (length(tree$Children) > 0) {
    for (child_name in names(tree$Children)) {
      child_entropys <- collect_all_entropy(tree$Children[[child_name]])
      entropy_list <- c(entropy_list, child_entropys)
    }
  }

  return(entropy_list)
}






