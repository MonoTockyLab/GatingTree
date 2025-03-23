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

#' Calculate Entropy
#'
#' This function calculates the entropy of a probability distribution.
#'
#' @param probs A numeric vector representing the probability distribution.
#'
#' @return The calculated entropy value.
#' @keywords internal
#' @examples
#' \dontrun{
#' entropy_result <- calculate_entropy(probs)
#'}
#' @family GatingTree
#' @noRd
calculate_entropy <- function(probs) {
    -sum(probs * log2(probs + 1e-10))
}


#' Calculate Enrichment Between Two Groups
#'
#' This function calculates the log2-ratio of the sum of percentages normalized by group size
#' between an experimental group and a control group. It adds a small constant to avoid
#' division by zero.
#'
#' @param percentage_data A character to specify the column name for the percentage data.
#' @param expr_group The name of the experimental group as a single string which must
#'   match exactly with one of the groups in `sampledef`.
#' @param ctrl_group The name of the control group as a single string which must
#'   match exactly with one of the groups in `sampledef`.
#'
#' @return A single numeric value representing the log2-transformed ratio of
#'   normalized sums of percentages between the experimental and control groups.
#'
#' @examples
#' \dontrun{
#' percentage <- c(0.1, 0.2, 0.15, 0.05)
#' calculate_enrichment(df, "exp", "ctrl")
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd
calculate_enrichment <- function(df, percentage_data = 'node_percentage', expr_group = NULL, ctrl_group = NULL) {
    percentage <- df[,percentage_data]
    group1_logic <- df[['group']]==expr_group
    group2_logic <- df[['group']]==ctrl_group
    p_g1 <- sum(percentage[group1_logic]/length(percentage[group1_logic]))
    p_g2 <-sum(percentage[group2_logic]/length(percentage[group2_logic]))
    output <- log2((p_g1 + 1e-10)/ (p_g2 + 1e-10))
    return(output)
}


#' Calculate Baseline Entropy
#'
#' This function calculates the baseline entropy of the groups in the data.
#'
#' @param data A data frame containing group information.
#'
#' @return The calculated baseline entropy value.
#' @keywords internal
#' @examples
#' \dontrun{
#' baseline_entropy(sampledef)
#'}
#' @family GatingTree
#' @noRd
baseline_entropy <- function(data) {
    freqs <- table(data[['group']])
    probs <- freqs / sum(freqs)
    calculate_entropy(probs)
}


#' Apply Gating Conditions
#'
#' This function applies a set of gating conditions to flow cytometry data.
#'
#' @param data A data frame or matrix containing flow cytometry data.
#' @param combinations_matrix A matrix specifying the combinations of gating conditions.
#' @param thresholds A vector of thresholds for the gating conditions.
#'
#' @return A matrix indicating whether each cell passes the gating conditions.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @keywords internal
#' @examples
#' \dontrun{
#'     gating_results <- apply_gating_conditions(marker_data, combinations_matrix, thresholds)
#'}
#' @family GatingTree
#' @noRd
apply_gating_conditions <- function(data, combinations_matrix, thresholds) {
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }
    cat("Analyzing combinatorial gates...\n")
    thresholds_matrix <- matrix(rep(thresholds, each = nrow(data)), nrow = nrow(data), byrow = TRUE)

    positive_masks <- data > thresholds_matrix
    negative_masks <- data <= thresholds_matrix
    
    results <- matrix(FALSE, nrow = nrow(data), ncol = nrow(combinations_matrix))
    pb <- txtProgressBar(min = 0, max = nrow(combinations_matrix), style = 3)
    on.exit(close(pb))
    
    for (i in seq_len(nrow(combinations_matrix))) {
        condition <- combinations_matrix[i, ]
        
        current_mask <- matrix(TRUE, nrow = nrow(data), ncol = length(condition))
        for (j in seq_along(condition)) {
            if (condition[j] == 2) {
                current_mask[, j] <- positive_masks[, j]
            } else if (condition[j] == 1) {
                current_mask[, j] <- negative_masks[, j]
            } else {
                current_mask[, j] <- TRUE
            }
        }
        results[, i] <- rowSums(current_mask, na.rm = TRUE) == ncol(current_mask)
        setTxtProgressBar(pb, i)
    }
    
    return(results)
}


#' Generate Concatenated Marker Names Based on States
#'
#' Constructs a string that represents a combination of markers and their states
#' (negative or positive). This function appends ".neg" or ".pos" to each marker
#' based on the state provided in `marker_states`. The function returns a single
#' string with these names concatenated, separated by underscores. If no markers
#' are specified as included (non-empty), it returns "all_unassigned".
#'
#' @param marker_states A numeric vector indicating the state of each marker where
#'   1 represents negative (".neg") and 2 represents positive (".pos").
#' @param markers A character vector of marker names corresponding to the states
#'   in `marker_states`.
#'
#' @return A character string of concatenated marker names with their states.
#'   If no markers are specified, returns "all_unassigned".
#'
#' @examples
#' \dontrun{
#' markers <- c("CD4", "CD8", "CD19")
#' states <- c(1, 2, 0)  # CD4.neg, CD8.pos, CD19 not included
#' generate_marker_names(states, markers)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd
generate_marker_names <- function(marker_states, markers) {
    state_labels <- c(".neg", ".pos")  # Labels for states
    names <- character(length(markers))
    for (i in seq_along(markers)) {
        if (marker_states[i] == 1) {
            names[i] <- paste0(markers[i], state_labels[1])  # Appends ".neg"
        } else if (marker_states[i] == 2) {
            names[i] <- paste0(markers[i], state_labels[2])  # Appends ".pos"
        }
    }
    names <- names[names != ""]  # Remove empty elements for markers not included
    names <- paste(names, collapse="_")
    if(names==""|is.na(names)) names <- "all_unassigned"
    return(names)
}

#' General Node Rule for Gating Strategy
#'
#' This function processes gating decisions based on given markers and their states,
#' managing data from a root node to determine optimal gating paths based on entropy and
#' enrichment calculations. It handles decision-making at each node of a gating tree,
#' determining whether to continue subdividing or to terminate based on statistical
#' thresholds and data availability.
#'
#' @param currentNode A list representing the current node in the gating tree, including
#'   markers used, current node indices, and gating history.
#' @param root_data A data frame containing the complete dataset for analysis.
#' @param sampledef A data frame specifying sample definitions and group assignments.
#' @param neg_gate A list containing thresholds for negative gating decisions.
#' @param expr_group The name of the experimental group within `sampledef`.
#' @param ctrl_group The name of the control group within `sampledef`.
#' @param total_cell_per_file A data frame mapping file names to total cell counts per file.
#' @param usedmarkers A vector of markers that have already been used in previous steps of gating.
#' @param min_cell_num The minimal number of cells allowed in nodes.
#'
#' @return A list structure describing the outcomes at the current node, including any child
#'   nodes created, or indicators if the node processing leads to termination.
#'
#' @examples
#' # Assuming the function is part of a larger framework and dependencies are met:
#' @keywords internal
#' @examples
#' \dontrun{
#' result <- general_node_rule(currentNode, markers, root_data, sampledef, neg_gate, 
#' expr_group, ctrl_group, total_cell_per_file, usedmarkers)
#'}
#'
#' @importFrom dplyr count full_join mutate coalesce
#' @importFrom rlang sym
#' @family GatingTree
#' @noRd
general_node_rule <- function(currentNode, root_data, sampledef, neg_gate, expr_group, ctrl_group, total_cell_per_file, usedmarkers, min_cell_num = 25) {
    
    if(currentNode$Leaf){
        return(currentNode)
    }
    all_indices <- root_data$indices
    history <- currentNode$History
    available_markers <- currentNode$AvailableMarkers
    current_indices <- currentNode$CurrentNodeIndices
    current_indices_logic <- all_indices %in% current_indices
    
    parent_marker_state <- currentNode$CurrentMarkerState
    node_data <- root_data[current_indices_logic, c(available_markers, 'file')]

    group_levels <- unique(sampledef[['group']])
    initial_entropy <- baseline_entropy(sampledef)
    thresholds <- setNames(neg_gate$negative.gate, neg_gate$variable)
    
    node_marker_names <- colnames(node_data)[colnames(node_data) != 'file']

    aligned_thresholds <- thresholds[node_marker_names]
    is_positive <- sweep(node_data[,node_marker_names], 2, aligned_thresholds, ">")
    is_negative <- !is_positive
    
   # total_cell_num <- table(node_data$file)
   # total_cell_num <- data.frame(file = names(total_cell_num), node_total = as.vector(total_cell_num))
   # total_cell_num <- merge(total_cell_num, data.frame( total_cell_per_file), by = 'file')
   # total_cell_num$node_percentage <- total_cell_num$node_total/total_cell_num$total
   total_cell_num <- node_data %>%
     count(file, name = "node_total") %>%
     full_join(data.frame(total_cell_per_file), by = "file") %>%
     full_join(sampledef, by = "file") %>%
     mutate(node_total = coalesce(!!sym("node_total"), 0),  # Replace NA values with 0
            node_percentage = !!sym("node_total") / !!sym("total"))
            
    min_node_total <- min(total_cell_num$node_total)
    
    if(min_node_total < min_cell_num){
        currentNode$Leaf <- TRUE
        currentNode$Exhausted <- TRUE
        currentNode$Terminated  <- TRUE

        return(currentNode)
        
    }else{
        initial_enrichment <- calculate_enrichment(total_cell_num, percentage_data = 'node_percentage', expr_group = expr_group, ctrl_group = ctrl_group)
        if(is.na(initial_enrichment)){
            initial_enrichment <- 0
        }
        negative_entropy_scores <- positive_entropy_scores <- rep(1, length(available_markers))
        positive_enrichment_scores <- negative_enrichment_scores <- rep(0, length(available_markers))
        names(negative_entropy_scores) <-names(positive_entropy_scores) <- names(positive_enrichment_scores) <-  names(negative_enrichment_scores) <- available_markers
        
        positive_average_proportion <- negative_average_proportion <- positive_indices <- negative_indices <- as.list(rep(0, length(available_markers)))
        names(positive_average_proportion) <- names(negative_average_proportion)  <- names(positive_indices) <-names(negative_indices) <-  available_markers
        
        
        for (marker in available_markers) {

            marker_response <- ifelse(is_positive[, marker], 'Positive', 'Negative')
            positive_indices[[marker]] <- all_indices[current_indices_logic][is_positive[, marker]]
            negative_indices[[marker]] <- all_indices[current_indices_logic][is_negative[, marker]]
            file_response <- table(node_data$file, marker_response)
            percentage_data <- as.data.frame(file_response); colnames(percentage_data) <- c('file', 'marker_response', 'gated_cell_num')
            percentage_data <- merge(percentage_data, sampledef, by = 'file'); percentage_data <- merge(percentage_data, total_cell_per_file, by = 'file')
            percentage_data$gated_percentage <- percentage_data$gated_cell_num / percentage_data$total
            positive_percentage <- percentage_data[percentage_data$marker_response =='Positive',]
            negative_percentage <- percentage_data[percentage_data$marker_response =='Negative',]
            
            if ((nrow(positive_percentage) == 0)|(nrow(negative_percentage) == 0)){ #| #(nrow(positive_percentage) < nrow(sampledef)) |  (nrow(negative_percentage) < nrow(sampledef)) ) {
                next
            }
            
            positive_group_mean_percentages <- aggregate(gated_percentage ~ group, data = positive_percentage, FUN = mean)
            negative_group_mean_percentages <- aggregate(gated_percentage ~ group, data = negative_percentage, FUN = mean)
            positive_enrichment <- positive_enrichment_scores[marker] <- calculate_enrichment(positive_percentage, percentage_data = 'gated_percentage', expr_group = expr_group, ctrl_group = ctrl_group)
            negative_enrichment <- negative_enrichment_scores[marker] <- calculate_enrichment(negative_percentage, percentage_data = 'gated_percentage', expr_group = expr_group, ctrl_group = ctrl_group)
            
            if((positive_enrichment <= initial_enrichment)&&(negative_enrichment <= initial_enrichment) ){
                next
            }
            
            if (is.null(positive_group_mean_percentages) || nrow(positive_group_mean_percentages) == 0 || is.null(negative_group_mean_percentages) || nrow(negative_group_mean_percentages) == 0) {
                next
            }
            positive_entropy_scores[marker] <- gating_entropy(positive_percentage)
            negative_entropy_scores[marker] <- gating_entropy(negative_percentage)
            positive_average_proportion[[marker]] <- mean(positive_percentage$gated_percentage)#mean(total_cell_num$node_percentage)
            negative_average_proportion[[marker]] <- mean(negative_percentage$gated_percentage)#mean(total_cell_num$node_percentage)

        }
        
        enrichment_df <- data.frame(positive_enrichment_scores = positive_enrichment_scores, negative_enrichment_scores = negative_enrichment_scores); enrichment_df <- t(enrichment_df)
        entropy_df <- data.frame(positive_entropy_scores = positive_entropy_scores, negative_entropy_scores = negative_entropy_scores); entropy_df <- t(entropy_df)
        colnames(entropy_df) <- colnames(enrichment_df)  <- available_markers
        
        currentNode$NodeEnrichmentDf <- enrichment_df
        currentNode$NodeEntropyDf <- entropy_df
        
        if(is.null(positive_enrichment_scores)|is.null(negative_enrichment_scores)){
            currentNode$Leaf <- TRUE
            currentNode$Terminated <- TRUE
            cat(paste('Enrichment score could not be calculated: for', paste(currentNode$CurrentPath, collapse='_'),'\n'))
            return(currentNode)
        }else{
            
            if (all(positive_enrichment_scores < initial_enrichment) && (all(negative_enrichment_scores < initial_enrichment))) {
                currentNode$Leaf <- TRUE
                currentNode$Terminated <- TRUE
                currentNode$Peaked <- TRUE
                return(currentNode)
            }
        }
        
        positive_enrichment_scores <- positive_enrichment_scores[!is.na(positive_enrichment_scores)]
        negative_enrichment_scores <- negative_enrichment_scores[!is.na(negative_enrichment_scores)]
        
        positive_rule <- (positive_enrichment_scores > max(initial_enrichment, 0)) & (positive_entropy_scores < initial_entropy)
        negative_rule <- (negative_enrichment_scores >  max(initial_enrichment, 0)) & (negative_entropy_scores < initial_entropy)
        positive_candidates <- available_markers[positive_rule]
        negative_candidates <- available_markers[negative_rule]
        
        children <- list()
        
        if(length(positive_candidates)!=0){
            positive_childNode <- as.list(positive_candidates)
            names(positive_childNode) <- positive_candidates
  
            for(positive_candidate in positive_candidates){
                positive_history <- history
                positive_history$entropy <- c(positive_history$entropy, positive_entropy_scores[positive_candidate])
                names(positive_history$entropy)[length(positive_history$entropy)]<- paste(positive_candidate, 'pos', sep='.')
                positive_history$enrichment <- c(positive_history$enrichment, positive_enrichment_scores[positive_candidate])
                names(positive_history$enrichment)[length(positive_history$enrichment)]<- paste(positive_candidate, 'pos', sep='.')
                usedmarkers <- c(currentNode$UsedMarkers, positive_candidate)
                childAvailableMarkers <- setdiff(available_markers, positive_candidate)
                current_entropy <- positive_entropy_scores[positive_candidate]
                current_enrichment <-  positive_enrichment_scores[positive_candidate]
                positive_current_marker_state <- parent_marker_state
                positive_current_marker_state[positive_candidate] <- 2
                
                positive_childNode[[positive_candidate]]  <- createChildNode(marker = positive_candidate,
                current_marker_state = positive_current_marker_state, indices = positive_indices[[positive_candidate]],
                available_markers = childAvailableMarkers,current_entropy, current_enrichment = current_enrichment,
                entropy_scores = entropy_df, enrichment_scores = enrichment_df, average_proportion = positive_average_proportion[[positive_candidate]], history = positive_history,
                isPositive = TRUE, depth= currentNode$Depth + 1, usedmarkers = usedmarkers, path = currentNode$CurrentPath)
                
                
            }
            names(positive_childNode) <- paste(positive_candidates, '.pos', sep='')
  
            children <- positive_childNode
            
        }
        if(length(negative_candidates)!=0){
            negative_childNode <- as.list(negative_candidates)
            names(negative_childNode) <- negative_candidates
            
            for(negative_candidate in negative_candidates){
                negative_history <- history
                negative_history$entropy <- c(negative_history$entropy, negative_entropy_scores[negative_candidate])
                names(negative_history$entropy)[length(negative_history$entropy)]<- paste(negative_candidate, 'neg', sep='.')
                negative_history$enrichment <- c(negative_history$enrichment, negative_enrichment_scores[negative_candidate])
                names(negative_history$enrichment)[length(negative_history$enrichment)]<- paste(negative_candidate, 'neg', sep='.')
                current_entropy <- negative_entropy_scores[negative_candidate]
                current_enrichment <-  negative_enrichment_scores[negative_candidate]
                negative_current_marker_state <- parent_marker_state
                negative_current_marker_state[negative_candidate] <- 1
                usedmarkers <- c(currentNode$UsedMarkers, negative_candidate)
                childAvailableMarkers <- setdiff(available_markers, negative_candidate)
                
                negative_childNode[[negative_candidate]] <- createChildNode(marker = negative_candidate, current_marker_state = negative_current_marker_state,
                indices = negative_indices[[negative_candidate]], average_proportion = negative_average_proportion[[negative_candidate]], available_markers = childAvailableMarkers, current_entropy,
                current_enrichment = current_enrichment, entropy_scores = entropy_df, enrichment_scores = enrichment_df, history = negative_history,
                isPositive = FALSE, depth= currentNode$Depth + 1, usedmarkers = usedmarkers, path = currentNode$CurrentPath)
            }
            names(negative_childNode) <- paste(negative_candidates, '.neg', sep='')
            children <- c(children, negative_childNode)
            
        }
        
        if(length(positive_candidates)==0 && length(negative_candidates)==0){
            
            currentNode$Leaf <- TRUE
            currentNode$Terminated  <-  TRUE
            currentNode$PercentageData  <-  total_cell_num
            currentNode$AverageProportion  <-  mean(total_cell_num$node_percentage)

            return(currentNode)
        }
        currentNode$Children <- children
        return(currentNode)
    }
}




#' Calculate Gating Entropy Based on Positive Percentages
#'
#' Computes the entropy associated with gating, categorizing the data into 'High' or 'Low'
#' based on whether the gated percentage of a group is above or below the overall mean percentage.
#' This function is useful for analyzing the variability or predictability of gating outcomes
#' across different groups.
#'
#' @param positive_percentage A data frame with at least two columns: 'group' which
#'   categorizes the data, and 'gated_percentage' which represents the percentage
#'   of some criterion being met (e.g., cells positive for a marker).
#'
#' @return A single numeric value representing the entropy, which quantifies the
#'   unpredictability associated with the gating classifications across the groups.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' gating_entropy(data)
#'}
#' @family GatingTree
#' @noRd

gating_entropy <- function(positive_percentage){
    group_levels <- unique(positive_percentage$group)
    positive_group_mean_percentages <- aggregate(gated_percentage ~ group, data = positive_percentage, FUN = mean)
    mean_percentage <- sum(positive_group_mean_percentages$gated_percentage)/length(group_levels)
    class_response <- ifelse(positive_percentage$gated_percentage > mean_percentage, 'High', 'Low')
    table_combined <- table(class_response, positive_percentage$group)
    probs_combined <- prop.table(table_combined, 1)  # Row-wise proportions
    row_entropies <- -rowSums(probs_combined * log2(probs_combined + 1e-10))
    
    output <- sum(row_entropies * (rowSums(table_combined) / sum(table_combined)))
    output <- round(output, digits = 3)
    return(output)
}


#' Create a Gating Tree Object
#'
#' Initializes and populates a gating tree based on user-defined criteria,
#' handling decision-making through multiple layers of a gating hierarchy.
#' This function integrates gating rules, applies negative gating definitions,
#' and assesses cell population statistics to manage the flow cytometry data analysis process.
#'
#' @param x An object, expected to be of class 'FlowObject',
#'   containing the initial data and parameters for gating.
#' @param select_markers Logical; if TRUE, allows the user to select markers interactively.
#' @param graphics Logical; if TRUE, enables graphical selection of markers.
#' @param markers Optional; a vector of markers to be included if not choosing interactively.
#' @param expr_group Optional; a character to specify the experimental group. If NULL, this is determined interactively.
#' @param ctrl_group Optional; a character to specify the control group. If NULL, this is determined interactively.
#' @param maxDepth Integer; the maximum depth of the gating tree, controlling how many
#'   levels of decision nodes can be created.
#' @param min_cell_num The minimal number of cells allowed in nodes.
#' @param verbose Logical indicating whether to print progress messages and outputs. Default is \code{TRUE}.
#'
#' @return Modifies the input object by adding a 'GatingTreeObject' that contains
#'   the entire structure of gating decisions and nodes.
#'
#' @details
#' The function first checks for the type of the input object to ensure it matches
#' expected classes. It then extracts necessary data and parameters from the object,
#' such as negative gating thresholds and sample definitions. Depending on the options,
#' it may allow interactive selection of markers. The function constructs a hierarchical
#' tree where each node represents a gating decision based on statistical calculations
#' like entropy and enrichment, which are used to determine the next steps in the gating
#' process or to terminate the process.
#'
#' The gating process involves:
#' - Merging data with sample definitions.
#' - Calculating initial gating statistics.
#' - Recursively creating child nodes based on gating outcomes and thresholds.
#' - Dynamically managing markers and gating paths based on user-defined depth and available data.
#'
#' @examples
#' \dontrun{
#' # Assuming 'x' is properly instantiated and contains necessary gating setups:
#' x <- createGatingTreeObject(x, select_markers = TRUE, graphics = FALSE, maxDepth = 3)
#'}
#' @export
#' @family GatingTree

createGatingTreeObject <- function(x, select_markers = FALSE, graphics = FALSE, markers = NULL, maxDepth = 3, min_cell_num = 25, expr_group = NULL, ctrl_group = NULL, verbose = TRUE) {
    if(!inherits(x, "FlowObject")){
        stop("Use a FlowObject.")
    }
    if(verbose){
        progress <- txtProgressBar(min = 0, max = 100, style = 3)
    }
    
    
    X <- x@Data
    neg_gate <- x@QCdata$negative_gate_def
    sampledef <- x@sampledef$sampledef
    group_levels <- unique(sampledef[['group']])
    X <- merge(X, sampledef, by = 'file')
    
    if(is.null(x@sampledef$expr_group)&is.null(x@sampledef$ctrl_group)&(is.null(expr_group)|is.null(ctrl_group))){
        expr_group <- select.list(group_levels, graphics = FALSE, title = "Which group is the target population (experimental group)?", multiple =FALSE)
        cat(paste(expr_group, 'has been chosen as the target population. \n'))
        ctrl_group <- select.list(group_levels, graphics = FALSE, title = "Which group is the reference population (control group)?", multiple =FALSE)
        cat(paste(ctrl_group, 'has been chosen as the reference population. \n'))
    }else{
        if(is.null(expr_group)|is.null(ctrl_group)){
            expr_group <- x@sampledef$expr_group
            ctrl_group <- x@sampledef$ctrl_group
        }

    }
    
    lg <- grepl(pattern = 'logdata', colnames(X))
    neg_gate <- neg_gate[neg_gate$negative.gate !=0, ]
    choices <- cn <- neg_gate$variable
    cn <- c('file',intersect(colnames(X)[lg], cn))
    
    if(is.null(markers)){
        if(select_markers){
            markers <- select.list(choices, graphics = graphics, title = "Data to be modelled", multiple =TRUE)
              choices <- markers
        }else{#Unless select_markers = TRUE, neg_gate$variable is used
            markers <- choices
        }
    }else{#If markers is specified and it is a character vector, then the existing variables specified by the option marker are used
        if(is.character(markers)){
            markers <- intersect(choices, markers)
        }else{
            markers <- choices
        }
        
    }
    if(verbose){
        cat("The number of markers used: ", length(markers), "\n")
        
    }
    available_markers <-  markers
    all_indices <- 1:nrow(X)
    root_data <- X[, c(available_markers, 'file', 'group')]
    root_data$indices <- all_indices
    initial_entropy <- baseline_entropy(sampledef)

    thresholds <- setNames(neg_gate$negative.gate, neg_gate$variable)
    aligned_thresholds <- thresholds[available_markers]
    is_positive <- sweep(root_data[,available_markers], 2, aligned_thresholds, ">")
    is_negative <- !is_positive
    
    total_cell_num <- table(root_data$file)
    total_cell_num <- data.frame(file = names(total_cell_num), total = as.vector(total_cell_num))
    colnames(total_cell_num) <- c('file', 'total')
    total_cell_per_file <- total_cell_num
    initial_enrichment <- 0 #by definition
    
    negative_entropy_scores <- positive_entropy_scores <- rep(1, length(available_markers))
    positive_indices <- negative_indices <- as.list(rep(0, length(available_markers)))
    positive_average_proportions <- negative_average_proportions <- positive_enrichment_scores <- negative_enrichment_scores <- rep(0, length(available_markers))
    names(positive_average_proportions) <- names(negative_average_proportions) <- names(positive_indices) <-names(negative_indices) <- names(negative_entropy_scores) <-names(positive_entropy_scores) <- names(positive_enrichment_scores) <-  names(negative_enrichment_scores) <- available_markers
    
    for (marker in available_markers) {

        marker_response <- ifelse(is_positive[, marker], 'Positive', 'Negative')
        positive_indices[[marker]] <- all_indices[is_positive[, marker]]
        negative_indices[[marker]] <- all_indices[is_negative[, marker]]
        file_response <- table(root_data$file, marker_response)
        percentage_data <- as.data.frame(file_response)
        
        colnames(percentage_data) <- c('file', 'marker_response', 'gated_cell_num')
        percentage_data <- merge(percentage_data, sampledef, by = 'file')
        percentage_data <- merge(percentage_data, total_cell_num, by = 'file')
        percentage_data$gated_percentage <- percentage_data$gated_cell_num / percentage_data$total
        positive_percentage <- percentage_data[percentage_data$marker_response =='Positive',]
        negative_percentage <- percentage_data[percentage_data$marker_response =='Negative',]


        if ((nrow(positive_percentage) == 0)|(nrow(negative_percentage) == 0)
       # | (nrow(positive_percentage) < nrow(sampledef)) |  (nrow(negative_percentage) < nrow(sampledef))
        ) {
            
            next
        }
        
        positive_group_mean_percentages <- aggregate(gated_percentage ~ group, data = positive_percentage, FUN = mean)
        negative_group_mean_percentages <- aggregate(gated_percentage ~ group, data = negative_percentage, FUN = mean)

        positive_enrichment_scores[marker] <- calculate_enrichment(positive_percentage, percentage_data = 'gated_percentage', expr_group = expr_group, ctrl_group = ctrl_group)
        negative_enrichment_scores[marker] <- calculate_enrichment(negative_percentage, percentage_data = 'gated_percentage', expr_group = expr_group, ctrl_group = ctrl_group)

        positive_average_proportions[marker] <-  mean(positive_group_mean_percentages$gated_percentage)
        negative_average_proportions[marker] <- mean(negative_group_mean_percentages$gated_percentage)

        if (is.null(positive_group_mean_percentages) || nrow(positive_group_mean_percentages) == 0 || is.null(negative_group_mean_percentages) || nrow(negative_group_mean_percentages) == 0) {

            next
        }
        positive_entropy_scores[marker] <- gating_entropy(positive_percentage)
        negative_entropy_scores[marker] <- gating_entropy(negative_percentage)
        
        
    }
    
    enrichment_df <- data.frame(
    positive_enrichment_scores = positive_enrichment_scores,
    negative_enrichment_scores = negative_enrichment_scores
    )

    entropy_df <- data.frame(
    positive_entropy_scores = positive_entropy_scores,
    negative_entropy_scores = negative_entropy_scores
    )
    
    enrichment_df <- t(enrichment_df)
    entropy_df <- t(entropy_df)
    marker_state <- rep(0, length(markers))
    names(marker_state) <- available_markers
    
    colnames(entropy_df) <- colnames(enrichment_df)  <- available_markers
    history <- list(entropy = c(), enrichment = c())
    rootNode <- structure(list(
    RootData = root_data,
    CurrentNodeName = "rootNode",
    CurrentNodeGate = "Unassigned",
    CurrentPath = "rootNode",
    CurrentMarkerState = marker_state,
    CurrentNodeIndices = all_indices,
    RootNodeEnrichmentVector = enrichment_df,
    RootNodeEntropyVector = entropy_df,
    CurrentEntropy = initial_entropy,
    CurrentEnrichment = initial_enrichment,
    AverageProportion = 1,
    History = history,
    Children = list(),
    AvailableMarkers = available_markers,
    Depth = 0,
    Leaf = FALSE,
    Terminated = FALSE,
    UsedMarkers = c()
    ), class = "GatingTreeRootNode")

    for(marker in markers){

        positive_marker <- paste(marker, 'pos', sep='.')
        negative_marker <- paste(marker, 'neg', sep='.')
        positive_history <- list(
        entropy = positive_entropy_scores[marker],
        enrichment = positive_enrichment_scores[marker]
        )
        names(positive_history$entropy)[1]<- positive_marker
        names(positive_history$enrichment)[1]<- positive_marker

        negative_history <- list(
        entropy = negative_entropy_scores[marker],
        enrichment = negative_enrichment_scores[marker]
        )
        names(negative_history$entropy)[1]<- negative_marker
        names(negative_history$enrichment)[1]<- negative_marker
        usedmarkers <- marker
        
        current_entropy <- positive_entropy_scores[marker]
        current_enrichment <- positive_enrichment_scores[marker]
        negative_current_marker_state <- positive_current_marker_state <- marker_state
        positive_current_marker_state[marker] <- 2
        negative_current_marker_state[marker] <- 1
        
        
        positive_childNode <- createChildNode(marker, indices = positive_indices[[marker]], current_marker_state = positive_current_marker_state, available_markers = markers, current_entropy = current_entropy, current_enrichment = current_enrichment, entropy_scores = entropy_df, enrichment_scores = enrichment_df, average_proportion = positive_average_proportions[marker], history=positive_history, isPositive = TRUE, depth = 1, usedmarkers = usedmarkers, path = rootNode$CurrentPath)

        current_entropy_neg <- negative_entropy_scores[marker]
        current_enrichment_neg <- negative_enrichment_scores[marker]
        negative_average_proportion <- mean(negative_group_mean_percentages$gated_percentage)
        negative_childNode <- createChildNode(marker, indices = negative_indices[[marker]], current_marker_state = negative_current_marker_state, available_markers = markers, current_entropy = current_entropy_neg, current_enrichment = current_enrichment_neg, entropy_scores = entropy_df, enrichment_scores = enrichment_df, average_proportion = negative_average_proportions[marker], history=negative_history, isPositive = FALSE, depth = 1, usedmarkers = usedmarkers, path = rootNode$CurrentPath)
        
        rootNode <- addChildNode(rootNode, positive_childNode, path = 'rootNode')
        rootNode <- addChildNode(rootNode, negative_childNode, path = 'rootNode')
        
    }
    
    addChildNodeDynamically <- function(rootNode, childNode, parentPath) {
        parentNode <- findNodeByPath(rootNode, parentPath)
        parentNode$Children[[childNode$CurrentNodeName]] <- childNode
    }
    
    firstLevelNames <- names(rootNode$Children)
    
    if(verbose){
        progress_increment <- round(100/length(firstLevelNames))
        current_progress <- 0
        }

    
    for(firstLevelName in firstLevelNames){
        rootNode$Children[[firstLevelName]] <- recursiveAddChildNode(rootNode$Children[[firstLevelName]], root_data = root_data, sampledef, neg_gate= neg_gate, expr_group=expr_group, ctrl_group=ctrl_group, total_cell_per_file, maxDepth = maxDepth, usedmarkers, depth = 1, min_cell_num = min_cell_num)
        
        if(verbose){
            current_progress <- current_progress + progress_increment
            setTxtProgressBar(progress, current_progress)
        }
        
    }
    
    if(verbose){
        setTxtProgressBar(progress, 100)
        close(progress)
    }

    x@sampledef$expr_group <- expr_group
    x@sampledef$ctrl_group <- ctrl_group
    x@Gating[['GatingTreeObject']] <- rootNode
    
    return(x)
}


#' Find a Node by Path in a Gating Tree
#'
#' Traverses a gating tree starting from a specified root node to find and return a node
#' located at a given path. The path should be a sequence of node names indicating the
#' traversal route from the root to the target node.
#'
#' @param rootNode The root node of the gating tree from which to start the search.
#' @param path A character vector representing the path to the desired node. Each element
#'   of the vector should correspond to a node name at each level of the tree.
#'
#' @return Returns the node at the specified path if found.
#'
#' @examples
#' \dontrun{
#' rootNode <- createGatingTreeObject(...) # setup initial node, assumes function definition
#' path <- c("rootNode", "CD4pos", "CD8neg") # example path to find a specific node
#' tryCatch({
#'   targetNode <- findNodeByPath(rootNode, path)
#'   print(targetNode)
#' }, error = function(e) {
#'   cat("Error in findNodeByPath: ", e$message, "\n")
#' })
#' }
#'
#' @export
#' @family GatingTree


findNodeByPath <- function(rootNode, path) {
    currentNode <- rootNode
    #path <- path[-1]
    for (i in 1:length(path)) {
        pathPart <- path[i]
        if (!is.null(currentNode$Children[[pathPart]])) {
            currentNode <- currentNode$Children[[pathPart]]
        } else {
            stop("Path not found: ", paste(path, collapse = " -> "))
        }
    }
    return(currentNode)
}

#' Recursively Add Child Nodes in a Gating Tree
#'
#' This function recursively expands a gating tree by adding child nodes according to gating rules.
#' It is designed to iteratively apply gating decisions down to a specified depth or until no further
#' subdivisions are applicable (terminal nodes).
#'
#' @param currentNode The current node in the gating tree from which children will be generated.
#' @param root_data A data frame containing the data set used for gating decisions.
#' @param sampledef A data frame specifying sample definitions and group assignments.
#' @param neg_gate A list containing thresholds for negative gating decisions.
#' @param expr_group The name of the experimental group within `sampledef`.
#' @param ctrl_group The name of the control group within `sampledef`.
#' @param total_cell_per_file A data frame mapping file names to total cell counts, used for normalization.
#' @param maxDepth The maximum depth to which the tree can expand.
#' @param usedmarkers A vector of markers already used in the gating path up to the current node.
#' @param min_cell_num The minimal number of cells allowed in nodes.
#' @param depth Integer, the depth of the node in the tree.
#'
#' @return The modified current node with potentially new child nodes added, reflecting the gating tree expansion.
#'
#' @details
#' The function checks if the current node is terminated or if its depth equals the maximum allowed depth.
#' If not, it applies gating rules to decide how to expand the tree by adding child nodes. These decisions
#' are based on statistical measures such as enrichment and entropy, calculated for different marker states.
#' If a child node results in a terminal condition ('Leaf'), the tree expansion stops at that node.
#'
#' The function uses recursion to navigate and expand deeper levels of the tree, ensuring that all potential
#' gating paths are explored up to the `maxDepth` or until no further divisions are valid.
#'
#' @examples
#' \dontrun{
#' # Assuming currentNode and other required objects are predefined:
#' updatedNode <- recursiveAddChildNode(currentNode, root_data, sampledef,
#' neg_gate, expr_group, ctrl_group, total_cell_per_file, maxDepth = 3, usedmarkers)
#' }
#' @keywords internal
#' @family GatingTree
#' @noRd
recursiveAddChildNode <- function(currentNode, root_data, sampledef, neg_gate, expr_group, ctrl_group, total_cell_per_file, maxDepth = 3, usedmarkers, min_cell_num=25, depth = 1) {

    if (currentNode$Leaf |currentNode$Depth >= maxDepth) {
        return(currentNode)
    }
    
    currentNode <- general_node_rule(currentNode, root_data, sampledef, neg_gate, expr_group, ctrl_group, total_cell_per_file, usedmarkers, min_cell_num)
    
    if(length(currentNode$Children)!=0){
        for(child_name in names(currentNode$Children)){
            currentNode$Children[[child_name]] <- recursiveAddChildNode(currentNode = currentNode$Children[[child_name]], root_data, sampledef, neg_gate, expr_group, ctrl_group, total_cell_per_file, maxDepth = maxDepth, depth + 1, usedmarkers, min_cell_num=25)
            
        }
        
        
    }else{
        return(currentNode)
    }
    
    return(currentNode)
}

#' Create a Child Node in a Gating Tree
#'
#' Constructs a child node for a gating tree based on the specified gating marker
#' and its state, along with related gating statistics and history.
#'
#' @param marker Character string of the marker name.
#' @param current_marker_state Current state of all markers at the node.
#' @param indices Indices of the data used in this node.
#' @param available_markers Markers available for further gating.
#' @param current_entropy Current entropy score of the node.
#' @param current_enrichment Current enrichment score of the node.
#' @param average_proportion Average proportion of cells within the node.
#' @param entropy_scores Entropy scores dataframe.
#' @param enrichment_scores Enrichment scores dataframe.
#' @param history List object containing the history of previous steps.
#' @param isPositive Logical, indicating if the marker state is positive.
#' @param depth Integer, the depth of the node in the tree.
#' @param usedmarkers Vector of markers already used in previous nodes.
#' @param path The path from the root to the current node.
#'
#' @return An object of class 'GatingTreeNode'.
#'
#' @examples
#' \dontrun{
#' marker = "CD4"
#' state = c(1, 0, 2)
#' indices = 1:100
#' markers = c("CD4", "CD8")
#' entropy = 0.5
#' enrichment = 0.7
#' scores = data.frame(score=runif(3))
#' history = list()
#' depth = 1
#' usedmarkers = c("CD4")
#' path = "rootNode"
#' createChildNode(marker, state, indices, markers, entropy, enrichment,
#' scores, scores, history, TRUE, depth, usedmarkers, path)
#' }
#'
#' @keywords internal
#' @family GatingTree
#' @noRd

createChildNode <- function(marker, current_marker_state, indices, available_markers, current_entropy, current_enrichment, average_proportion,entropy_scores, enrichment_scores, history, isPositive = TRUE, depth, usedmarkers, path) {
    # Determine the suffix based on marker state
    suffix <- if(isPositive) "pos" else "neg"
    markerName <- paste(marker, suffix, sep='.')

    structure(list(
        CurrentNodeName = markerName,
        CurrentNodeGate = markerName,
        CurrentPath = c(path, markerName),
        CurrentMarkerState = current_marker_state,
        CurrentNodeIndices = indices,
        NodeEnrichmentDf = enrichment_scores,
        NodeEntropyDf = entropy_scores,
        CurrentEntropy = current_entropy,
        CurrentEnrichment = current_enrichment,
        AverageProportion = average_proportion,
        History = history,
        Children = list(),
        AvailableMarkers = setdiff(available_markers, marker),
        Depth = depth,
        Terminated = FALSE,
        Leaf = FALSE,
        UsedMarkers = usedmarkers
    ), class = "GatingTreeNode")
}

#' Add a Child Node to a Gating Tree
#'
#' This function adds a new child node to a specified location in a gating tree.
#'
#' @param rootNode The root node or any node acting as a root in a sub-tree.
#' @param childNode The child node to be added.
#' @param path The path where the child node should be added.
#'
#' @return The modified node with the new child added.
#'
#' @examples
#' \dontrun{
#' rootNode <- createGatingTreeObject(...) # Setup initial node
#' childNode <- createChildNode(...) # Create a child node
#' path <- "rootNode"
#' addChildNode(rootNode, childNode, path)
#' }
#'
#' @keywords internal
#' @family GatingTree
#' @noRd

addChildNode <- function(rootNode, childNode, path) {
    if (length(path) == 1) {
        # Base case: add the child node directly here
        rootNode$Children[[childNode$CurrentNodeName]] <- childNode
        return(rootNode)
    } else {
        # Recursive case: navigate down the path
        nextNodeName <- path[1]
        if (!is.null(rootNode$Children[[nextNodeName]])) {
            rootNode$Children[[nextNodeName]] <- addChildNode(rootNode$Children[[nextNodeName]], childNode, path[-1])
        } else {
            stop("Path does not exist in the tree")
        }
        return(rootNode)
    }
}

#' Retrieve a Node from a Gating Tree
#'
#' This function traverses a gating tree structure to retrieve a node at a specified path.
#'
#' @param gatingTreeObject The root object of the gating tree containing child nodes.
#' @param path A vector of node names specifying the path to the desired node.
#'
#' @return The node object at the specified path.
#'
#' @examples
#' \dontrun{
#' gatingTreeObject <- createGatingTreeObject(...) # assuming this function exists
#' path <- c("rootNode", "CD4pos")
#' getNode(gatingTreeObject, path)
#' }
#'
#' @export
#' @family GatingTree

getNode <- function(gatingTreeObject, path) {
    currentNode <- gatingTreeObject
    for (i in 1:length(path)) {
        nodeName <- path[i]
        currentNode <- currentNode$Children[[nodeName]]
        if (is.null(currentNode)) {
            stop(paste("No node found at path:", paste(nodeName, collapse = " -> ")))
        }
    }
    return(currentNode)
}


#' Convert Gating Tree Node to Data Tree Node
#'
#' This function converts a node from a gating tree into a data tree node format,
#' suitable for further processing or visualization.
#'
#' @param node A node from a gating tree, usually containing statistical information
#' and child nodes.
#' @param pathName A string representing the node path, defaults to "rootNode".
#'
#' @return A data tree node object structured for hierarchical data representations.
#'
#' @examples
#' \dontrun{
#' # Assuming `node` is an object from a gating analysis tree:
#' convertedNode <- convertToDataTree(node)
#' }
#'
#' @export
#' @importFrom data.tree Node
#' @family GatingTree Visualization

convertToDataTree <- function(node, pathName = "rootNode") {
    dataTreeNode <- Node$new(pathName)
    dataTreeNode$CurrentNodeName <- node$CurrentNodeName
    dataTreeNode$CurrentNodeGate <- node$CurrentNodeGate
    dataTreeNode$CurrentNodeStats <- node$CurrentNodeStats
    dataTreeNode$CurrentEnrichment <- node$CurrentEnrichment
    dataTreeNode$CurrentEntropy <- node$CurrentEntropy
    dataTreeNode$AverageProportion <- node$AverageProportion


    for(childName in names(node$Children)) {
        childNode <- node$Children[[childName]]
        childDataTreeNode <- convertToDataTree(childNode, childName)
        dataTreeNode$AddChildNode(childDataTreeNode)
    }
    
    return(dataTreeNode)
}


#' Collect Markers from a Gating Tree
#'
#' This function collects all markers used in a gating tree.
#'
#' @param tree The gating tree structure.
#'
#' @return A vector of unique markers used in the tree.
#' @keywords internal
#' @examples
#' \dontrun{
#' collect_markers(tree)
#'}
#' @family GatingTree
#' @noRd
collect_markers <- function(tree) {
    markers <- vector("character")
    if (!is.null(tree$Marker)) {
        markers <- c(markers, tree$Marker)
    }
    if (!is.null(tree$Left)) {
        markers <- c(markers, collect_markers(tree$Left))
    }
    if (!is.null(tree$Right)) {
        markers <- c(markers, collect_markers(tree$Right))
    }
    return(markers)
}


#' Collect History from Gating Tree
#'
#' This function recursively collects history from all nodes in a gating tree,
#' aggregating any historical data stored at each node.
#'
#' @param tree A tree object representing the root or any subtree in a gating tree.
#'
#' @return A list containing all historical data from the tree.
#'
#' @examples
#' \dontrun{
#' # Assuming `tree` is a root node of a gating tree with history data:
#' historyData <- collect_history(tree)
#' }
#'
#' @keywords internal
#' @family GatingTree
#' @noRd
collect_history <- function(tree) {
    history <- list()
    
    # Check if the current node has a History object and add it to the history list
    if (!is.null(tree$History)) {
        history <- c(history, tree$History)
    }
    
    # Recursively collect history from each child if any exist
    if (!is.null(tree$Children)) {
        for (child_name in names(tree$Children)) {
            # Recursive call must be to the child node, not the current node's History
            child_histories <- collect_history(tree$Children[[child_name]])
            history <- c(history, child_histories)
        }
    }

    return(history)
}

#' Collect Enrichment Data from Leaf Nodes
#'
#' This function traverses a gating tree to collect enrichment data specifically
#' from leaf nodes, where no further subdivisions occur.
#'
#' @param tree A tree object that may be a root node or any subtree in a gating tree.
#'
#' @return A list of enrichment data collected from the leaf nodes of the tree.
#'
#' @examples
#' \dontrun{
#' # Assuming `tree` is a root node of a gating tree with enrichment data in leaf nodes:
#' leafEnrichments <- collect_leaf_enrichment(tree)
#' }
#'
#' @keywords internal
#' @family GatingTree
#' @noRd

collect_leaf_enrichment <- function(tree) {
    enrichment_list <- list()
    
    # Check if there are any children; if not, this is a leaf node
    if (length(tree$Children)==0) {
        # Since this is a leaf, check and collect the enrichment data if it exists
        if (!is.null(tree$History) && !is.null(tree$History$enrichment)) {
            enrichment_list <- c(enrichment_list, list(tree$History$enrichment))
        }
    } else {
        # Recursively collect enrichment from each child
        for (child_name in names(tree$Children)) {
            child_enrichments <- collect_leaf_enrichment(tree$Children[[child_name]])
            enrichment_list <- c(enrichment_list, child_enrichments)
        }
    }

    return(enrichment_list)
}


#' Count Nodes
#'
#' This function counts the total number of nodes in a gating tree.
#'
#' @param tree The gating tree structure.
#'
#' @return The total number of nodes in the tree.
#' @keywords internal
#' @examples
#' \dontrun{
#' count_nodes(tree)
#'}
#' @family GatingTree
#' @noRd
count_nodes <- function(tree) {
    count <- 1
    if (!is.null(tree$Left)) {
        count <- count + count_nodes(tree$Left)
    }
    if (!is.null(tree$Right)) {
        count <- count + count_nodes(tree$Right)
    }
    return(count)
}

