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

#' Perform Random Forest Analysis Using GatingTree Node Data
#'
#' This function trains a Random Forest model using node-level percentage data
#' from a FlowObject that includes a GatingTree. The fitted Random Forest model
#' and associated metadata are then stored in the `@Model` slot of the same
#' FlowObject.
#'
#' @param train_x A FlowObject containing a GatingTree and node-level
#'   percentage data in `train_x@Gating$PrunedGatingTreeNodePercentages`.
#' @param node_paths A character vector of node paths to include as predictors.
#'   If `NULL`, all columns containing "logdata" are used.
#' @param ntree Number of trees to grow in the Random Forest. Defaults to 100.
#' @param mtry Number of variables randomly sampled at each split in the Random Forest. The defaul is NULL and uses the automatic selection implemented in randomForest.
#' @param filter A character to spercify how to filter nodes. The option 'all' will consider all enrichment score, entropy, and average proportion of node with equal weight. The option 'enrichment' will use enrichment score only.
#' @param q Number between 0 and 1. Quantile value for the filter parameter to select top nodes. For example, 0.95 will select nodes above the 95th percentile.
#' @param class_weight Logical. Whether to consider class weight. The default is FALSE.
#' @param expr_group Optional. A character vector to specify the experimental group when class weight is considered.
#' @param ctrl_group Optional. A character vector to specify the control group when class weight is considered.
#' @return The updated FlowObject with a new Random Forest model stored in the
#'   `@Model` slot under the name `"RandomForestGatingTree"`.
#' @export
#' @examples
#' \dontrun{
#'   train_x <- GatingTreeRandomForest(train_x)
#' }
#' @importFrom randomForest randomForest importance
#' @family GatingTree Random Forest Analysis

GatingTreeRandomForest <- function(train_x, node_paths = NULL, ntree = 100, mtry = NULL, expr_group = NULL, ctrl_group = NULL, class_weight = FALSE, filter = "all", q = 0.75) {
    
    if(is.null(train_x@Gating$PrunedGatingTreeNodePercentages)){
        stop("Perform PruneGatingTree. \n")
    }
    
    train_df <- train_x@Gating$PrunedGatingTreeNodePercentages
    
    if(!is.na(filter)){
        if(filter=='all'){
            prndf <- train_x@Gating$PrunedGatingTreeDF
            enrichment <- prndf$max_enrichment
            entropy <- prndf$entropy
            average_proportion <- prndf$average_proportion
            
            enrichment_std <- as.vector(scale(enrichment))
            entropy_std <- as.vector(scale((10^entropy)^(-1)))
            average_proportion_std <- as.vector(scale(average_proportion))
            composite_score <- enrichment_std + entropy_std +average_proportion_std
             prndf$composite_score <- composite_score
            top_combination <- prndf$markers_up_to_max[prndf$composite_score > quantile(prndf$composite_score, q)]
            #train_df <- train_df[,top_combination]
        }
        if(filter=='enrichment'){
            prndf <- train_x@Gating$PrunedGatingTreeDF
            enrichment <- prndf$max_enrichment
            top_combination <- prndf$markers_up_to_max[enrichment > quantile(enrichment, q)]

        }
    }else{
        prndf <- train_x@Gating$PrunedGatingTreeDF
        top_combination <- prndf$markers_up_to_max
    }

    
    node_paths <- top_combination
    train_sampledef <- train_x@sampledef$sampledef
    train_sampledef <- train_sampledef[,c('file', 'group')]
    
    train_exprs <- train_x@Data
    neg_df <- train_x@QCdata$negative_gate_def
    
    exprs_matrix <- as.matrix(train_exprs[, neg_df$variable])
    marker_cols <- colnames(exprs_matrix)#[grepl(pattern = 'logdata', colnames(exprs_matrix))]
    thresholds_vec <- sapply(marker_cols, function(m) neg_df[neg_df$variable == m, "negative.gate"])
    
    make_combinations_matrix <- function(node_paths, marker_cols) {
        mat <- matrix(0, nrow = length(node_paths), ncol = length(marker_cols))
        for (i in seq_along(node_paths)) {
            markers <- strsplit(node_paths[i], "_")[[1]]
            for (m in markers) {
                parts <- strsplit(m, "\\.")[[1]]
                marker_base <- paste(parts[1], "logdata", sep=".")
                condition <- parts[3]
                marker_idx <- match(marker_base, marker_cols)  # find the column index
                if (!is.na(marker_idx)) {
                    if (condition == "pos") {
                        mat[i, marker_idx] <- 2
                    } else {  # "neg"
                        mat[i, marker_idx] <- 1
                    }
                }
            }
        }
        mat
    }
    
    combinations_matrix <- make_combinations_matrix(node_paths, marker_cols)
    
    results_matrix <- apply_gating_conditions(
        data = exprs_matrix,
        combinations_matrix = combinations_matrix,
        thresholds = thresholds_vec
    )
    results_matrix <- as.data.frame(results_matrix)
    results_df <- data.frame(file = train_exprs$file, results_matrix)
    
    results_summary <- results_df %>%
    group_by(file) %>%
    summarise(across(where(is.logical), ~ mean(.x), .names = "{.col}"))
    results_summary[is.na(results_summary)] <- 0

    logic <- colnames(results_summary) %in% 'file'
    train_matrix <- results_summary[, !logic]
    train_matrix[is.na(train_matrix)] <- 0
    
    original_names <- node_paths
    simplified_names <- paste("V", seq_along(original_names), sep = "")
    name_mapping <- setNames(simplified_names, original_names)

    colnames(train_matrix) <- name_mapping
    train_matrix$file <- as.character(results_summary$file)
    train_sampledef$file <- as.character(train_sampledef$file)

    train_matrix_df <- merge(train_matrix, train_sampledef, by = "file", all.x = TRUE)

    train_matrix_df$group <- as.factor(train_matrix_df$group)
    train_matrix_df_file <- train_matrix_df$file
    train_matrix_df <- train_matrix_df[, colnames(train_matrix_df)!='file']
    group_name <- "group"
    features <- simplified_names
    rf_formula <- as.formula(paste(group_name, "~", paste(features, collapse = " + ")))

    if(class_weight){
        total_samples <- nrow(train_matrix_df)
        weight_expr <- 1 / (sum(train_matrix_df$group == expr_group) / total_samples)
        weight_ctrl <- 1 / (sum(train_matrix_df$group == ctrl_group) / total_samples)
        
        classwt_vec <- setNames(c(weight_expr, weight_ctrl), c(expr_group, ctrl_group))

        if(!is.null(mtry)){
            rf_model <- randomForest(group ~ ., data = train_matrix_df, ntree = ntree,
                                     classwt = classwt_vec, mtry = mtry)
        }else{
            rf_model <- randomForest(group ~ ., data = train_matrix_df, ntree = ntree,
                                     classwt = classwt_vec)
            
        }

    }else{
        
        if(!is.null(mtry)){
            rf_model <- randomForest(group ~ ., data = train_matrix_df, ntree = ntree, importance = TRUE, mtry = mtry)
        }else{
            rf_model <- randomForest(group ~ ., data = train_matrix_df, ntree = ntree, importance = TRUE)
        }
        
    }
    
    reversed_name_mapping <- setNames(names(name_mapping), name_mapping)
    colnames(train_matrix_df[,name_mapping]) <- reversed_name_mapping[colnames(train_matrix_df[,name_mapping])]
    train_matrix_df$file <- train_matrix_df_file

    current_names <- colnames(train_matrix_df)
    idx <- which(current_names %in% simplified_names)
    current_names[idx] <- original_names[match(current_names[idx], simplified_names)]
    colnames(train_matrix_df) <- current_names
    train_matrix_df$file <- train_matrix_df_file
    importances_df <- importance(rf_model)
    importances <- importances_df[,'MeanDecreaseGini']
    names(importances) <- rownames(importances_df)
    
    names(importances) <- original_names[match(names(importances), simplified_names)]
    importances <- importances[order(importances, decreasing = TRUE)]

    out <- list( formula = rf_formula, importance_score = importances, name_mapping = name_mapping, train_data = train_matrix_df, node_paths = node_paths)

    
    train_x@Model[['GatingTreeRandomForest']] <- rf_model
    train_x@Model[['GatingTreeRandomForestData']] <- out
    
    return(train_x)
    
}
    
#' Predict Using GatingTreeRandomForest Model
#'
#' This function uses the GatingTreeRandomForest model stored in a `FlowObject` to make predictions on a new dataset.
#' It prepares the test dataset, applies the gating conditions specified in the model, and then uses the Random Forest model to predict outcomes.
#'
#' @param train_x A `FlowObject` that contains a previously fitted `GatingTreeRandomForest` model.
#' @param test_y A `FlowObject` containing the test dataset for prediction.
#'
#' @return A list containing:
#' \itemize{
#'   \item predicted_scores: The probabilities predicted by Random Forest.
#'   \item test_data: The test data used for Random Forest.
#'   \item results_pred: A data frame for predicted scores.
#' }

#'
#' @examples
#' \dontrun{
#'   train_x <- GatingTreeRandomForest(train_x)
#'   predictions <- predictGatingTreeRandomForest(train_x, test_y)
#'   head(predictions$predicted_scores)
#' }
#'
#' @export
#' @importFrom dplyr group_by summarise across where %>%
#' @family GatingTree Random Forest Analysis

predictGatingTreeRandomForest <- function(train_x, test_y) {
    
    if(is.null(train_x@Model[['GatingTreeRandomForest']])){
        stop("Perform GatingTreeRandomForest. \n")
    }
    
    rf_model <- train_x@Model[['GatingTreeRandomForest']]
    trained_model_data <- train_x@Model[['GatingTreeRandomForestData']]
    name_mapping <- trained_model_data[['name_mapping']]
    node_paths <- trained_model_data[['node_paths']]
    reversed_name_mapping <- setNames(names(name_mapping), name_mapping)
    
    test_exprs <- test_y@Data
    neg_df <- test_y@QCdata$negative_gate_def
    test_sampledef <- test_y@sampledef$sampledef
    
    exprs_matrix <- as.matrix(test_exprs[, neg_df$variable])
    marker_cols <- colnames(exprs_matrix)#[grepl(pattern = 'logdata', colnames(test_exprs))]
    thresholds_vec <- sapply(marker_cols, function(m) neg_df[neg_df$variable == m, "negative.gate"])
    
    make_combinations_matrix <- function(node_paths, marker_cols) {
        mat <- matrix(0, nrow = length(node_paths), ncol = length(marker_cols))
        for (i in seq_along(node_paths)) {
            markers <- strsplit(node_paths[i], "_")[[1]]
            for (m in markers) {
                parts <- strsplit(m, "\\.")[[1]]
                marker_base <- paste(parts[1], "logdata", sep=".")
                condition <- parts[3]
                marker_idx <- match(marker_base, marker_cols)  # find the column index
                if (!is.na(marker_idx)) {
                    if (condition == "pos") {
                        mat[i, marker_idx] <- 2
                    } else {  # "neg"
                        mat[i, marker_idx] <- 1
                    }
                }
            }
        }
        mat
    }
    
    combinations_matrix <- make_combinations_matrix(node_paths, marker_cols)
    
    results_matrix <- apply_gating_conditions(
        data = exprs_matrix,
        combinations_matrix = combinations_matrix,
        thresholds = thresholds_vec
    )
    results_matrix <- as.data.frame(results_matrix)
    colnames(results_matrix) <- name_mapping[node_paths]
    results_df <- data.frame(file = test_exprs$file, results_matrix)
    
    results_summary <- results_df %>%
    group_by(file) %>%
    summarise(across(where(is.logical), ~ mean(.x), .names = "{.col}"))
    results_summary[is.na(results_summary)] <- 0

    test_features <- setdiff(colnames(results_summary), "file")
    test_data <- results_summary[, test_features]
    predicted_scores <- predict(rf_model, newdata = test_data, type = "prob")
    predicted_scores <- as.data.frame(predicted_scores)
    predicted_scores$file <- results_summary$file
    test_data$file <- results_summary$file
    
    lg <- colnames(results_summary)=='file'

    colnames(results_summary)[!lg] <- reversed_name_mapping[colnames(results_summary)[!lg]]
    
    results_summary$predicted_scores <- predicted_scores[, 1]
    results_pred <- results_summary[, c('file','predicted_scores')]
    results_pred <- merge(test_sampledef, results_pred, by = 'file')
    
    out <- list(predicted_scores = predicted_scores,
    test_data = test_data, results_pred = results_pred, positive_class = colnames(predicted_scores)[1])
    
}
