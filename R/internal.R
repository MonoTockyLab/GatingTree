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

#' Logarithmically Transforms a Vector
#'
#' Applies a logarithmic transformation to each element of a numeric vector.
#' Elements greater than 1 are transformed using log10, while elements less than
#' or equal to 1 are set to 0.
#'
#' @param x A numeric vector to be transformed.
#'
#' @return A numeric vector with the transformed values.
#' @export
#' @examples
#' logSingleData(c(10, 1, 0.5, 20))

logSingleData <- function(x) {
  sapply(x, function(value) {
    if (value > 1) log10(value) else 0
  })
}


#' Logarithmically Transforms Each Column of a Matrix
#'
#' Transforms each column of a numeric matrix by applying a logarithmic
#' transformation adjusted by a quantile and a safety margin.
#'
#' @param data A numeric matrix whose columns are to be transformed.
#' @param quant_val A quantile value used for the transformation adjustment.
#' @param safety_margin A numeric value added for safety to ensure positive logarithm domain.
#'
#' @return A numeric matrix with each column transformed.
#' @export
#' @examples
#' data <- matrix(c(1, 2, 3, 4, 5, 6), nrow=2)
#' logTransformData(data)
logTransformData <- function(data, quant_val = 0.001, safety_margin = 100) {
    logSingleData <- function(x) {
      x_log <- numeric(length(x))
      logic <- x > 1
      x_log[logic] <- log10(x[logic])
      x_log[!logic] <- 0
      return(x_log)
    }

    
  cols <- ncol(data)
  
  for (i in 1:cols) {
    current_column <- data[, i]
    s <- quantile(current_column, probs = quant_val, type = 1)
    data[, i] <- logSingleData(current_column - s + safety_margin)
  }
  data
}



#' Moderately Log Transforms a Vector with Adjustments
#'
#' Applies a moderated logarithmic transformation to a numeric vector.
#' Adjustments are made based on a specified quantile and a floor value.
#' All transformed values are guaranteed to be in the domain suitable for log10.
#'
#' @param x A numeric vector to be transformed.
#' @param q The quantile adjustment value.
#' @param f The minimum value of the adjustment floor, used to avoid negative or zero logarithm domain.
#'
#' @return A numeric vector with the logarithmically transformed values.
#' @export
#' @examples
#' moderate_log_transform(c(1, 2, 3, 4, 5), 0.5, 0.1)
moderate_log_transform <- function(x, q, f) {
  e <- max(f, 0.1 * abs(q))
  adjusted_val <- x - q + e
  
  adjusted_val[adjusted_val <= 1] <- 1
  
  x_log <- log10(adjusted_val)
  
  x_log
}


