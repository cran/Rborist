# Copyright (C)  2012-2016   Mark Seligman
##
## This file is part of ArboristBridgeR.
##
## ArboristBridgeR is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## ArboristBridgeR is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

PreFormat.default <- function(x) {
  # Argument checking:
  if (any(is.na(x)))
    stop("NA not supported in design matrix")

  predBlock <- PredBlock(x)
  rowRank <- .Call("RcppRowRank", predBlock)

  preTrain <- list(
    predBlock = predBlock,
    rowRank = rowRank
  )
  class(preTrain) <- "PreFormat"

  preTrain
}