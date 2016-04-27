// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppMonolith.cc

   @brief Auto-generated monolith of 'rcpp' files.

   @author Mark Seligman
 */

// This abomination was created by concatening several 'Rcpp' files
// into a single module, in an effort to reduce the size of the
// generated object.
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppExport.cc

   @brief C++ interface to R entry for export methods.

   @author Mark Seligman
 */

#include <Rcpp.h>
using namespace Rcpp;

using namespace std;
//#include <iostream>

#include "rcppPredblock.h"
#include "rcppForest.h"
#include "rcppLeaf.h"
#include "forest.h"
#include "bv.h"
#include "leaf.h"


/**
   @brief Recasts 'pred' field of nonterminals to front-end facing values.

   @return void.
 */
void PredTree(const int predMap[], std::vector<unsigned int> &pred, std::vector<unsigned int> &bump) {
  for (unsigned int i = 0; i < pred.size(); i++) {
    if (bump[i] > 0) { // terminal 'pred' values do not reference predictors.
      unsigned int predCore = pred[i];
      pred[i] = predMap[predCore];
    }
  }
}


/**
   @brief Prepares predictor field for export by remapping to front-end indices.
 */
void PredExport(const int predMap[], std::vector<std::vector<unsigned int> > &predTree, std::vector<std::vector<unsigned int> > &bumpTree) {
  for (unsigned int tIdx = 0; tIdx < predTree.size(); tIdx++) {
    PredTree(predMap, predTree[tIdx], bumpTree[tIdx]);
  }
}


/**
   @brief Exports core data structures as vector of per-tree vectors.

   @return List with common and regression-specific members.
 */
RcppExport SEXP ExportReg(SEXP sForest, SEXP sLeaf, IntegerVector predMap) {

  // Instantiates the forest-wide data structures as long vectors, then
  // distributes per tree.
  //
  std::vector<unsigned int> nodeOrigin, facOrigin, splitBV;
  std::vector<ForestNode> *forestNode;
  ForestUnwrap(sForest, nodeOrigin, facOrigin, splitBV, forestNode);

  unsigned int nTree = nodeOrigin.size();
  std::vector<std::vector<unsigned int> > predTree(nTree), bumpTree(nTree);
  std::vector<std::vector<double > > splitTree(nTree);
  ForestNode::Export(nodeOrigin, *forestNode, predTree, bumpTree, splitTree);
  PredExport(predMap.begin(), predTree, bumpTree);
  
  std::vector<std::vector<unsigned int> > facSplitTree(nTree);
  BVJagged::Export(facOrigin, splitBV, facSplitTree);

  std::vector<double> yRanked;
  std::vector<unsigned int> leafOrigin;
  std::vector<LeafNode> *leafNode;
  std::vector<BagRow> *bagRow;
  unsigned int rowTrain;
  std::vector<unsigned int> rank;
  LeafUnwrapReg(sLeaf, yRanked, leafOrigin, leafNode, bagRow, rowTrain, rank);

  std::vector<std::vector<unsigned int> > rowTree(nTree), sCountTree(nTree);
  std::vector<std::vector<double> > scoreTree(nTree);
  std::vector<std::vector<unsigned int> > extentTree(nTree);
  std::vector<std::vector<unsigned int> > rankTree(nTree);
  LeafReg::Export(leafOrigin, *leafNode, *bagRow, rank, rowTree, sCountTree, scoreTree, extentTree, rankTree);

  List outBundle = List::create(
				_["rowTrain"] = rowTrain,
				_["pred"] = predTree,
				_["bump"] = bumpTree,
				_["split"] = splitTree,
				_["facSplit"] = facSplitTree,
				_["row"] = rowTree,
				_["sCount"] = sCountTree,
				_["score"] = scoreTree,
				_["extent"] = extentTree,
				_["rank"] = rankTree
				);
  outBundle.attr("class") = "ExportReg";

  return outBundle;
}


/**
   @brief Exports core data structures as vector of per-tree vectors.

   @return List with common and classification-specific members.
 */
RcppExport SEXP ExportCtg(SEXP sForest, SEXP sLeaf, IntegerVector predMap) {
  std::vector<unsigned int> nodeOrigin, facOrigin, splitBV;
  std::vector<ForestNode> *forestNode;
  ForestUnwrap(sForest, nodeOrigin, facOrigin, splitBV, forestNode);

  unsigned int nTree = nodeOrigin.size();
  std::vector<std::vector<unsigned int> > predTree(nTree), bumpTree(nTree);
  std::vector<std::vector<double > > splitTree(nTree);
  ForestNode::Export(nodeOrigin, *forestNode, predTree, bumpTree, splitTree);
  PredExport(predMap.begin(), predTree, bumpTree);
  
  std::vector<std::vector<unsigned int> > facSplitTree(nTree);
  BVJagged::Export(facOrigin, splitBV, facSplitTree);

  std::vector<unsigned int> leafOrigin;
  std::vector<LeafNode> *leafNode;
  std::vector<BagRow> *bagRow;
  unsigned int rowTrain;
  std::vector<double> weight;
  CharacterVector yLevel;
  LeafUnwrapCtg(sLeaf, leafOrigin, leafNode, bagRow, rowTrain, weight, yLevel);

  std::vector<std::vector<unsigned int> > rowTree(nTree), sCountTree(nTree);
  std::vector<std::vector<double> > scoreTree(nTree);
  std::vector<std::vector<unsigned int> > extentTree(nTree);
  std::vector<std::vector<double> > weightTree(nTree);
  LeafCtg::Export(leafOrigin, *leafNode, *bagRow, weight, yLevel.length(), rowTree, sCountTree, scoreTree, extentTree, weightTree);

  List outBundle = List::create(
				_["rowTrain"] = rowTrain,
				_["pred"] = predTree,
				_["bump"] = bumpTree,
				_["split"] = splitTree,
				_["facSplit"] = facSplitTree,
				_["row"] = rowTree,
				_["sCount"] = sCountTree,
				_["score"] = scoreTree,
				_["extent"] = extentTree,
				_["yLevel"] = yLevel,
				_["weight"] = weightTree
				);
  outBundle.attr("class") = "ExportCtg";

  return outBundle;
}


unsigned int NTree(SEXP sExp) {
  List exp(sExp);
  if (!exp.inherits("ExportCtg") && !exp.inherits("ExportReg"))
    stop("Unrecognized export object");

  std::vector<std::vector<unsigned int> > pred = exp["pred"];

  return pred.size();
}


/**
   @brief Only the scores are of interest to ForestFloor.

   @param forestCore is the exported core image.

   @param tIdx is the tree index.

   @return Vector of score values.
 */
RcppExport SEXP FFloorLeafReg(SEXP sForestCore, unsigned int tIdx) {
  List forestCore(sForestCore);
  std::vector<std::vector<double> > score = forestCore["score"];
  List ffLeaf = List::create(
     _["score"] = score[tIdx]
    );
  
  ffLeaf.attr("class") = "FFloorLeafReg";
  return ffLeaf;
}


/**
   @brief Only the scores and weights are of interest to ForestFloor.

   @param forestCore is the exported core image.

   @param tIdx is the tree index.

   @return Vector of score values.
 */
RcppExport SEXP FFloorLeafCtg(SEXP sForestCore, unsigned int tIdx) {
  List forestCore(sForestCore);
  std::vector<std::vector<double> > score = forestCore["score"];
  std::vector<std::vector<double> > weight = forestCore["weight"];
  unsigned int leafCount = score[tIdx].size();
  NumericMatrix weightOut = NumericMatrix(weight[tIdx].size() / leafCount, leafCount, weight[tIdx].begin());
  List ffLeaf = List::create(
     _["score"] = score[tIdx],
     _["weight"] = transpose(weightOut)
     );

  ffLeaf.attr("class") = "FFloorLeafCtg";
  return ffLeaf;
}


/**
 */
RcppExport SEXP FFloorInternal(SEXP sForestCore, unsigned int tIdx) {
  List forestCore(sForestCore);
  std::vector<std::vector<unsigned int> > predTree = forestCore["pred"];
  std::vector<std::vector<unsigned int> > bumpTree = forestCore["bump"];
  std::vector<std::vector<double > > splitTree = forestCore["split"];
  std::vector<std::vector<unsigned int> > facSplitTree = forestCore["facSplit"];
  IntegerVector incrL(bumpTree[tIdx].begin(), bumpTree[tIdx].end());
  IntegerVector predIdx(predTree[tIdx].begin(), predTree[tIdx].end());
  List ffTree = List::create(
     _["pred"] = ifelse(incrL == 0, -(predIdx + 1), predIdx),
     _["daughterL"] = incrL,
     _["daughterR"] = ifelse(incrL == 0, 0, incrL + 1),
     _["split"] = splitTree[tIdx],
     _["facSplit"] = facSplitTree[tIdx]
     );

  ffTree.attr("class") = "FFloorTree";
  return ffTree;
}


RcppExport SEXP FFloorBag(SEXP sForestCore, int tIdx) {
  List forestCore(sForestCore);
  unsigned int rowTrain = as<unsigned int>(forestCore["rowTrain"]);
  std::vector<std::vector<unsigned int> > rowTree = forestCore["row"];
  std::vector<std::vector<unsigned int> > sCountTree = forestCore["sCount"];
  IntegerVector bag = IntegerVector(rowTrain);
  IntegerVector row(rowTree[tIdx].begin(), rowTree[tIdx].end());
  IntegerVector sCount(sCountTree[tIdx].begin(), sCountTree[tIdx].end());
  bag[row] = sCount;

  return bag;
}


/**
 */
RcppExport SEXP FFloorTreeReg(SEXP sCoreReg, unsigned int tIdx) {
  List ffReg = List::create(
    _["internal"] = FFloorInternal(sCoreReg, tIdx),
    _["leaf"] = FFloorLeafReg(sCoreReg, tIdx),
    _["bag"] = FFloorBag(sCoreReg, tIdx)
  );

  ffReg.attr("class") = "FFloorTreeReg";
  return ffReg;
}


/**
 */
RcppExport SEXP FFloorTreeCtg(SEXP sCoreCtg, unsigned int tIdx) {
  List ffCtg = List::create(
    _["internal"] = FFloorInternal(sCoreCtg, tIdx),
    _["leaf"] = FFloorLeafCtg(sCoreCtg, tIdx),
    _["bag"] = FFloorBag(sCoreCtg, tIdx)
  );

  ffCtg.attr("class") = "FFloorTreeCtg";
  return ffCtg;
}


/**
 */
RcppExport SEXP FFloorReg(SEXP sForest, SEXP sLeaf, IntegerVector predMap, List predLevel) {
  SEXP sCoreReg = ExportReg(sForest, sLeaf, predMap);
  unsigned int nTree = NTree(sCoreReg);
  List trees(nTree);
  for (unsigned int tIdx = 0; tIdx < nTree; tIdx++) {
     trees[tIdx] = FFloorTreeReg(sCoreReg, tIdx);
  }

  int facCount = predLevel.length();
  IntegerVector facMap(predMap.end() - facCount, predMap.end());
  List ffe = List::create(
    _["facMap"] = facMap,
    _["predLevel"] = predLevel,
    _["tree"] = trees
  );
  ffe.attr("class") = "ForestFloorReg";
  return ffe;
}


/**
 */
RcppExport SEXP FFloorCtg(SEXP sForest, SEXP sLeaf, IntegerVector predMap, List predLevel) {
  SEXP sCoreCtg = ExportCtg(sForest, sLeaf, predMap);
  unsigned int nTree = NTree(sCoreCtg);
  List trees(nTree);
  for (unsigned int tIdx = 0; tIdx < nTree; tIdx++) {
    trees[tIdx] = FFloorTreeCtg(sCoreCtg, tIdx);
  }

  int facCount = predLevel.length();
  IntegerVector facMap(predMap.end() - facCount, predMap.end());
  List coreCtg(sCoreCtg);
  List ffe = List::create(
   _["facMap"] = facMap,
   _["predLevel"] = predLevel,
   _["yLevel"] = as<CharacterVector>(coreCtg["yLevel"]),
   _["tree"] = trees
  );
  ffe.attr("class") = "ForestFloorCtg";
  return ffe;
}


/**
   @brief Structures forest summary for analysis by ForestFloor package.

   @param sForest is the Forest summary.

   @return ForestFloorExport as List.
 */
RcppExport SEXP RcppForestFloorExport(SEXP sArbOut) {
  List arbOut(sArbOut);
  if (!arbOut.inherits("Rborist")) {
    warning("Expecting an Rborist object");
    return List::create(0);
  }

  IntegerVector predMap;
  List predLevel;
  SignatureUnwrap(arbOut["signature"], predMap, predLevel);

  List leaf((SEXP) arbOut["leaf"]);
  if (leaf.inherits("LeafReg"))  {
    return FFloorReg(arbOut["forest"], arbOut["leaf"], predMap, predLevel);
  }
  else if (leaf.inherits("LeafCtg")) {
    return FFloorCtg(arbOut["forest"], arbOut["leaf"], predMap, predLevel);
  }
  else {
    warning("Unrecognized forest type.");
    return List::create(0);
  }
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppForest.cc

   @brief C++ interface to R entry for Forest methods.

   @author Mark Seligman
 */


#include <RcppCommon.h>

#include "forest.h"

namespace Rcpp {
  template <> SEXP wrap(const std::vector<ForestNode> &);
  template <> std::vector<ForestNode>* as(SEXP);
}

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#include "rcppForest.h"


template<> SEXP Rcpp::wrap(const std::vector<ForestNode> &forestNode) {
  XPtr<const std::vector<ForestNode> > extWrap(new std::vector<ForestNode>(forestNode), true);
 return wrap(extWrap);
}


template <> std::vector<ForestNode>* Rcpp::as(SEXP sFN) {
  Rcpp::XPtr<std::vector<ForestNode> > xp(sFN);
  return (std::vector<ForestNode>*) xp;
}


//#include <iostream>

RcppExport SEXP ForestWrap(const std::vector<unsigned int> &origin, const std::vector<unsigned int> &facOrigin, const std::vector<unsigned int> &facSplit, const std::vector<ForestNode> &forestNode) {
  List forest = List::create(
     _["forestNode"] = forestNode,
     _["origin"] = origin,
     _["facOrig"] = facOrigin,
     _["facSplit"] = facSplit);
  forest.attr("class") = "Forest";

  return forest;
}


/**
   @brief Exposes front-end Forest fields for transmission to core.

   @return void.
 */
void ForestUnwrap(SEXP sForest, std::vector<unsigned int> &_origin, std::vector<unsigned int> &_facOrig, std::vector<unsigned int> &_facSplit, std::vector<ForestNode> *&_forestNode) {
  List forest(sForest);
  if (!forest.inherits("Forest"))
    stop("Expecting Forest");

  _origin = as<std::vector<unsigned int> >(forest["origin"]);
  _facOrig = as<std::vector<unsigned int> >(forest["facOrig"]);
  _facSplit = as<std::vector<unsigned int> >(forest["facSplit"]);
  _forestNode = as<std::vector<ForestNode> *>(forest["forestNode"]);
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppLeaf.cc

   @brief C++ interface to R entry for Leaf methods.

   @author Mark Seligman
 */


#include <RcppCommon.h>

#include "leaf.h"

namespace Rcpp {
  template <> SEXP wrap(const std::vector<LeafNode> &);
  template <> SEXP wrap(const std::vector<BagRow> &);
  template <> std::vector<LeafNode>* as(SEXP);
  template <> std::vector<BagRow>* as(SEXP);
}

#include "rcppLeaf.h"

template <> SEXP Rcpp::wrap(const std::vector<LeafNode> &leafNode) {
  XPtr<const std::vector<LeafNode> > extWrap(new std::vector<LeafNode>(leafNode), true);

  return wrap(extWrap);
}


template <> SEXP Rcpp::wrap(const std::vector<BagRow> &bagRow) {
  XPtr<const std::vector<BagRow> > extWrap(new std::vector<BagRow>(bagRow), true);

  return wrap(extWrap);
}


template <> std::vector<LeafNode>* Rcpp::as(SEXP sLNReg) {
  Rcpp::XPtr<std::vector<LeafNode> > xp(sLNReg);
  return (std::vector<LeafNode>*) xp;
}


template <> std::vector<BagRow>* Rcpp::as(SEXP sLNReg) {
  Rcpp::XPtr<std::vector<BagRow> > xp(sLNReg);
  return (std::vector<BagRow>*) xp;
}


/**
   @brief Wraps core (regression) Leaf vectors for reference by front end.
 */
RcppExport SEXP LeafWrapReg(const std::vector<unsigned int> &leafOrigin, const std::vector<LeafNode> &leafNode, const std::vector<BagRow> &bagRow, unsigned int rowTrain, const std::vector<unsigned int> &rank, const std::vector<double> &yRanked) {
  List leaf = List::create(
   _["origin"] = leafOrigin,
   _["node"] = leafNode,
   _["bagRow"] = bagRow,
   _["rowTrain"] = rowTrain,
   _["rank"] = rank,
   _["yRanked"] = yRanked
  );
  leaf.attr("class") = "LeafReg";
  
  return leaf;
}


/**
   @brief Wraps core (classification) Leaf vectors for reference by front end.
 */
RcppExport SEXP LeafWrapCtg(const std::vector<unsigned int> &leafOrigin, const std::vector<LeafNode> &leafNode, const std::vector<BagRow> &bagRow, unsigned int rowTrain, const std::vector<double> &weight, const CharacterVector &levels) {
  List leaf = List::create(
   _["origin"] = leafOrigin,	
   _["node"] = leafNode,
   _["bagRow"] = bagRow,
   _["rowTrain"] = rowTrain,
   _["weight"] = weight,
   _["levels"] = levels
   );
  leaf.attr("class") = "LeafCtg";

  return leaf;
}


/**
   @brief Exposes front-end (regression) Leaf fields for transmission to core.

   @param sLeaf is the R object containing the leaf (list) data.

   @param _yRanked outputs the sorted response.

   @param _leafInfoReg outputs the sample ranks and counts, organized by leaf.

   @return void, with output reference parameters.
 */
void LeafUnwrapReg(SEXP sLeaf, std::vector<double> &_yRanked, std::vector<unsigned int> &_leafOrigin, std::vector<LeafNode> *&_leafNode, std::vector<BagRow> *&_bagRow, unsigned int &_rowTrain, std::vector<unsigned int> &_rank) {
  List leaf(sLeaf);
  if (!leaf.inherits("LeafReg"))
    stop("Expecting LeafReg");

  _yRanked = as<std::vector<double> >(leaf["yRanked"]);
  _leafOrigin = as<std::vector<unsigned int>>(leaf["origin"]);
  _leafNode = as<std::vector<LeafNode> *>(leaf["node"]);
  _bagRow = as<std::vector<BagRow> *>(leaf["bagRow"]);
  _rowTrain = as<unsigned int>(leaf["rowTrain"]);
  _rank = as<std::vector<unsigned int> >(leaf["rank"]);
}


/**
   @brief Exposes front-end (classification) Leaf fields for transmission to core.

   @param sLeaf is the R object containing the leaf (list) data.

   @param _weight outputs the sample weights.

   @param _levels outputs the category levels; retains as front-end object.

   @return void, with output reference parameters.
 */
void LeafUnwrapCtg(SEXP sLeaf, std::vector<unsigned int> &_leafOrigin, std::vector<LeafNode> *&_leafNode, std::vector<BagRow> *&_bagRow, unsigned int &_rowTrain, std::vector<double> &_weight, CharacterVector &_levels) {
  List leaf(sLeaf);
  if (!leaf.inherits("LeafCtg")) {
    stop("Expecting LeafCtg");
  }
  _leafOrigin = as<std::vector<unsigned int> >(leaf["origin"]);
  _leafNode = as<std::vector<LeafNode> *>(leaf["node"]);
  _bagRow = as<std::vector<BagRow> *>(leaf["bagRow"]);
  _rowTrain = as<unsigned int>(leaf["rowTrain"]);
  _weight = as<std::vector<double> >(leaf["weight"]);
  _levels = as<CharacterVector>((SEXP) leaf["levels"]);
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppPredBlock.cc

   @brief C++ interface to R entries for maintaining predictor data structures.

   @author Mark Seligman
*/

#include "rowrank.h"
// Testing only:
//#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;


void FactorRemap(IntegerMatrix &xFac, List &level, List &levelTrain);

/**
  @brief Extracts contents of a data frame into numeric and (zero-based) factor blocks.  Can be quite slow for large predictor counts, as a linked list is being walked.

  @param sX is the raw data frame, with columns assumed to be either factor or numeric.

  @param sNRow is the number of rows.

  @param sNCol is the number of columns.

  @param sFacCol is the number of factor-valued columns.

  @param sNumCol is the number of numeric-valued columns.

  @param sLevel is a vector of level counts for each column.

  @return PredBlock with separate numeric and integer matrices.
*/
RcppExport SEXP RcppPredBlockFrame(SEXP sX, SEXP sNumElt, SEXP sFacElt, SEXP sLevels, SEXP sSigTrain) {
  DataFrame xf(sX);
  IntegerVector numElt = IntegerVector(sNumElt) - 1;
  IntegerVector facElt = IntegerVector(sFacElt) - 1;
  IntegerVector levels(sLevels);
  int nRow = xf.nrows();
  int nPredFac = facElt.length();
  int nPredNum = numElt.length();
  int nPred = nPredFac + nPredNum;

  IntegerVector predMap(nPred);
  IntegerVector facCard;
  IntegerMatrix xFac;
  NumericMatrix xNum;
  if (nPredNum > 0) {
    xNum = NumericMatrix(nRow, nPredNum);
  }
  else
    xNum = NumericMatrix(0, 0);
  if (nPredFac > 0) {
    facCard = IntegerVector(nPredFac); // Compressed factor vector.
    xFac = IntegerMatrix(nRow, nPredFac);
  }
  else {
    facCard = IntegerVector(0);
    xFac = IntegerMatrix(0);
  }

  int numIdx = 0;
  int facIdx = 0;
  List level(nPredFac);
  for (int feIdx = 0; feIdx < nPred; feIdx++) {
    int card = levels[feIdx];
    if (card == 0) {
      xNum(_, numIdx) = as<NumericVector>(xf[feIdx]);
      predMap[numIdx++] = feIdx;
    }
    else {
      facCard[facIdx] = card;
      level[facIdx] = as<CharacterVector>(as<IntegerVector>(xf[feIdx]).attr("levels"));
      xFac(_, facIdx) = as<IntegerVector>(xf[feIdx]) - 1;
      predMap[nPredNum + facIdx++] = feIdx;
    }
  }

  // Factor positions must match those from training and values must conform.
  //
  if (!Rf_isNull(sSigTrain) && nPredFac > 0) {
    List sigTrain(sSigTrain);
    IntegerVector predTrain(as<IntegerVector>(sigTrain["predMap"]));
    if (!is_true(all( predMap == predTrain)))
      stop("Signature mismatch");

    List levelTrain(as<List>(sigTrain["level"]));
    FactorRemap(xFac, level, levelTrain);
  }
  List signature = List::create(
        _["predMap"] = predMap,
        _["level"] = level
	);
  signature.attr("class") = "Signature";
  
  List predBlock = List::create(
      _["colNames"] = colnames(xf),
      _["rowNames"] = rownames(xf),
      _["blockNum"] = xNum,
      _["nPredNum"] = nPredNum,
      _["blockFac"] = xFac,
      _["nPredFac"] = nPredFac,
      _["nRow"] = nRow,
      _["facCard"] = facCard,
      _["signature"] = signature
      );
  predBlock.attr("class") = "PredBlock";

  return predBlock;
}


void FactorRemap(IntegerMatrix &xFac, List &levelTest, List &levelTrain) {
  for (int col = 0; col < xFac.ncol(); col++) {
    CharacterVector colTest(as<CharacterVector>(levelTest[col]));
    CharacterVector colTrain(as<CharacterVector>(levelTrain[col]));
    if (is_true(any(colTest != colTrain))) {
      IntegerVector colMatch = match(colTest, colTrain);
      IntegerVector sq = seq(0, colTest.length() - 1);
      IntegerVector idxNonMatch = sq[is_na(colMatch)];
      if (idxNonMatch.length() > 0) {
	warning("Factor levels not observed in training:  employing proxy");
	int proxy = colTrain.length() + 1;
	colMatch[idxNonMatch] = proxy;
      }

      colMatch = colMatch - 1;  // match() is one-based.
      IntegerMatrix::Column xCol = xFac(_, col);
      IntegerVector colT(xCol);
      IntegerVector colRemap = colMatch[colT];
      xFac(_, col) = colRemap;
    }
  }
}


RcppExport SEXP RcppPredBlockNum(SEXP sX) {
  NumericMatrix blockNum(as<NumericMatrix>(sX));
  int nPred = blockNum.ncol();
  List dimnames = blockNum.attr("dimnames");
  List signature = List::create(
      _["predMap"] = seq_len(nPred) - 1,
      _["level"] = List::create(0)
  );
  signature.attr("class") = "Signature";
  
  List predBlock = List::create(
	_["colNames"] = colnames(blockNum),
	_["rowNames"] = rownames(blockNum),
	_["blockNum"] = blockNum,
	_["nPredNum"] = nPred,
        _["blockFac"] = IntegerMatrix(0),
	_["nPredFac"] = 0,
	_["nRow"] = blockNum.nrow(),
        _["facCard"] = IntegerVector(0),
	_["signature"] = signature
      );
  predBlock.attr("class") = "PredBlock";

  return predBlock;
}


/**
   @brief Unwraps field values useful for prediction.
 */
void PredblockUnwrap(SEXP sPredBlock, int &_nRow, int &_nPredNum, int &_nPredFac, NumericMatrix &_blockNum, IntegerMatrix &_blockFac) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");
  
  _nRow = as<int>((SEXP) predBlock["nRow"]);
  _nPredFac = as<int>((SEXP) predBlock["nPredFac"]);
  _nPredNum = as<int>((SEXP) predBlock["nPredNum"]);
  _blockNum = as<NumericMatrix>((SEXP) predBlock["blockNum"]);
  _blockFac = as<IntegerMatrix>((SEXP) predBlock["blockFac"]);
}


/**
   @brief Unwraps field values useful for export.
 */
void SignatureUnwrap(SEXP sSignature, IntegerVector &_predMap, List &_level) {
  List signature(sSignature);
  if (!signature.inherits("Signature"))
    stop("Expecting Signature");
  _predMap = as<IntegerVector>((SEXP) signature["predMap"]);
  _level = as<List>(signature["level"]);
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppPredict.cc

   @brief C++ interface to R entry for prediction methods.

   @author Mark Seligman
 */

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#include "rcppPredblock.h"
#include "rcppForest.h"
#include "rcppLeaf.h"
#include "predict.h"
#include "forest.h"
#include "leaf.h"

#include <algorithm>
//#include <iostream>


/**
   @brief Utility for computing mean-square error of prediction.
   
   @param yValid is the test vector.

   @param y is the vector of predictions.

   @param rsq outputs the r-squared statistic.

   @return mean-square error, with output parameter.
 */
double MSE(const double yValid[], NumericVector y, double &rsq) {
  double sse = 0.0;
  for (int i = 0; i < y.length(); i++) {
    double error = yValid[i] - y[i];
    sse += error * error;
  }

  // TODO:  Repair assumption that every row sampled.

  double mse = sse / y.length();
  rsq = 1.0 - (mse * y.length()) / (var(y) * (y.length() - 1.0));

  return mse;
}


/**
   @brief Predction for regression.

   @return Wrapped zero, with copy-out parameters.
 */
RcppExport SEXP RcppPredictReg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest, bool bag) {
  int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  PredblockUnwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac);

  std::vector<unsigned int> origin, facOrig, facSplit;
  std::vector<ForestNode> *forestNode;
  ForestUnwrap(sForest, origin, facOrig, facSplit, forestNode);
  
  std::vector<double> yRanked;
  std::vector<unsigned int> leafOrigin;
  std::vector<LeafNode> *leafNode;
  std::vector<BagRow> *bagRow;
  unsigned int rowTrain;
  std::vector<unsigned int> rank;
  LeafUnwrapReg(sLeaf, yRanked, leafOrigin, leafNode, bagRow, rowTrain, rank);

  std::vector<double> yPred(nRow);
  Predict::Regression(nPredNum > 0 ? transpose(blockNum).begin() : 0, nPredFac > 0 ? transpose(blockFac).begin() : 0, nPredNum, nPredFac, *forestNode, origin, facOrig, facSplit, leafOrigin, *leafNode, *bagRow, rank, yPred, bag ? rowTrain : 0);

  List prediction;
  if (Rf_isNull(sYTest)) { // Prediction
    prediction = List::create(
			 _["yPred"] = yPred,
			 _["qPred"] = NumericMatrix(0)
		     );
    prediction.attr("class") = "PredictReg";
  }
  else { // Validation
    NumericVector yTest(sYTest);
    double rsq;
    double mse = MSE(&yPred[0], yTest, rsq);
    prediction = List::create(
			 _["yPred"] = yPred,
			 _["mse "]= mse,
			 _["rsq"] = rsq,
			 _["qPred"] = NumericMatrix(0)
		     );
    prediction.attr("class") = "ValidReg";
  }

  return prediction;
}


RcppExport SEXP RcppValidateReg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictReg(sPredBlock, sForest, sLeaf, sYTest, true);
}


RcppExport SEXP RcppTestReg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictReg(sPredBlock, sForest, sLeaf, sYTest, false);
}


/**
   @brief Prediction for classification.

   @return Prediction list.
 */
RcppExport SEXP RcppPredictCtg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest, bool bag, bool doProb) {
  int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  PredblockUnwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac);
    
  std::vector<unsigned int> origin, facOrig, facSplit;
  std::vector<ForestNode> *forestNode;
  ForestUnwrap(sForest, origin, facOrig, facSplit, forestNode);

  std::vector<unsigned int> leafOrigin;
  std::vector<LeafNode> *leafNode;
  std::vector<BagRow> *bagRow;
  unsigned int rowTrain;
  std::vector<double> weight;
  CharacterVector levelsTrain;
  LeafUnwrapCtg(sLeaf, leafOrigin, leafNode, bagRow, rowTrain, weight, levelsTrain);

  unsigned int ctgWidth = levelsTrain.length();
  bool validate = !Rf_isNull(sYTest);
  IntegerVector yTest = validate ? IntegerVector(sYTest) - 1 : IntegerVector(0);
  CharacterVector levelsTest = validate ? as<CharacterVector>(IntegerVector(sYTest).attr("levels")) : CharacterVector(0);
  IntegerVector levelMatch = validate ? match(levelsTest, levelsTrain) : IntegerVector(0);
  unsigned int testWidth;
  unsigned int testLength = yTest.length();
  bool dimFixup = false;
  if (validate) {
    if (is_true(any(levelsTest != levelsTrain))) {
      dimFixup = true;
      IntegerVector sq = seq(0, levelsTest.length() - 1);
      IntegerVector idxNonMatch = sq[is_na(levelMatch)];
      if (idxNonMatch.length() > 0) {
	warning("Unreachable test levels not encountered in training");
	int proxy = ctgWidth + 1;
	for (int i = 0; i < idxNonMatch.length(); i++) {
	  int idx = idxNonMatch[i];
	  levelMatch[idx] = proxy++;
        }
      }

    // Matches are one-based.
      for (unsigned int i = 0; i < testLength; i++) {
        yTest[i] = levelMatch[yTest[i]] - 1;
      }
      testWidth = max(yTest) + 1;
    }
    else {
      testWidth = levelsTest.length();
    }
  }
  else {
    testWidth = 0;
  }
  std::vector<unsigned int> testCore(testLength);
  for (unsigned int i = 0; i < testLength; i++) {
    testCore[i] = yTest[i];
  }

  IntegerVector confCore(testWidth * ctgWidth);
  std::vector<double> misPredCore(testWidth);
  IntegerVector censusCore = IntegerVector(nRow * ctgWidth);
  std::vector<int> yPred(nRow);
  NumericVector probCore = doProb ? NumericVector(nRow * ctgWidth) : NumericVector(0);
  Predict::Classification(nPredNum > 0 ? transpose(blockNum).begin() : 0, nPredFac > 0 ? transpose(blockFac).begin() : 0, nPredNum, nPredFac, *forestNode, origin, facOrig, facSplit, leafOrigin, *leafNode, *bagRow, weight, yPred, censusCore.begin(), testCore, validate ? confCore.begin() : 0, misPredCore, doProb ? probCore.begin() : 0, bag ? rowTrain : 0);

  List predBlock(sPredBlock);
  IntegerMatrix census = transpose(IntegerMatrix(ctgWidth, nRow, censusCore.begin()));
  census.attr("dimnames") = List::create(predBlock["rowNames"], levelsTrain);
  NumericMatrix prob = doProb ? transpose(NumericMatrix(ctgWidth, nRow, probCore.begin())) : NumericMatrix(0);
  if (doProb) {
    prob.attr("dimnames") = List::create(predBlock["rowNames"], levelsTrain);
  }

  for (int i = 0; i < nRow; i++) // Bases to unity for front end.
    yPred[i] = yPred[i] + 1;

  List prediction;
  if (validate) {
    IntegerMatrix conf = transpose(IntegerMatrix(ctgWidth, testWidth, confCore.begin()));
    NumericVector misPred(levelsTest.length());
    if (dimFixup) {
      IntegerMatrix confOut(levelsTest.length(), ctgWidth);
      for (int i = 0; i < levelsTest.length(); i++) {
        confOut(i, _) = conf(levelMatch[i] - 1, _);
	misPred[i] = misPredCore[levelMatch[i] - 1];
      }
      conf = confOut;
    }
    else {
      for (int i = 0; i < levelsTest.length(); i++)
	misPred[i] = misPredCore[i];
    }
    misPred.attr("names") = levelsTest;
    conf.attr("dimnames") = List::create(levelsTest, levelsTrain);
    prediction = List::create(
      _["misprediction"] = misPred,
      _["confusion"] = conf,
      _["yPred"] = yPred,
      _["census"] = census,
      _["prob"] = prob
    );
    prediction.attr("class") = "ValidCtg";
  }
  else {
    prediction = List::create(
      _["yPred"] = yPred,
      _["census"] = census,
      _["prob"] = prob
   );
   prediction.attr("class") = "PredictCtg";
  }

  return prediction;
}


RcppExport SEXP RcppValidateVotes(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictCtg(sPredBlock, sForest, sLeaf, sYTest, true, false);
}


RcppExport SEXP RcppValidateProb(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictCtg(sPredBlock, sForest, sLeaf, sYTest, true, true);
}


/**
   @brief Predicts with class votes.

   @param sPredBlock contains the blocked observations.

   @param sForest contains the trained forest.

   @param sLeaf contains the trained leaves.

   @param sVotes outputs the vote predictions.

   @return Prediction object.
 */
RcppExport SEXP RcppTestVotes(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictCtg(sPredBlock, sForest, sLeaf, sYTest, false, false);
}


/**
   @brief Predicts with class votes.

   @param sPredBlock contains the blocked observations.

   @param sForest contains the trained forest.

   @param sLeaf contains the trained leaves.

   @param sVotes outputs the vote predictions.

   @return Prediction object.
 */
RcppExport SEXP RcppTestProb(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest) {
  return RcppPredictCtg(sPredBlock, sForest, sLeaf, sYTest, false, true);
}


/**
   @brief Prediction with quantiles.

   @param sPredBlock contains the blocked observations.

   @param sForest contains the trained forest.

   @param sLeaf contains the trained leaves.

   @param sVotes outputs the vote predictions.

   @param sQuantVec is a vector of quantile training data.
   
   @param sQBin is the bin parameter.

   @param sYTest is the test vector.

   @param bag is true iff validating.

   @return Prediction list.
*/
RcppExport SEXP RcppPredictQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sQuantVec, SEXP sQBin, SEXP sYTest, bool bag) {
  int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  PredblockUnwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac);
    
  std::vector<unsigned int> origin, facOrig, facSplit;
  std::vector<ForestNode> *forestNode;
  ForestUnwrap(sForest, origin, facOrig, facSplit, forestNode);

  std::vector<double> yRanked;
  std::vector<unsigned int> leafOrigin;
  std::vector<LeafNode> *leafNode;
  std::vector<BagRow> *bagRow;
  unsigned int rowTrain;
  std::vector<unsigned int> rank;
  LeafUnwrapReg(sLeaf, yRanked, leafOrigin, leafNode, bagRow, rowTrain, rank);

  std::vector<double> yPred(nRow);
  std::vector<double> quantVecCore(as<std::vector<double> >(sQuantVec));
  std::vector<double> qPredCore(nRow * quantVecCore.size());
  Predict::Quantiles(nPredNum > 0 ? transpose(blockNum).begin() : 0, nPredFac > 0 ? transpose(blockFac).begin() : 0, nPredNum, nPredFac, *forestNode, origin, facOrig, facSplit, leafOrigin, *leafNode, *bagRow, rank, yRanked, yPred, quantVecCore, as<int>(sQBin), qPredCore,  bag ? rowTrain : 0);

  NumericMatrix qPred(transpose(NumericMatrix(quantVecCore.size(), nRow, qPredCore.begin())));
  List prediction;
  if (!Rf_isNull(sYTest)) {
    double rsq;
    double mse = MSE(&yPred[0], NumericVector(sYTest), rsq);
    prediction = List::create(
 	 _["yPred"] = yPred,
	 _["qPred"] = qPred,
	 _["mse"] = mse,
	 _["rsq"] = rsq
	  );
    prediction.attr("class") = "ValidReg";
  }
  else {
    prediction = List::create(
		 _["yPred"] = yPred,
		 _["qPred"] = qPred
	     );
    prediction.attr("class") = "PredictReg";
  }

  return prediction;
}


RcppExport SEXP RcppValidateQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sQuantVec, SEXP sQBin, SEXP sYTest) {
  return RcppPredictQuant(sPredBlock, sForest, sLeaf, sQuantVec, sQBin, sYTest, true);
}


RcppExport SEXP RcppTestQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sQuantVec, SEXP sQBin, SEXP sYTest) {
  return RcppPredictQuant(sPredBlock, sForest, sLeaf, sQuantVec, sQBin, sYTest, false);
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppRowrank.cc

   @brief C++ interface to R entries for maintaining predictor data structures.

   @author Mark Seligman
*/

#include <Rcpp.h>

#include "rowrank.h"
// Testing only:
//#include <iostream>


using namespace Rcpp;
using namespace std;


/**
   @brief Builds row/rank maps as parallel arrays.

   @param sPredBlock is an (S3) PredBlock object.

   @return parallel row and rank arrays and the inverse numeric mapping.
 */
RcppExport SEXP RcppRowRank(SEXP sPredBlock) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");

  unsigned int nRow = as<unsigned int>(predBlock["nRow"]);
  unsigned int nPredNum = as<unsigned int>(predBlock["nPredNum"]);
  unsigned int nPredFac = as<unsigned int>(predBlock["nPredFac"]);
  unsigned int nPred = nPredNum + nPredFac;
  IntegerVector rank = IntegerVector(nRow * nPred);
  IntegerVector row = IntegerVector(nRow * nPred);
  IntegerVector invNum = 0;
  if (nPredNum > 0) {
    invNum = IntegerVector(nRow * nPredNum);
    NumericMatrix blockNum(as<NumericMatrix>(predBlock["blockNum"]));
    RowRank::PreSortNum(blockNum.begin(), nPredNum, nRow, row.begin(), rank.begin(), invNum.begin());
  }
  if (nPredFac > 0) {
    IntegerMatrix blockFac(as<IntegerMatrix>(predBlock["blockFac"]));
    RowRank::PreSortFac(blockFac.begin(), nPredNum, nPredFac, nRow, row.begin(), rank.begin());
  }
  
  List rowRank = List::create(
      _["row"] = row,			      
      _["rank"] = rank,
      _["invNum"] = invNum
    );
  rowRank.attr("class") = "RowRank";

  return rowRank;
}
// Copyright (C)  2012-2016   Mark Seligman
//
// This file is part of ArboristBridgeR.
//
// ArboristBridgeR is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// ArboristBridgeR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ArboristBridgeR.  If not, see <http://www.gnu.org/licenses/>.

/**
   @file rcppTrain.cc

   @brief C++ interface to R entry for training.

   @author Mark Seligman
 */

#include <Rcpp.h>
using namespace Rcpp;

#include "rcppForest.h"
#include "rcppLeaf.h"
#include "train.h"
#include "forest.h"
#include "leaf.h"

//#include <iostream>
using namespace std;


/**
   @brief R-language interface to response caching.

   @parm sY is the response vector.

   @return Wrapped value of response cardinality, if applicable.
 */
void RcppProxyCtg(IntegerVector y, NumericVector classWeight, std::vector<double> &proxy) {
  // Class weighting constructs a proxy response from category frequency.
  // The response is then jittered to diminish the possibility of ties
  // during scoring.  The magnitude of the jitter, then, should be scaled
  // so that no combination of samples can "vote" themselves into a
  // false plurality.
  //
  if (is_true(all(classWeight == 0.0))) { // Place-holder for balancing.
    NumericVector tb(table(y));
    for (unsigned int i = 0; i < classWeight.length(); i++) {
      classWeight[i] = tb[i] == 0.0 ? 0.0 : 1.0 / tb[i];
    }
  }
  classWeight = classWeight / sum(classWeight);

  unsigned int nRow = y.length();
  double recipLen = 1.0 / nRow;
  NumericVector yWeighted = classWeight[y];
  RNGScope scope;
  NumericVector rn(runif(nRow));
  for (unsigned int i = 0; i < nRow; i++)
    proxy[i] = yWeighted[i] + (rn[i] - 0.5) * 0.5 * (recipLen * recipLen);
}


/**
   @brief Constructs classification forest.

   @param sNTree is the number of trees requested.

   @param sMinNode is the smallest index node width allowed for splitting.

   @param sMinRatio is a threshold ratio of information measures between an index node and its offspring, below which the node does not split.

   @param sTotLevels is an upper bound on the number of levels to construct for each tree.

   @param sTotLevels is an upper bound on the number of levels to construct for each tree.

   @return Wrapped length of forest vector, with output parameters.
 */
RcppExport SEXP RcppTrainCtg(SEXP sPredBlock, SEXP sRowRank, SEXP sYOneBased, SEXP sNTree, SEXP sNSamp, SEXP sSampleWeight, SEXP sWithRepl, SEXP sTrainBlock, SEXP sMinNode, SEXP sMinRatio, SEXP sTotLevels, SEXP sPredFixed, SEXP sProbVec, SEXP sClassWeight) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");

  List rowRank(sRowRank);
  if (!rowRank.inherits("RowRank"))
    stop("Expecting RowRank");

  unsigned int nRow = as<unsigned int>(predBlock["nRow"]);
  int nPredNum = as<int>(predBlock["nPredNum"]);
  int nPredFac = as<int>(predBlock["nPredFac"]);
  double *feNum = 0;
  int *feInvNum = 0;
  if (nPredNum > 0) {
    NumericMatrix xNum(as<NumericMatrix>(predBlock["blockNum"]));
    feNum = xNum.begin();
    IntegerVector invNum(as<IntegerVector>(rowRank["invNum"]));
    feInvNum = invNum.begin();
  }
  int *feFacCard = 0;
  int cardMax = 0;
  if (nPredFac > 0) {
    IntegerVector facCard(as<IntegerVector>(predBlock["facCard"]));
    feFacCard = facCard.begin();
    cardMax = max(facCard);
  }
  List signature(as<List>(predBlock["signature"]));
  IntegerVector predMap(as<IntegerVector>(signature["predMap"]));

  IntegerVector feRow(as<IntegerVector>(rowRank["row"]));
  IntegerVector feRank(as<IntegerVector>(rowRank["rank"]));

  IntegerVector yOneBased(sYOneBased);
  CharacterVector levels(yOneBased.attr("levels"));
  unsigned int ctgWidth = levels.length();

  IntegerVector y = yOneBased - 1;
  std::vector<double> proxy(y.length());
  NumericVector classWeight(as<NumericVector>(sClassWeight));
  RcppProxyCtg(y, classWeight, proxy);

  int nTree = as<int>(sNTree);
  NumericVector sampleWeight(as<NumericVector>(sSampleWeight));

  int nPred = nPredNum + nPredFac;
  NumericVector predProb = NumericVector(sProbVec)[predMap];

  Train::Init(feNum, feFacCard, cardMax, nPredNum, nPredFac, nRow, nTree, as<int>(sNSamp), sampleWeight.begin(), as<bool>(sWithRepl), as<int>(sTrainBlock), as<int>(sMinNode), as<double>(sMinRatio), as<int>(sTotLevels), ctgWidth, as<int>(sPredFixed), predProb.begin());

  std::vector<unsigned int> origin(nTree);
  std::vector<unsigned int> facOrig(nTree);
  std::vector<unsigned int> leafOrigin(nTree);
  NumericVector predInfo(nPred);

  std::vector<ForestNode> forestNode;

  std::vector<unsigned int> facSplit;
  std::vector<LeafNode> leafNode;
  std::vector<BagRow> bagRow;
  std::vector<double> weight;

  Train::Classification(feRow.begin(), feRank.begin(), feInvNum, as<std::vector<unsigned int> >(y), ctgWidth, proxy, origin, facOrig, predInfo.begin(), forestNode, facSplit, leafOrigin, leafNode, bagRow, weight);

  return List::create(
      _["forest"] = ForestWrap(origin, facOrig, facSplit, forestNode),
      _["leaf"] = LeafWrapCtg(leafOrigin, leafNode, bagRow, nRow, weight, CharacterVector(yOneBased.attr("levels"))),
      _["predInfo"] = predInfo[predMap] // Maps back from core order.
  );
}


RcppExport SEXP RcppTrainReg(SEXP sPredBlock, SEXP sRowRank, SEXP sY, SEXP sNTree, SEXP sNSamp, SEXP sSampleWeight, SEXP sWithRepl, SEXP sTrainBlock, SEXP sMinNode, SEXP sMinRatio, SEXP sTotLevels, SEXP sPredFixed, SEXP sProbVec, SEXP sRegMono) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");

  List rowRank(sRowRank);
  if (!rowRank.inherits("RowRank"))
    stop("Expecting RowRank");

  unsigned int nRow = as<unsigned int>(predBlock["nRow"]);
  int nPredNum = as<int>(predBlock["nPredNum"]);
  int nPredFac = as<int>(predBlock["nPredFac"]);
  double *feNum = 0;
  int *feInvNum = 0;
  if (nPredNum > 0) {
    NumericMatrix xNum(as<NumericMatrix>(predBlock["blockNum"]));
    feNum = xNum.begin();
    IntegerVector invNum(as<IntegerVector>(rowRank["invNum"]));
    feInvNum = invNum.begin();
  }
  int *feFacCard = 0;
  int cardMax = 0;
  if (nPredFac > 0) {
    IntegerVector facCard(as<IntegerVector>(predBlock["facCard"]));
    feFacCard = facCard.begin();
    cardMax = max(facCard);
  }
  List signature(as<List>(predBlock["signature"]));
  IntegerVector predMap(as<IntegerVector>(signature["predMap"]));
  
  int nTree = as<int>(sNTree);
  NumericVector sampleWeight(as<NumericVector>(sSampleWeight));

  int nPred = nPredNum + nPredFac;
  NumericVector predProb = NumericVector(sProbVec)[predMap];
  NumericVector regMono = NumericVector(sRegMono)[predMap];
  
  Train::Init(feNum, feFacCard, cardMax, nPredNum, nPredFac, nRow, nTree, as<int>(sNSamp), sampleWeight.begin(), as<bool>(sWithRepl), as<int>(sTrainBlock), as<int>(sMinNode), as<double>(sMinRatio), as<int>(sTotLevels), 0, as<int>(sPredFixed), predProb.begin(), regMono.begin());

  IntegerVector feRow(as<IntegerVector>(rowRank["row"]));
  IntegerVector feRank(as<IntegerVector>(rowRank["rank"]));

  NumericVector y(sY);
  NumericVector yRanked = clone(y).sort();
  IntegerVector row2Rank = match(y, yRanked) - 1;

  std::vector<unsigned int> origin(nTree);
  std::vector<unsigned int> facOrig(nTree);
  std::vector<unsigned int> leafOrigin(nTree);
  NumericVector predInfo(nPred);

  std::vector<ForestNode> forestNode;
  std::vector<LeafNode> leafNode;
  std::vector<BagRow> bagRow;
  std::vector<unsigned int> rank;
  std::vector<unsigned int> facSplit;

  Train::Regression(feRow.begin(), feRank.begin(), feInvNum, as<std::vector<double> >(y), as<std::vector<unsigned int> >(row2Rank), origin, facOrig, predInfo.begin(), forestNode, facSplit, leafOrigin, leafNode, bagRow, rank);

  return List::create(
      _["forest"] = ForestWrap(origin, facOrig, facSplit, forestNode),
      _["leaf"] = LeafWrapReg(leafOrigin, leafNode, bagRow, nRow, rank, as<std::vector<double> >(yRanked)),
      _["predInfo"] = predInfo[predMap] // Maps back from core order.
    );
}