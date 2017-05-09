// Copyright (C)  2012-2017   Mark Seligman
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
// Copyright (C)  2012-2017   Mark Seligman
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
   @file rcppSample.cc

   @brief Interface to front-end methods implementing (response) sampling.

   @author Mark Seligman
 */

#include <RcppArmadilloExtensions/sample.h>
#include "rcppSample.h"

//#include <iostream>
//using namespace std;

unsigned int RcppSample::nRow = 0;
bool RcppSample::withRepl = false;

NumericVector weightNull(0);
NumericVector &RcppSample::weight = weightNull;

/**
   @brief Caches row sampling parameters as static values.

   @param _nRow is length of the response vector.

   @param _weight is user-specified weighting of row samples.

   @param _withRepl is true iff sampling with replacement.

   @return void.
 */
void RcppSample::Init(unsigned int _nRow, const double feWeight[], bool _withRepl) {
  nRow = _nRow;
  NumericVector _weight(nRow);
  weight = _weight;
  for (unsigned int i = 0; i < nRow; i++)
    weight[i] = feWeight[i];
  withRepl = _withRepl;
}


/**
   @brief Samples row indices either with or without replacement using methods from RccpArmadillo.

   @param nSamp is the number of samples to draw.

   @param out[] is an output vector of sampled row indices.

   @return void, with output vector.
 */
void RcppSample::SampleRows(unsigned int nSamp, int out[]) {
  RNGScope scope;
  IntegerVector rowVec(seq_len(nRow)-1);
  IntegerVector samp = RcppArmadillo::sample(rowVec, nSamp, withRepl, weight);

  for (unsigned int i = 0; i < nSamp; i++)
    out[i] = samp[i];
}
// Copyright (C)  2012-2017   Mark Seligman
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

//using namespace std;
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
  unsigned int *nodeOrigin, *facOrigin, *facSplit;
  ForestNode *forestNode;
  unsigned int nTree, nFac, nodeEnd;
  size_t facLen;
  RcppForest::Unwrap(sForest, nodeOrigin, nTree, facSplit, facLen, facOrigin, nFac, forestNode, nodeEnd);

  std::vector<std::vector<unsigned int> > predTree(nTree), bumpTree(nTree);
  std::vector<std::vector<double > > splitTree(nTree);
  ForestNode::Export(nodeOrigin, nTree, forestNode, nodeEnd, predTree, bumpTree, splitTree);
  PredExport(predMap.begin(), predTree, bumpTree);
  
  std::vector<std::vector<unsigned int> > facSplitTree(nTree);
  BVJagged::Export(facSplit, facLen, facOrigin, nTree, facSplitTree);

  std::vector<double> yTrain;
  std::vector<unsigned int> leafOrigin;
  LeafNode *leafNode;
  unsigned int leafCount;
  BagLeaf *bagLeaf;
  unsigned int bagLeafTot;
  unsigned int *bagBits;
  RcppLeaf::UnwrapReg(sLeaf, yTrain, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, true);
  unsigned int rowTrain = yTrain.size();

  std::vector<std::vector<unsigned int> > rowTree(nTree), sCountTree(nTree);
  std::vector<std::vector<double> > scoreTree(nTree);
  std::vector<std::vector<unsigned int> > extentTree(nTree);
  LeafReg::Export(leafOrigin, leafNode, leafCount, bagLeaf, bagBits, rowTrain, rowTree, sCountTree, scoreTree, extentTree);

  List outBundle = List::create(
				_["rowTrain"] = rowTrain,
				_["pred"] = predTree,
				_["bump"] = bumpTree,
				_["split"] = splitTree,
				_["facSplit"] = facSplitTree,
				_["row"] = rowTree,
				_["sCount"] = sCountTree,
				_["score"] = scoreTree,
				_["extent"] = extentTree
				);
  outBundle.attr("class") = "ExportReg";

  return outBundle;
}


/**
   @brief Exports core data structures as vector of per-tree vectors.

   @return List with common and classification-specific members.
 */
RcppExport SEXP ExportCtg(SEXP sForest, SEXP sLeaf, IntegerVector predMap) {
  unsigned int *nodeOrigin, *facOrigin, *facSplit;
  ForestNode *forestNode;
  unsigned int nTree, nFac, nodeEnd;
  size_t facLen;
  RcppForest::Unwrap(sForest, nodeOrigin, nTree, facSplit, facLen, facOrigin, nFac, forestNode, nodeEnd);

  std::vector<std::vector<unsigned int> > predTree(nTree), bumpTree(nTree);
  std::vector<std::vector<double > > splitTree(nTree);
  ForestNode::Export(nodeOrigin, nTree, forestNode, nodeEnd, predTree, bumpTree, splitTree);
  PredExport(predMap.begin(), predTree, bumpTree);
  
  std::vector<std::vector<unsigned int> > facSplitTree(nTree);
  BVJagged::Export(facSplit, facLen, facOrigin, nTree, facSplitTree);

  std::vector<unsigned int> leafOrigin;
  LeafNode *leafNode;
  unsigned int leafCount;
  BagLeaf *bagLeaf;
  unsigned int bagLeafTot;
  unsigned int *bagBits;
  double *weight;
  unsigned int rowTrain;
  CharacterVector yLevel;
  RcppLeaf::UnwrapCtg(sLeaf, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, weight, rowTrain, yLevel, true);

  std::vector<std::vector<unsigned int> > rowTree(nTree), sCountTree(nTree);
  std::vector<std::vector<double> > scoreTree(nTree);
  std::vector<std::vector<unsigned int> > extentTree(nTree);
  std::vector<std::vector<double> > weightTree(nTree);
  LeafCtg::Export(leafOrigin, leafNode, leafCount, bagLeaf, bagBits, rowTrain, weight, yLevel.length(), rowTree, sCountTree, scoreTree, extentTree, weightTree);

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
  std::vector<std::vector<unsigned int> > rowTree = forestCore["row"];
  std::vector<std::vector<unsigned int> > sCountTree = forestCore["sCount"];
  IntegerVector row(rowTree[tIdx].begin(), rowTree[tIdx].end());
  IntegerVector sCount(sCountTree[tIdx].begin(), sCountTree[tIdx].end());
  IntegerVector bag = IntegerVector(as<unsigned int>(forestCore["rowTrain"]), 0);
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
  RcppPredblock::SignatureUnwrap(arbOut["signature"], predMap, predLevel);

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
// Copyright (C)  2012-2017   Mark Seligman
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


#include "forest.h"
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#include "rcppForest.h"

//#include <iostream>

SEXP RcppForest::Wrap(const std::vector<unsigned int> &origin, const std::vector<unsigned int> &facOrigin, const std::vector<unsigned int> &facSplit, const std::vector<ForestNode> &forestNode) {
  size_t forestSize = forestNode.size() * sizeof(ForestNode);
  RawVector forestRaw(forestSize);
  for (size_t i = 0; i < forestSize; i++) {
    forestRaw[i] = ((unsigned char *) &forestNode[0])[i];
  }

  size_t facSize = facSplit.size() * sizeof(unsigned int);
  RawVector facRaw(facSize);
  for (size_t i = 0; i < facSize; i++) {
    facRaw[i] = ((unsigned char*) &facSplit[0])[i];
  }
  
  List forest = List::create(
     _["forestNode"] = forestRaw,
     _["origin"] = origin,
     _["facOrig"] = facOrigin,
     _["facSplit"] = facRaw);
  forest.attr("class") = "Forest";

  return forest;
}

RawVector RcppForest::rv1 = RawVector(0);
RawVector RcppForest::rv2 = RawVector(0);
IntegerVector RcppForest::iv1 = IntegerVector(0);
IntegerVector RcppForest::iv2 = IntegerVector(0);

/**
   @brief Exposes front-end Forest fields for transmission to core.

   @return void.
 */
void RcppForest::Unwrap(SEXP sForest, unsigned int *&_origin, unsigned int &_nTree, unsigned int *&_facSplit, size_t &_facLen, unsigned int *&_facOrig, unsigned int &_nFac, ForestNode *&_forestNode, unsigned int &_nodeEnd) {
  List forest(sForest);
  if (!forest.inherits("Forest"))
    stop("Expecting Forest");

  // Alignment should be sufficient to guarantee safety of
  // the casted loads.
  //
  iv1 = IntegerVector((SEXP) forest["origin"]);
  _origin = (unsigned int*) &iv1[0];
  _nTree = iv1.length();

  rv1 = RawVector((SEXP) forest["facSplit"]);
  _facSplit = (unsigned int*) &rv1[0];
  _facLen = rv1.length() / sizeof(unsigned int);

  iv2 = IntegerVector((SEXP) forest["facOrig"]);
  _facOrig = (unsigned int*) &iv2[0];
  _nFac = iv2.length();

  rv2 = RawVector((SEXP) forest["forestNode"]);
  _forestNode = (ForestNode*) &rv2[0];
  _nodeEnd = rv2.length() / sizeof(ForestNode);
}


void RcppForest::Clear() {
  rv1 = RawVector(0);
  rv2 = RawVector(0);
  iv1 = IntegerVector(0);
  iv2 = IntegerVector(0);
}
// Copyright (C)  2012-2017  Mark Seligman
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
   @file rcppInit.cc

   @brief C++ interface to R entry for symbol registration.

   @author Mark Seligman
 */
void R_init_Rborist(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
// Copyright (C)  2012-2017   Mark Seligman
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


#include "leaf.h"
#include "rcppLeaf.h"


/**
   @brief Wraps core (regression) Leaf vectors for reference by front end.
 */
SEXP RcppLeaf::WrapReg(const std::vector<unsigned int> &leafOrigin, std::vector<LeafNode> &leafNode, const std::vector<BagLeaf> &bagLeaf, const std::vector<unsigned int> &bagBits, const std::vector<double> &yTrain) {
  RawVector leafRaw(leafNode.size() * sizeof(LeafNode));
  RawVector blRaw(bagLeaf.size() * sizeof(BagLeaf));
  RawVector bbRaw(bagBits.size() * sizeof(unsigned int));
  Serialize(leafNode, bagLeaf, bagBits, leafRaw, blRaw, bbRaw);
  List leaf = List::create(
   _["origin"] = leafOrigin,
   _["node"] = leafRaw,
   _["bagLeaf"] = blRaw,
   _["bagBits"] = bbRaw,
   _["yTrain"] = yTrain
  );
  leaf.attr("class") = "LeafReg";
  
  return leaf;
}


RawVector RcppLeaf::rv1 = RawVector(0);
RawVector RcppLeaf::rv2 = RawVector(0);
RawVector RcppLeaf::rv3 = RawVector(0);
NumericVector RcppLeaf::nv1 = NumericVector(0);

/**
   @brief Exposes front-end (regression) Leaf fields for transmission to core.

   @param sLeaf is the R object containing the leaf (list) data.

   @param _yTrain outputs the training response.

   @param _leafInfoReg outputs the sample counts, organized by leaf.

   @param bag indicates whether to include bagging information.

   @return void, with output reference parameters.
 */
void RcppLeaf::UnwrapReg(SEXP sLeaf, std::vector<double> &_yTrain, std::vector<unsigned int> &_leafOrigin, LeafNode *&_leafNode, unsigned int &_leafCount, BagLeaf *&_bagLeaf, unsigned int &_bagLeafTot, unsigned int *&_bagBits, bool bag) {
  List leaf(sLeaf);
  if (!leaf.inherits("LeafReg"))
    stop("Expecting LeafReg");

  rv1 = RawVector((SEXP) leaf["bagBits"]);
  _bagBits = bag ? (unsigned int *) &rv1[0] : 0;
  
  rv2 = RawVector((SEXP) leaf["bagLeaf"]);
  _bagLeaf = bag ? (BagLeaf *) &rv2[0] : 0;
  _bagLeafTot = bag ? rv2.length() / sizeof(BagLeaf) : 0;
  
  _leafOrigin = as<std::vector<unsigned int> >(leaf["origin"]);

  rv3 = RawVector((SEXP) leaf["node"]);
  _leafNode = (LeafNode*) &rv3[0];
  _leafCount = rv3.length() / sizeof(LeafNode);

  _yTrain = as<std::vector<double> >(leaf["yTrain"]);
}


/**
   @brief Wraps core (classification) Leaf vectors for reference by front end.
 */
SEXP RcppLeaf::WrapCtg(const std::vector<unsigned int> &leafOrigin, const std::vector<LeafNode> &leafNode, const std::vector<BagLeaf> &bagLeaf, const std::vector<unsigned int> &bagBits, const std::vector<double> &weight, unsigned int rowTrain, const CharacterVector &levels) {
  RawVector leafRaw(leafNode.size() * sizeof(LeafNode));
  RawVector blRaw(bagLeaf.size() * sizeof(BagLeaf));
  RawVector bbRaw(bagBits.size() * sizeof(unsigned int));
  Serialize(leafNode, bagLeaf, bagBits, leafRaw, blRaw, bbRaw);
  List leaf = List::create(
   _["origin"] = leafOrigin,	
   _["node"] = leafRaw,
   _["bagLeaf"] = blRaw,
   _["bagBits"] = bbRaw,
   _["weight"] = weight,
   _["rowTrain"] = rowTrain,
   _["levels"] = levels
   );
  leaf.attr("class") = "LeafCtg";

  return leaf;
}


/** 
    @brief Serializes the internally-typed objects, 'LeafNode', as well
    as the unsigned integer (packed bit) vector, "bagBits".
*/
void RcppLeaf::Serialize(const std::vector<LeafNode> &leafNode, const std::vector<BagLeaf> &bagLeaf, const std::vector<unsigned int> &bagBits, RawVector &leafRaw, RawVector &blRaw, RawVector &bbRaw) {
  for (size_t i = 0; i < leafNode.size() * sizeof(LeafNode); i++) {
    leafRaw[i] = ((unsigned char*) &leafNode[0])[i];
  }

  for (size_t i = 0; i < bagLeaf.size() * sizeof(BagLeaf); i++) {
    blRaw[i] = ((unsigned char*) &bagLeaf[0])[i];
  }

  for (size_t i = 0; i < bagBits.size() * sizeof(unsigned int); i++) {
    bbRaw[i] = ((unsigned char*) &bagBits[0])[i];
  }
}


/**
   @brief Exposes front-end (classification) Leaf fields for transmission to core.

   @param sLeaf is the R object containing the leaf (list) data.

   @param _weight outputs the sample weights.

   @param _levels outputs the category levels; retains as front-end object.

   @param bag indicates whether to include bagging information.

   @return void, with output reference parameters.
 */
void RcppLeaf::UnwrapCtg(SEXP sLeaf, std::vector<unsigned int> &_leafOrigin, LeafNode *&_leafNode, unsigned int &_leafCount, BagLeaf *&_bagLeaf, unsigned int &_bagLeafTot, unsigned int *&_bagBits, double *&_weight, unsigned int &_rowTrain, CharacterVector &_levels, bool bag) {
  List leaf(sLeaf);
  if (!leaf.inherits("LeafCtg")) {
    stop("Expecting LeafCtg");
  }

  rv1 = RawVector((SEXP) leaf["bagBits"]);
  _bagBits = bag ? (unsigned int *) &rv1[0] : 0;
  
  rv2 = RawVector((SEXP) leaf["bagLeaf"]);
  _bagLeaf = bag ? (BagLeaf *) &rv2[0] : 0;
  _bagLeafTot = bag ? rv2.length() / sizeof(BagLeaf) : 0;

  _leafOrigin = as<std::vector<unsigned int> >(leaf["origin"]);

  rv3 = RawVector((SEXP) leaf["node"]);
  _leafNode = (LeafNode*) &rv3[0];
  _leafCount = rv3.length() / sizeof(LeafNode);

  nv1 = NumericVector((SEXP) leaf["weight"]);
  _weight = &nv1[0];

  _rowTrain = as<unsigned int>((SEXP) leaf["rowTrain"]);
  _levels = as<CharacterVector>((SEXP) leaf["levels"]);
}


void RcppLeaf::Clear() {
  rv1 = RawVector(0);
  rv2 = RawVector(0);
  rv3 = RawVector(0);
  nv1 = NumericVector(0);
}
// Copyright (C)  2012-2017   Mark Seligman
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

// Testing only:
//#include <iostream>
//using namespace std;

#include "rcppPredblock.h"
#include "rowrank.h"


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
  std::vector<unsigned int> levels = as<std::vector<unsigned int> >(sLevels);
  unsigned int nRow = xf.nrows();
  unsigned int nPredFac = facElt.length();
  unsigned int nPredNum = numElt.length();
  unsigned int nPred = nPredFac + nPredNum;

  IntegerVector predMap(nPred);
  IntegerVector facCard(0);
  IntegerMatrix xFac;
  NumericMatrix xNum;
  if (nPredNum > 0) {
    xNum = NumericMatrix(nRow, nPredNum);
  }
  else
    xNum = NumericMatrix(0, 0);
  if (nPredFac > 0) {
    facCard = IntegerVector(nPredFac);
    xFac = IntegerMatrix(nRow, nPredFac);
  }
  else {
    xFac = IntegerMatrix(0);
  }

  int numIdx = 0;
  int facIdx = 0;
  List level(nPredFac);
  for (unsigned int feIdx = 0; feIdx < nPred; feIdx++) {
    unsigned int card = levels[feIdx];
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
    RcppPredblock::FactorRemap(xFac, level, levelTrain);
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
      _["blockNumRLE"] = R_NilValue, // For now.
      _["blockFacRLE"] = R_NilValue, // For now.
      _["blockFac"] = xFac,
      _["nPredFac"] = nPredFac,
      _["nRow"] = nRow,
      _["facCard"] = facCard,
      _["signature"] = signature
      );
  predBlock.attr("class") = "PredBlock";

  return predBlock;
}


void RcppPredblock::FactorRemap(IntegerMatrix &xFac, List &levelTest, List &levelTrain) {
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

  IntegerVector facCard(0);
  List predBlock = List::create(
	_["colNames"] = colnames(blockNum),
	_["rowNames"] = rownames(blockNum),
	_["blockNum"] = blockNum,
	_["blockNumRLE"] = R_NilValue, // For now.
	_["blockFacRLE"] = R_NilValue, // For now.
	_["nPredNum"] = nPred,
        _["blockFac"] = IntegerMatrix(0),
	_["nPredFac"] = 0,
	_["nRow"] = blockNum.nrow(),
        _["facCard"] = facCard,
	_["signature"] = signature
      );
  predBlock.attr("class") = "PredBlock";

  return predBlock;
}


/**
   @brief Reads an S4 object containing (sparse) dgCMatrix.
 */
RcppExport SEXP RcppPredBlockSparse(SEXP sX) {
  S4 spNum(sX);

  IntegerVector i;
  if (R_has_slot(sX, Rf_mkString("i"))) {
    i = spNum.slot("i");
  }
  IntegerVector j;
  if (R_has_slot(sX, Rf_mkString("j"))) {
    j = spNum.slot("j");
  }
  IntegerVector p;
  if (R_has_slot(sX, Rf_mkString("p"))) {
    p = spNum.slot("p");
  }

  if (!R_has_slot(sX, Rf_mkString("Dim"))) {
    stop("Expecting dimension slot");
  }
  IntegerVector dim = spNum.slot("Dim");
  unsigned int nRow = dim[0];
  unsigned int nPred = dim[1];

  // 'eltsNZ' holds the nonzero elements.
  NumericVector eltsNZ;
  if (R_has_slot(sX, Rf_mkString("x"))) {
    eltsNZ = spNum.slot("x");
  }
  else {
    stop("Pattern matrix:  NYI");
  }


  std::vector<double> valNum;
  std::vector<unsigned int> rowStart;
  std::vector<unsigned int> runLength;
  std::vector<unsigned int> predStart;

  // Divines the encoding format and packs appropriately.
  //
  if (i.length() == 0) {
    RcppPredblock::SparseJP(eltsNZ, j, p, nRow, valNum, rowStart, runLength);
  }
  else if (j.length() == 0) {
    RcppPredblock::SparseIP(eltsNZ, i, p, nRow, nPred, valNum, rowStart, runLength, predStart);
  }
  else if (p.length() == 0) {
    RcppPredblock::SparseIJ(eltsNZ, i, j, nRow, valNum, rowStart, runLength);
   }
  else {
    stop("Indeterminate sparse matrix format");
  }

  List blockNumRLE = List::create(
	  _["valNum"] = valNum,
	  _["rowStart"] = rowStart,
	  _["runLength"] = runLength,
	  _["predStart"] = predStart);
  blockNumRLE.attr("class") = "BlockNumRLE";

  List dimNames;
  CharacterVector rowName, colName;
  if (R_has_slot(sX, Rf_mkString("Dimnames"))) {
    dimNames = spNum.slot("Dimnames");
    if (!Rf_isNull(dimNames[0])) {
      rowName = dimNames[0];
    }
    if (!Rf_isNull(dimNames[1])) {
      colName = dimNames[1];
    }
  }

  List signature = List::create(
      _["predMap"] = seq_len(nPred) - 1,
      _["level"] = List::create(0)
  );
  signature.attr("class") = "Signature";
  IntegerVector facCard(0);

  List predBlock = List::create(
	_["colNames"] = colName,
	_["rowNames"] = rowName,
	_["blockNum"] = NumericMatrix(0),
	_["nPredNum"] = nPred,
	_["blockNumRLE"] = blockNumRLE,
	_["blockFacRLE"] = R_NilValue, // For now.
        _["blockFac"] = IntegerMatrix(0),
	_["nPredFac"] = 0,
	_["nRow"] = nRow,
        _["facCard"] = facCard,
	_["signature"] = signature
      );

  predBlock.attr("class") = "PredBlock";

  return predBlock;
}


// 'i' in [0, nRow-1] list rows with nonzero elements.
// 'p' holds the starting offset for each column in 'eltsNZ'.
//    Repeated values indicate full-zero columns. 
//
void RcppPredblock::SparseIP(const NumericVector &eltsNZ, const IntegerVector &i, const IntegerVector &p, unsigned int nRow, unsigned int nCol, std::vector<double> &valNum, std::vector<unsigned int> &rowStart, std::vector<unsigned int> &runLength, std::vector<unsigned int> &predStart) {
  // Pre-scans column heights. 'p' has length one greater than number
  // of columns, providing ready access to heights.
  const double zero = 0.0;
  std::vector<unsigned int> nzHeight(p.length());
  unsigned int idxStart = p[0];
  for (R_len_t colIdx = 1; colIdx < p.length(); colIdx++) {
    nzHeight[colIdx - 1] = p[colIdx] - idxStart;
    idxStart = p[colIdx];
  }
  
  for (unsigned int colIdx = 0; colIdx < nCol; colIdx++) {
    unsigned int height = nzHeight[colIdx];
    predStart.push_back(valNum.size());
    if (height == 0) {
      valNum.push_back(zero);
      runLength.push_back(nRow);
      rowStart.push_back(0);
    }
    else {
      unsigned int nzPrev = nRow; // Inattainable row value.
      // Row indices into 'i' and 'x' are zero-based.
      unsigned int idxStart = p[colIdx];
      unsigned int idxEnd = idxStart + height;
      for (unsigned int rowIdx = idxStart; rowIdx < idxEnd; rowIdx++) {
        unsigned int nzRow = i[rowIdx];
        if (nzPrev == nRow && nzRow > 0) { // Zeroes lead.
	  valNum.push_back(zero);
	  runLength.push_back(nzRow);
	  rowStart.push_back(0);
	}
	else if (nzRow > nzPrev + 1) { // Zeroes precede.
	  valNum.push_back(zero);
	  runLength.push_back(nzRow - (nzPrev + 1));
	  rowStart.push_back(nzPrev + 1);
	}
	valNum.push_back(eltsNZ[rowIdx]);
	runLength.push_back(1);
	rowStart.push_back(nzRow);
	nzPrev = nzRow;
      }
      if (nzPrev + 1 < nRow) { // Zeroes trail.
	valNum.push_back(zero);
	runLength.push_back(nRow - nzPrev - 1);
	rowStart.push_back(nzPrev + 1);
      }
    }
  }
}



void RcppPredblock::SparseJP(NumericVector &eltsNZ, IntegerVector &j, IntegerVector &p, unsigned int nRow, std::vector<double> &valNum, std::vector<unsigned int> &rowStart, std::vector<unsigned int> &runLength) {
  stop("Sparse form j/p:  NYI");
}


    // 'i' holds row indices of nonzero elements.
    // 'j' " column " "
void RcppPredblock::SparseIJ(NumericVector &eltsNZ, IntegerVector &i, IntegerVector &j, unsigned int nRow, std::vector<double> &valNum, std::vector<unsigned int> &rowStart, std::vector<unsigned int> &runLength) {
  stop("Sparse form i/j:  NYI");
}


/**
   @brief Unwraps field values useful for prediction.
 */
void RcppPredblock::Unwrap(SEXP sPredBlock, unsigned int &_nRow, unsigned int &_nPredNum, unsigned int &_nPredFac, NumericMatrix &_blockNum, IntegerMatrix &_blockFac, std::vector<double> &_valNum, std::vector<unsigned int> &_rowStart, std::vector<unsigned int> &_runLength, std::vector<unsigned int> &_predStart) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");
  
  _nRow = as<unsigned int>((SEXP) predBlock["nRow"]);
  _nPredFac = as<unsigned int>((SEXP) predBlock["nPredFac"]);
  _nPredNum = as<unsigned int>((SEXP) predBlock["nPredNum"]);
  if (!Rf_isNull(predBlock["blockNumRLE"])) {
    List blockNumRLE((SEXP) predBlock["blockNumRLE"]);
    _valNum = as<std::vector<double> >((SEXP) blockNumRLE["valNum"]);
    _rowStart = as<std::vector<unsigned int> >((SEXP) blockNumRLE["rowStart"]);
    _runLength = as<std::vector<unsigned int> >((SEXP) blockNumRLE["runLength"]);
    _predStart = as<std::vector<unsigned int> >((SEXP) blockNumRLE["predStart"]);
  }
  else {
    _blockNum = as<NumericMatrix>((SEXP) predBlock["blockNum"]);
  }

  if (!Rf_isNull(predBlock["blockFacRLE"])) {
    stop("Sparse factors:  NYI");
  }
  else {
    _blockFac = as<IntegerMatrix>((SEXP) predBlock["blockFac"]);
  }
}


/**
   @brief Unwraps field values useful for export.
 */
void RcppPredblock::SignatureUnwrap(SEXP sSignature, IntegerVector &_predMap, List &_level) {
  List signature(sSignature);
  if (!signature.inherits("Signature"))
    stop("Expecting Signature");
  _predMap = as<IntegerVector>((SEXP) signature["predMap"]);
  _level = as<List>(signature["level"]);
}
// Copyright (C)  2012-2017  Mark Seligman
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

using namespace Rcpp;

#include "rcppPredblock.h"
#include "rcppForest.h"
#include "rcppLeaf.h"
#include "predict.h"

#include "forest.h"
#include "leaf.h"

#include <algorithm>

//#include <iostream>
//using namespace std;

/**
   @brief Utility for computing mean-square error of prediction.
   
   @param yValid is the test vector.

   @param y is the vector of predictions.

   @param rsq outputs the r-squared statistic.

   @return mean squared error, with output parameter.
 */
double MSE(const double yValid[], NumericVector y, double &rsq, double &mae) {
  double sse = 0.0;
  mae = 0.0;
  for (R_len_t i = 0; i < y.length(); i++) {
    double error = yValid[i] - y[i];
    sse += error * error;
    mae += abs(error);
  }
  rsq = 1.0 - sse / (var(y) * (y.length() - 1.0));
  mae /= y.length();

  return sse / y.length();
}


/**
   @brief Predction for regression.

   @return Wrapped zero, with copy-out parameters.
 */
RcppExport SEXP RcppPredictReg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest, bool validate) {
  unsigned int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  std::vector<double> valNum;
  std::vector<unsigned int> rowStart;
  std::vector<unsigned int> runLength;
  std::vector<unsigned int> predStart;
  RcppPredblock::Unwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac, valNum, rowStart, runLength, predStart);

  unsigned int *origin, *facOrig, *facSplit;
  ForestNode *forestNode;
  unsigned int nTree, nFac, nodeEnd;
  size_t facLen;
  RcppForest::Unwrap(sForest, origin, nTree, facSplit, facLen, facOrig, nFac, forestNode, nodeEnd);
  
  std::vector<double> yTrain;
  std::vector<unsigned int> leafOrigin;
  LeafNode *leafNode;
  unsigned int leafCount;
  BagLeaf *bagLeaf;
  unsigned int bagLeafTot;
  unsigned int *bagBits;
  RcppLeaf::UnwrapReg(sLeaf, yTrain, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, validate);

  std::vector<double> yPred(nRow);
  Predict::Regression(valNum, rowStart, runLength, predStart, (valNum.size() == 0 && nPredNum > 0) ? transpose(blockNum).begin() : 0, nPredFac > 0 ? (unsigned int *) transpose(blockFac).begin() : 0, nPredNum, nPredFac, forestNode, origin, nTree, facSplit, facLen, facOrig, nFac, leafOrigin, leafNode, leafCount, bagBits, yTrain, yPred);

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
    double rsq, mae;
    double mse = MSE(&yPred[0], yTest, rsq, mae);
    prediction = List::create(
			 _["yPred"] = yPred,
			 _["mse"] = mse,
			 _["mae"] = mae,
			 _["rsq"] = rsq,
			 _["qPred"] = NumericMatrix(0)
		     );
    prediction.attr("class") = "ValidReg";
  }
  RcppLeaf::Clear();
  RcppForest::Clear();

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
RcppExport SEXP RcppPredictCtg(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest, bool validate, bool doProb) {
  unsigned int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  std::vector<double> valNum;
  std::vector<unsigned int> rowStart;
  std::vector<unsigned int> runLength;
  std::vector<unsigned int> predStart;
  RcppPredblock::Unwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac, valNum, rowStart, runLength, predStart);
    
  unsigned int *origin, *facOrig, *facSplit;
  ForestNode *forestNode;
  unsigned int nTree, nFac, nodeEnd;
  size_t facLen;
  RcppForest::Unwrap(sForest, origin, nTree, facSplit, facLen, facOrig, nFac, forestNode, nodeEnd);

  std::vector<unsigned int> leafOrigin;
  LeafNode *leafNode;
  unsigned int leafCount;
  BagLeaf *bagLeaf;
  unsigned int bagLeafTot;
  unsigned int *bagBits;
  double *weight;
  unsigned int rowTrain;
  CharacterVector levelsTrain;
  RcppLeaf::UnwrapCtg(sLeaf, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, weight, rowTrain, levelsTrain, validate);

  unsigned int ctgWidth = levelsTrain.length();
  bool test = !Rf_isNull(sYTest);
  IntegerVector yTest = test ? IntegerVector(sYTest) - 1 : IntegerVector(0);
  CharacterVector levelsTest = test ? as<CharacterVector>(IntegerVector(sYTest).attr("levels")) : CharacterVector(0);
  IntegerVector levelMatch = test ? match(levelsTest, levelsTrain) : IntegerVector(0);
  unsigned int testWidth;
  unsigned int testLength = yTest.length();
  bool dimFixup = false;
  if (test) {
    if (is_true(any(levelsTest != levelsTrain))) {
      dimFixup = true;
      IntegerVector sq = seq(0, levelsTest.length() - 1);
      IntegerVector idxNonMatch = sq[is_na(levelMatch)];
      if (idxNonMatch.length() > 0) {
	warning("Unreachable test levels not encountered in training");
	int proxy = ctgWidth + 1;
	for (R_len_t i = 0; i < idxNonMatch.length(); i++) {
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

  std::vector<unsigned int> confCore(testWidth * ctgWidth);
  std::vector<double> misPredCore(testWidth);
  std::vector<unsigned int> censusCore(nRow * ctgWidth);
  std::vector<unsigned int> yPred(nRow);
  NumericVector probCore = doProb ? NumericVector(nRow * ctgWidth) : NumericVector(0);
  Predict::Classification(valNum, rowStart, runLength, predStart, (valNum.size() == 0 && nPredNum > 0) ? transpose(blockNum).begin() : 0, nPredFac > 0 ? (unsigned int*) transpose(blockFac).begin() : 0, nPredNum, nPredFac, forestNode, origin, nTree, facSplit, facLen, facOrig, nFac, leafOrigin, leafNode, leafCount, bagBits, rowTrain, weight, ctgWidth, yPred, &censusCore[0], testCore, test ? &confCore[0] : 0, misPredCore, doProb ? probCore.begin() : 0);

  List predBlock(sPredBlock);
  IntegerMatrix census = transpose(IntegerMatrix(ctgWidth, nRow, &censusCore[0]));
  census.attr("dimnames") = List::create(predBlock["rowNames"], levelsTrain);
  NumericMatrix prob = doProb ? transpose(NumericMatrix(ctgWidth, nRow, probCore.begin())) : NumericMatrix(0);
  if (doProb) {
    prob.attr("dimnames") = List::create(predBlock["rowNames"], levelsTrain);
  }

  // OOB error is mean(prediction != training class)
  unsigned int missed = 0;
  for (unsigned int i = 0; i < nRow; i++) { // Bases to unity for front end.
    if (test) {
      missed += (unsigned int) yTest[i] != yPred[i];
    }
    yPred[i] = yPred[i] + 1;
  }
  double oobError = double(missed) / nRow;


  List prediction;
  if (test) {
    IntegerMatrix conf = transpose(IntegerMatrix(ctgWidth, testWidth, &confCore[0]));
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
      for (R_len_t i = 0; i < levelsTest.length(); i++) {
	misPred[i] = misPredCore[i];
      }
    }

    misPred.attr("names") = levelsTest;
    conf.attr("dimnames") = List::create(levelsTest, levelsTrain);
    prediction = List::create(
      _["misprediction"] = misPred,
      _["oobError"] = oobError,
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

  RcppLeaf::Clear();
  RcppForest::Clear();
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
RcppExport SEXP RcppPredictQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sQuantVec, SEXP sQBin, SEXP sYTest, bool validate) {
  unsigned int nPredNum, nPredFac, nRow;
  NumericMatrix blockNum;
  IntegerMatrix blockFac;
  std::vector<double> valNum;
  std::vector<unsigned int> rowStart;
  std::vector<unsigned int> runLength;
  std::vector<unsigned int> predStart;
  RcppPredblock::Unwrap(sPredBlock, nRow, nPredNum, nPredFac, blockNum, blockFac, valNum, rowStart, runLength, predStart);
    
  unsigned int *origin, *facOrig, *facSplit;
  ForestNode *forestNode;
  unsigned int nTree, nFac, nodeEnd;
  size_t facLen;
  RcppForest::Unwrap(sForest, origin, nTree, facSplit, facLen, facOrig, nFac, forestNode, nodeEnd);

  std::vector<double> yTrain;
  std::vector<unsigned int> leafOrigin;
  LeafNode *leafNode;
  unsigned int leafCount;
  BagLeaf *bagLeaf;
  unsigned int bagLeafTot;
  unsigned int *bagBits;

  // Quantile prediction requires full bagging information regardless
  // whether validating.
  RcppLeaf::UnwrapReg(sLeaf, yTrain, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, true);

  std::vector<double> yPred(nRow);
  std::vector<double> quantVecCore(as<std::vector<double> >(sQuantVec));
  std::vector<double> qPredCore(nRow * quantVecCore.size());
  Predict::Quantiles(valNum, rowStart, runLength, predStart, (valNum.size() == 0 && nPredNum > 0) ? transpose(blockNum).begin() : 0, nPredFac > 0 ? (unsigned int*) transpose(blockFac).begin() : 0, nPredNum, nPredFac, forestNode, origin, nTree, facSplit, facLen, facOrig, nFac, leafOrigin, leafNode, leafCount, bagLeaf, bagLeafTot, bagBits, yTrain, yPred, quantVecCore, as<unsigned int>(sQBin), qPredCore, validate);
  
  NumericMatrix qPred(transpose(NumericMatrix(quantVecCore.size(), nRow, qPredCore.begin())));
  List prediction;
  if (!Rf_isNull(sYTest)) {
    double rsq, mae;
    double mse = MSE(&yPred[0], NumericVector(sYTest), rsq, mae);
    prediction = List::create(
 	 _["yPred"] = yPred,
	 _["qPred"] = qPred,
	 _["mse"] = mse,
	 _["mae"] = mae,
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

  RcppLeaf::Clear();
  RcppForest::Clear();
  return prediction;
}


RcppExport SEXP RcppValidateQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sYTest, SEXP sQuantVec, SEXP sQBin) {
  return RcppPredictQuant(sPredBlock, sForest, sLeaf, sQuantVec, sQBin, sYTest, true);
}


RcppExport SEXP RcppTestQuant(SEXP sPredBlock, SEXP sForest, SEXP sLeaf, SEXP sQuantVec, SEXP sQBin, SEXP sYTest) {
  return RcppPredictQuant(sPredBlock, sForest, sLeaf, sQuantVec, sQBin, sYTest, false);
}
// Copyright (C)  2012-2017   Mark Seligman
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


#include "rcppRowrank.h"
#include "rowrank.h"

// Testing only:
//#include <iostream>
//using namespace std;


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

  std::vector<unsigned int> rank;
  std::vector<unsigned int> row;
  std::vector<unsigned int> runLength;
  std::vector<unsigned int> numOff(nPredNum);
  std::vector<double> numVal;
  if (nPredNum > 0) {
    if (!Rf_isNull(predBlock["blockNumRLE"])) {
      List blockNumRLE((SEXP) predBlock["blockNumRLE"]);
      if (!blockNumRLE.inherits("BlockNumRLE"))
	stop("Expecting BlockNumRLE");
      RowRank::PreSortNumRLE(NumericVector((SEXP) blockNumRLE["valNum"]).begin(), (unsigned int *) IntegerVector((SEXP) blockNumRLE["rowStart"]).begin(), (unsigned int*) IntegerVector((SEXP) blockNumRLE["runLength"]).begin(), nPredNum, nRow, row, rank, runLength, numOff, numVal);
    }
    else {
      NumericMatrix blockNum = predBlock["blockNum"];
      RowRank::PreSortNum(blockNum.begin(), nPredNum, nRow, row, rank, runLength, numOff, numVal);
    }
  }

  if (nPredFac > 0) {
    IntegerMatrix blockFac = predBlock["blockFac"];
    RowRank::PreSortFac((unsigned int*) blockFac.begin(), nPredFac, nRow, row, rank, runLength);
  }

  List rowRank = List::create(
      _["row"] = row,			      
      _["rank"] = rank,
      _["runLength"] = runLength,
      _["numOff"] = numOff,
      _["numVal"] = numVal
    );
  rowRank.attr("class") = "RowRank";

  return rowRank;
}


IntegerVector RcppRowrank::iv1 = IntegerVector(0);
IntegerVector RcppRowrank::iv2 = IntegerVector(0);
IntegerVector RcppRowrank::iv3 = IntegerVector(0);
IntegerVector RcppRowrank::iv4 = IntegerVector(0);
NumericVector RcppRowrank::nv1 = NumericVector(0);

void RcppRowrank::Unwrap(SEXP sRowRank, unsigned int *&feNumOff, double *&feNumVal, unsigned int *&feRow, unsigned int *&feRank, unsigned int *&feRLE, unsigned int &rleLength) {
  List rowRank(sRowRank);
  if (!rowRank.inherits("RowRank"))
    stop("Expecting RowRank");

  iv1 = IntegerVector((SEXP) rowRank["numOff"]);
  feNumOff = (unsigned int*) &iv1[0];

  nv1 = NumericVector((SEXP) rowRank["numVal"]);
  feNumVal = (double *) &nv1[0];
  
  iv2 = IntegerVector((SEXP) rowRank["row"]);
  feRow = (unsigned int *) &iv2[0];

  iv3 = IntegerVector((SEXP) rowRank["rank"]);
  feRank = (unsigned int *) &iv3[0];

  iv4 = IntegerVector((SEXP) rowRank["runLength"]);
  feRLE = (unsigned int *) &iv4[0];
  rleLength = iv4.length();
}


void RcppRowrank::Clear() {
  iv1 = IntegerVector(0);
  iv2 = IntegerVector(0);
  iv3 = IntegerVector(0);
  iv4 = IntegerVector(0);
  nv1 = NumericVector(0);
}
// Copyright (C)  2012-2017   Mark Seligman
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

////#include <iostream>
//using namespace std;


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
    for (R_len_t i = 0; i < classWeight.length(); i++) {
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
RcppExport SEXP RcppTrainCtg(SEXP sPredBlock, SEXP sRowRank, SEXP sYOneBased, SEXP sNTree, SEXP sNSamp, SEXP sSampleWeight, SEXP sWithRepl, SEXP sTrainBlock, SEXP sMinNode, SEXP sMinRatio, SEXP sTotLevels, SEXP sPredFixed, SEXP sSplitQuant, SEXP sProbVec, SEXP sThinLeaves, SEXP sClassWeight) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");

  IntegerVector yOneBased(sYOneBased);
  CharacterVector levels(yOneBased.attr("levels"));
  unsigned int ctgWidth = levels.length();

  IntegerVector y = yOneBased - 1;
  std::vector<double> proxy(y.length());
  NumericVector classWeight(as<NumericVector>(sClassWeight));
  RcppProxyCtg(y, classWeight, proxy);

  unsigned int nTree = as<unsigned int>(sNTree);
  std::vector<double> sampleWeight(as<std::vector<double> >(sSampleWeight));

  unsigned int nPredNum = as<unsigned int>(predBlock["nPredNum"]);
  unsigned int nPredFac = as<unsigned int>(predBlock["nPredFac"]);
  unsigned int nPred = nPredNum + nPredFac;

  List signature(as<List>(predBlock["signature"]));
  IntegerVector predMap(as<IntegerVector>(signature["predMap"]));

  NumericVector predProb = NumericVector(sProbVec)[predMap];
  NumericVector splitQuant = NumericVector(sSplitQuant)[predMap];

  Train::Init(nPred, nTree, as<unsigned int>(sNSamp), sampleWeight, as<bool>(sWithRepl), as<unsigned int>(sTrainBlock), as<unsigned int>(sMinNode), as<double>(sMinRatio), as<unsigned int>(sTotLevels), ctgWidth, as<unsigned int>(sPredFixed), splitQuant.begin(), predProb.begin(), as<bool>(sThinLeaves));

  std::vector<unsigned int> facCard(as<std::vector<unsigned int> >(predBlock["facCard"]));
  std::vector<unsigned int> origin(nTree);
  std::vector<unsigned int> facOrig(nTree);
  std::vector<unsigned int> leafOrigin(nTree);
  std::vector<double> predInfo(nPred);

  std::vector<ForestNode> forestNode;
  std::vector<unsigned int> facSplit;
  std::vector<LeafNode> leafNode;
  std::vector<BagLeaf> bagLeaf;
  std::vector<unsigned int> bagBits;
  std::vector<double> weight;

  double *feNumVal;
  unsigned int *feNumOff, *feRow, *feRank, *feRLE, rleLength;
  RcppRowrank::Unwrap(sRowRank, feNumOff, feNumVal, feRow, feRank, feRLE, rleLength);

  Train::Classification(feRow, feRank, feNumOff, feNumVal, feRLE, rleLength, as<std::vector<unsigned int> >(y), ctgWidth, proxy, origin, facOrig, predInfo, facCard, forestNode, facSplit, leafOrigin, leafNode, bagLeaf, bagBits, weight);

  RcppRowrank::Clear();
  
  NumericVector infoOut(predInfo.begin(), predInfo.end());
  return List::create(
      _["forest"] = RcppForest::Wrap(origin, facOrig, facSplit, forestNode),
      _["leaf"] = RcppLeaf::WrapCtg(leafOrigin, leafNode, bagLeaf, bagBits, weight, yOneBased.length(), CharacterVector(yOneBased.attr("levels"))),
      _["predInfo"] = infoOut[predMap] // Maps back from core order.
  );
}


RcppExport SEXP RcppTrainReg(SEXP sPredBlock, SEXP sRowRank, SEXP sY, SEXP sNTree, SEXP sNSamp, SEXP sSampleWeight, SEXP sWithRepl, SEXP sTrainBlock, SEXP sMinNode, SEXP sMinRatio, SEXP sTotLevels, SEXP sPredFixed, SEXP sSplitQuant, SEXP sProbVec, SEXP sThinLeaves, SEXP sRegMono) {
  List predBlock(sPredBlock);
  if (!predBlock.inherits("PredBlock"))
    stop("Expecting PredBlock");

  List signature(as<List>(predBlock["signature"]));
  IntegerVector predMap(as<IntegerVector>(signature["predMap"]));
  
  unsigned int nTree = as<unsigned int>(sNTree);
  std::vector<double> sampleWeight(as<std::vector<double> >(sSampleWeight));

  unsigned int nPredNum = as<unsigned int>(predBlock["nPredNum"]);
  unsigned int nPredFac = as<unsigned int>(predBlock["nPredFac"]);
  unsigned int nPred = nPredNum + nPredFac;

  NumericVector predProb = NumericVector(sProbVec)[predMap];
  NumericVector regMono = NumericVector(sRegMono)[predMap];
  NumericVector splitQuant = NumericVector(sSplitQuant)[predMap];
  
  Train::Init(nPred, nTree, as<unsigned int>(sNSamp), sampleWeight, as<bool>(sWithRepl), as<unsigned int>(sTrainBlock), as<unsigned int>(sMinNode), as<double>(sMinRatio), as<unsigned int>(sTotLevels), 0, as<unsigned int>(sPredFixed), splitQuant.begin(), predProb.begin(), as<bool>(sThinLeaves), regMono.begin());

  double *feNumVal;
  unsigned int *feRow, *feNumOff, *feRank, *feRLE, rleLength;
  RcppRowrank::Unwrap(sRowRank, feNumOff, feNumVal, feRow, feRank, feRLE, rleLength);

  NumericVector y(sY);
  NumericVector yOrdered = clone(y).sort();
  IntegerVector row2Rank = match(y, yOrdered) - 1;

  std::vector<unsigned int> origin(nTree);
  std::vector<unsigned int> facOrig(nTree);
  std::vector<unsigned int> leafOrigin(nTree);
  std::vector<double> predInfo(nPred);

  std::vector<ForestNode> forestNode;
  std::vector<LeafNode> leafNode;
  std::vector<BagLeaf> bagLeaf;
  std::vector<unsigned int> bagBits;
  std::vector<unsigned int> facSplit;

  const std::vector<unsigned int> facCard(as<std::vector<unsigned int> >(predBlock["facCard"]));
  Train::Regression(feRow, feRank, feNumOff, feNumVal, feRLE, rleLength, as<std::vector<double> >(y), as<std::vector<unsigned int> >(row2Rank), origin, facOrig, predInfo, facCard, forestNode, facSplit, leafOrigin, leafNode, bagLeaf, bagBits);

  RcppRowrank::Clear();

  // Temporary copy for subscripted access by IntegerVector.
  NumericVector infoOut(predInfo.begin(), predInfo.end()); 
  return List::create(
      _["forest"] = RcppForest::Wrap(origin, facOrig, facSplit, forestNode),
      _["leaf"] = RcppLeaf::WrapReg(leafOrigin, leafNode, bagLeaf, bagBits, as<std::vector<double> >(y)),
      _["predInfo"] = infoOut[predMap] // Maps back from core order.
    );
}
