// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// @file HoughTransformSeeder.hpp
// @author Max Goblirsch, based on work by Riley Xu and Jahred Adelman
// @brief Implements a set of tools to
/// implement a hough transform.

#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <optional>
#include <set>
#include <map>
#include <unordered_set>

#include "HoughVectors.hpp"

#pragma once

namespace Acts {
namespace HoughTransformUtils {

/// this type is responsible for encoding the parameters of our hough space
using coordType = double;

// this type is used to encode hit counts.
// Floating point to allow hit weights to be applied
using yieldType = float;

/// @brief this function represents a mapping of a coordinate point in detector space to a line in
/// hough space. Given the value of the first hough coordinate, it shall return
/// the corresponding second coordinate according to the line parametrisation.
/// Should be implemented by the user.
/// @tparam PointType: The class representing a point in detector space (can differ between implementations)
template <class PointType>
using lineParametrisation =
    std::function<coordType(coordType, const PointType&)>;

/// @brief struct to define the ranges of the hough histogram.
/// Used to move between parameter and bin index coordinates.
/// Disconnected from the hough plane binning to be able to re-use
/// a plane with a given binning for several parameter ranges
struct houghAxisRanges {
  coordType xMin = 0.0f;  // minimum value of the first coordinate
  coordType xMax = 0.0f;  // maximim value of the first coordinate
  coordType yMin = 0.0f;  // minimum value of the second coordinate
  coordType yMax = 0.0f;  // maximum value of the second coordinate
};

/// convenience functions to link bin indices to axis coordinates

/// @brief find the bin index corresponding to a certain abscissa
/// of the coordinate axis, based on the axis limits and binning.
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param val: value to find the corresponding bin for
/// @return the bin number.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
int binIndex(double min, double max, unsigned nSteps, double val) {
  return static_cast<int>((val - min) / (max - min) * nSteps);
}
// Returns the lower bound of the bin specified by step
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param binIndex: The index of the bin
/// @return the parameter value at the lower bin edge.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
double lowerBinEdge(double min, double max, unsigned nSteps, size_t binIndex) {
  return min + (max - min) * binIndex / nSteps;
}
// Returns the lower bound of the bin specified by step
/// @param min: Start of axis range
/// @param max: End of axis range
/// @param nSteps: Number of bins in axis
/// @param binIndex: The index of the bin
/// @return the parameter value at the bin center.
/// No special logic to prevent over-/underflow, checking these is
/// left to the caller
double binCenter(double min, double max, unsigned nSteps, size_t binIndex) {
  return min + (max - min) * 0.5 * (2 * binIndex + 1) / nSteps;
}

/// @brief data class for the information to store for each
/// cell of the hough histogram.
/// @tparam identifier_t: Type of the identifier to associate to the hits
///                       Should be sortable. Used to uniquely identify each
///                       hit and to eventually return the list of hits per cell
template <class identifier_t>
class HoughCell {
 public:
  /// @brief construct the cell as empty
  HoughCell() = default;
  /// @brief add an entry to this cell
  /// @param identifier: Identifier of the hit (used to distinguish hits from another)
  /// @param layer: Layer of the hit (used when counting layers)
  /// @param weight: Optional weight to assign to the hit
  void fill(const identifier_t& identifier, unsigned layer,
            yieldType weight = 1.);
  /// @brief access the number of layers with hits compatible with this cell
  yieldType nLayers() const { return m_nLayers; }
  /// @brief access the number of unique hits compatible with this cell
  yieldType nHits() const { return m_nHits; }
  /// @brief access the set of layers compatible with this cell
  const std::unordered_set<unsigned>& layers() const { return m_layers; }
  /// @brief access the set of unique hits compatible with this cell
  const std::unordered_set<identifier_t>& hits() const { return m_hits; }
  /// @brief reset this cell, removing any existing content.
  void reset();

 private:
  /// data members

  yieldType m_nLayers =
      0;                  // (weighted) number of layers with hits on this cell
  yieldType m_nHits = 0;  // (weighted) number of unique hits on this cell
  std::unordered_set<unsigned> m_layers =
      {};  // set of layers with hits on this cell
  std::unordered_set<identifier_t> m_hits =
      {};  // set of unique hits on this cell
};

/// @brief Configuration - number of bins in each axis.
/// The Hough plane is agnostic of how the bins map to
/// coordinates, allowing to re-use a plane for several
/// (sub) detectors of different dimensions if the bin number
/// remains applicable
struct HoughPlaneConfig {
  size_t nBinsX = 0;  // number of bins in the first coordinate
  size_t nBinsY = 0;  // number of bins in the second coordinate
};

/// @brief Representation of the hough plane - the histogram used
/// for the hough transform with methods to fill and evaluate
/// the histogram. Templated to a class used as identifier for the hits
template <class identifier_t>
class HoughPlane {
 public:
  /// @brief hough histogram representation as a 2D-indexable vector of hough cells
  using HoughHist = vector2D<HoughCell<identifier_t>>;

  /// @brief instantiate the (empty) hough plane
  /// @param cfg: configuration
  HoughPlane(const HoughPlaneConfig& cfg);

  /// fill and reset methods to modify the grid content

  /// @brief add one measurement to the hough plane
  /// @tparam PointType: Type of the objects to use when adding measurements (e.g. experiment EDM object)
  /// @param measurement: The measurement to add
  /// @param axisRanges: Ranges of the hough axes, used to map the bin numbers to parameter values
  /// @param linePar: The function y(x) parametrising the hough space line for a given measurement
  /// @param widthPar: The function dy(x) parametrising the width of the y(x) curve
  ///                   for a given measurement
  /// @param identifier: The unique identifier for the given hit
  /// @param layer: A layer index for this hit
  /// @param weight: An optional weight to assign to this hit
  template <class PointType>
  void fill(const PointType& measurement, const houghAxisRanges& axisRanges,
            lineParametrisation<PointType> linePar,
            lineParametrisation<PointType> widthPar,
            const identifier_t& identifier, unsigned layer = 0,
            yieldType weight = 1.0f);
  /// @brief resets the contents of the grid. Can be used to avoid reallocating the histogram
  /// when switching regions / (sub)detectors
  void reset();

  //// user-facing accessors

  /// @brief get the layers with hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param y: bin index in the second coordinate
  /// @return the set of layer indices that have hits for this cell
  const std::unordered_set<unsigned>& layers(size_t xBin, size_t yBin) const {
    return m_houghHist(xBin, yBin).layers();
  }

  /// @brief get the (weighted) number of layers  with hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the (weighed) number of layers that have hits for this cell
  yieldType nLayers(size_t xBin, size_t yBin) const {
    return m_houghHist(xBin, yBin).nLayers();
  }

  /// @brief get the identifiers of all hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the set of identifiers of the hits for this cell
  const std::unordered_set<identifier_t>& hitIds(size_t xBin,
                                                 size_t yBin) const {
    return m_houghHist(xBin, yBin).hits();
  }
  /// @brief get the (weighted) number of hits in one cell of the histogram
  /// @param xBin: bin index in the first coordinate
  /// @param yBin: bin index in the second coordinate
  /// @return the (weighted) number of hits for this cell
  yieldType nHits(size_t xBin, size_t yBin) const {
    return m_houghHist(xBin, yBin).nHits();
  }

  /// @brief get the number of bins on the first coordinate
  size_t nBinsX() const { return m_cfg.nBinsX; }
  /// @brief get the number of bins on the second coordinate
  size_t nBinsY() const { return m_cfg.nBinsY; }

  /// @brief get the maximum number of (weighted) hits seen in a single
  /// cell across the entire histrogram.
  yieldType maxHits() const { return m_maxHits; }

  /// @brief get the list of cells with non-zero content.
  /// Useful for peak-finders in sparse data
  /// to avoid looping over all cells
  const std::set<std::pair<size_t, size_t>>& getNonEmptyBins() const {
    return m_touchedBins;
  }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) hits across the entire histogram
  std::pair<size_t, size_t> locMaxHits() const { return m_maxLocHits; }

  /// @brief get the maximum number of (weighted) layers with hits  seen
  /// in a single cell across the entire histrogram.
  yieldType maxLayers() const { return m_maxLayers; }

  /// @brief get the bin indices of the cell containing the largest number
  /// of (weighted) layers with hits across the entire histogram
  std::pair<size_t, size_t> locMaxLayers() const { return m_maxLocLayers; }

 private:
  /// @brief Helper method to fill a bin of the hough histogram.
  /// Updates the internal helper data structures (maximum tracker etc).
  /// @param binX: bin number along x
  /// @param binY: bin number along y
  /// @param identifier: hit identifier
  /// @param layer: layer index
  /// @param w: optional hit weight
  void fillBin(size_t binX, size_t binY, const identifier_t& identifier,
               unsigned layer, double w = 1.0f);

  yieldType m_maxHits = 0.0f;    // track the maximum number of hits seen
  yieldType m_maxLayers = 0.0f;  // track the maximum number of layers seen
  std::pair<size_t, size_t> m_maxLocHits = {
      0, 0};  // track the location of the maximum in hits
  std::pair<size_t, size_t> m_maxLocLayers = {
      0, 0};  // track the location of the maximum in layers
  std::set<std::pair<size_t, size_t>> m_touchedBins =
      {};                  // track the bins with non-trivial content
  HoughPlaneConfig m_cfg;  // the configuration object
  HoughHist m_houghHist;   // the histogram data object
};

/// example peak finders.
namespace PeakFinders {
/// configuration for the LayerGuidedCombinatoric peak finder
struct LayerGuidedCombinatoricConfig {
  yieldType threshold = 3.0f;  // min number of layers to obtain a maximum
  int localMaxWindowSize = 0;  // Only create candidates from a local maximum
};

/// @brief Peak finder inspired by ATLAS ITk event filter developments.
/// Builds peaks based on layer counts and allows for subsequent resolution
/// of the combinatorics by building multiple candidates per peak if needed.
/// In flat regions, peak positions are moved towards smaller values of the
/// second and first coordinate.
/// @tparam identifier_t: The identifier type to use. Should match the one used in the hough plane.
template <class identifier_t>
class LayerGuidedCombinatoric {
 public:
  /// @brief data class representing the found maxima.
  /// Here, we just have a list of cluster identifiers
  struct Maximum {
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  /// @brief constructor
  /// @param cfg: Configuration object
  LayerGuidedCombinatoric(const LayerGuidedCombinatoricConfig& cfg);

  /// @brief main peak finder method.
  /// @param plane: Filled hough plane to search
  /// @return vector of found maxima
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane) const;

 private:
  /// @brief check if a given bin is a local maximum.
  /// @param plane: The filled hough plane
  /// @param xBin: x bin index
  /// @param yBin: y bin index
  /// @return true if a maximum, false otherwise
  bool passThreshold(const HoughPlane<identifier_t>& plane, size_t xBin,
                     size_t yBin) const;  // did we pass extensions?

  LayerGuidedCombinatoricConfig m_cfg;  // configuration data object
};
/// @brief Configuration for the IslandsAroundMax peak finder
struct IslandsAroundMaxConfig {
  yieldType threshold =
      3.0f;  // min number of weigted hits required in a maximum
  yieldType fractionCutoff =
      0;  // Fraction of the global maximum at which to cut off maxima
  std::pair<coordType, coordType> minSpacingBetweenPeaks = {
      0.0f, 0.0f};  // minimum distance of a new peak from existing peaks in
                    // parameter space
};
/// @brief Peak finder inspired by ATLAS muon reconstruction.
/// Looks for regions above a given fraction of the global maximum
/// hit count and connects them into islands comprising adjacent bins
/// above the threshold. Peak positions are averaged across cells in the island,
/// weighted by hit counts
/// @tparam identifier_t: The identifier type to use. Should match the one used in the hough plane.
template <class identifier_t>
class IslandsAroundMax {
 public:
  /// @brief data struct representing a local maximum.
  /// Comes with a position estimate and a list of hits within the island
  struct Maximum {
    coordType x = 0;   // x value of the maximum
    coordType y = 0;   // y value of the maximum
    coordType wx = 0;  // x width of the maximum
    coordType wy = 0;  // y width of the maximum
    std::unordered_set<identifier_t> hitIdentifiers =
        {};  // identifiers of contributing hits
  };
  /// @brief constructor.
  /// @param cfg: configuration object
  IslandsAroundMax(const IslandsAroundMaxConfig& cfg);

  /// @brief main peak finder method.
  /// @param plane: The filled hough plane to search
  /// @ranges: The axis ranges used for mapping between parameter space and bins.
  /// @return List of the found maxima
  std::vector<Maximum> findPeaks(const HoughPlane<identifier_t>& plane,
                                 const houghAxisRanges& ranges);

 private:
  /// @brief method to incrementally grow an island by adding adjacent cells
  /// Performs a breadth-first search for neighbours above threshold and adds
  /// them to candidate. Stops when no suitable neighbours are left.
  /// @param inMaximum: List of cells found in the island. Incrementally populated by calls to the method
  /// @param toExplore: List of neighbour cell candidates left to explore. Method will not do anything once this is empty
  /// @param threshold: the threshold to apply to check if a cell should be added to an island
  /// @param yieldMap: A map of the hit content of above-threshold cells. Used cells will be set to empty content to avoid re-use by subsequent calls
  void extendMaximum(std::vector<std::pair<size_t, size_t>>& inMaximum,
                     std::vector<std::pair<size_t, size_t>>& toExplore,
                     yieldType threshold, std::map<std::pair<size_t,size_t>, yieldType>& yieldMap);

  IslandsAroundMaxConfig m_cfg;  // configuration data object

  /// @brief array of steps to consider when exploring neighbouring cells.
  const std::array<std::pair<int, int>, 8> m_stepDirections{
      std::make_pair(-1, -1), std::make_pair(0, -1), std::make_pair(1, -1),
      std::make_pair(-1, 0),  std::make_pair(1, 0),  std::make_pair(-1, 1),
      std::make_pair(0, 1),   std::make_pair(1, 1)};
};
}  // namespace PeakFinders
}  // namespace HoughTransformUtils
}  // namespace Acts

#include "HoughTransformUtils.ipp"
