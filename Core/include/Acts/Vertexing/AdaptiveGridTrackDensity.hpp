// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {

template <int trkGridSize = 15>
class AdaptiveGridTrackDensity {
  // Assert odd trkGridSize
  static_assert(trkGridSize % 2);

 public:
  /// @struct Config The configuration struct
  struct Config {
    /// @param zMinMax The minimum and maximum z-values (in mm) that
    /// should be covered by the main 1-dim density grid along the z-axis
    /// Note: This value together with 'mainGridSize' determines the
    /// overall bin size to be used as seen below
    Config(float binSize_) : binSize(binSize_) {}

    // Z size of one single bin in grid
    float binSize;  // mm

    // Do NOT use just the z-bin with the highest
    // track density, but instead check the (up to)
    // first three density maxima (only those that have
    // a maximum relative deviation of 'relativeDensityDev'
    // from the main maximum) and take the z-bin of the
    // maximum with the highest surrounding density sum
    bool useHighestSumZPosition = false;

    // The maximum relative density deviation from the main
    // maximum to consider the second and third maximum for
    // the highest-sum approach from above
    float maxRelativeDensityDev = 0.01;
  };

  AdaptiveGridTrackDensity(const Config& cfg) : m_cfg(cfg) {}


  /// @brief Returns the z position of maximum track density
  ///
  /// @param mainGrid The main 1-dim density grid along the z-axis
  ///
  /// @return The z position of maximum track density
  Result<float>
  getMaxZPosition(const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const;

  /// @brief Returns the z position of maximum track density and
  /// the estimated width
  ///
  /// @param mainGrid The main 1-dim density grid along the z-axis
  ///
  /// @return The z position of maximum track density and width
  Result<std::pair<float, float>>
    getMaxZPositionAndWidth(const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues) const;


  
  /// @brief Adds a single track to the overall grid density
  ///
  /// @param trk The track to be added
  /// @param mainGrid The main 1-dim density grid along the z-axis
  ///
  /// @return A pair storing information about the z-bin position
  /// the track was added (int) and the 1-dim density contribution
  /// of the track itself
  void addTrack(
      const BoundTrackParameters& trk,
      std::vector<float>& mainGridDensity,
      std::vector<int>& mainGridZValues) const;

  /// @brief Removes a track from the overall grid density
  ///
  /// @param zBin The center z-bin position the track needs to be
  /// removed from
  /// @param trkGrid The 1-dim density contribution of the track
  /// @param mainGrid The main 1-dim density grid along the z-axis
  // void removeTrackGridFromMainGrid(int zBin,
  //                                  const ActsVectorF<trkGridSize>& trkGrid,
  //                                  ActsVectorF<mainGridSize>& mainGrid) const;

 //private:
  /// @brief Helper function that acutally adds the track to the
  /// main density grid
  ///
  /// @param zBin The center z-bin position the track
  /// @param trkGrid The 1-dim density contribution of the track
  /// @param mainGrid The main 1-dim density grid along the z-axis
  // void addTrackGridToMainGrid(int zBin, const ActsVectorF<trkGridSize>& trkGrid,
  //                             ActsVectorF<mainGridSize>& mainGrid) const;

  /// @brief Helper function that modifies the main density grid
  /// (either adds or removes a track)
  ///
  /// @param zBin The center z-bin position the track
  /// @param trkGrid The 1-dim density contribution of the track
  /// @param mainGrid The main 1-dim density grid along the z-axis
  /// @param modifyModeSign Sign that determines the mode of modification,
  /// +1 for adding a track, -1 for removing a track
  // void modifyMainGridWithTrackGrid(int zBin,
  //                                  const ActsVectorF<trkGridSize>& trkGrid,
  //                                  ActsVectorF<mainGridSize>& mainGrid,
  //                                  int modifyModeSign) const;

  /// @brief Function that creates a 1-dim track grid (i.e. a vector)
  /// with the correct density contribution of a track along the z-axis
  ///
  /// @param offset Offset in d0 direction, to account for the 2-dim part
  /// of the Gaussian track distribution
  /// @param cov The track covariance matrix
  /// @param distCtrD The distance in d0 from the track position to its
  /// bin center in the 2-dim grid
  /// @param distCtrZ The distance in z0 from the track position to its
  /// bin center in the 2-dim grid
  ActsVectorF<trkGridSize> createTrackGrid(int offset, const SymMatrix2D& cov,
                                           float distCtrD,
                                           float distCtrZ) const;

  /// @brief Function that estimates the seed width based on the FWHM of
  /// the maximum density peak
  ///
  /// @param mainGrid The main 1-dim density grid along the z-axis
  /// @maxZ z-position of the maximum density value
  ///
  /// @return The width
  Result<float>
  estimateSeedWidth(
    const std::vector<float>& mainGridDensity,
    const std::vector<int>& mainGridZValues, float maxZ) const;

  /// @brief Helper to retrieve values according to a 2-dim normal distribution
   float normal2D(float d, float z, const SymMatrix2D& cov) const;

  /// @brief Checks the (up to) first three density maxima (only those that have
  /// a maximum relative deviation of 'relativeDensityDev' from the main
  /// maximum) and take the z-bin of the maximum with the highest surrounding
  /// density
  ///
  /// @param mainGrid The main 1-dim density grid along the z-axis
  ///
  /// @return The z-bin position
  // int getHighestSumZPosition(ActsVectorF<mainGridSize>& mainGrid) const;

  /// @brief Calculates the density sum of a z-bin and its two neighboring bins
  /// as needed for 'getHighestSumZPosition'
  ///
  /// @parem mainGrid The main 1-dim density grid along the z-axis
  /// @param pos The center z-bin positon
  ///
  /// @return The sum
  // double getDensitySum(const ActsVectorF<mainGridSize>& mainGrid,
  //                      int pos) const;

  Config m_cfg;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveGridTrackDensity.ipp"
