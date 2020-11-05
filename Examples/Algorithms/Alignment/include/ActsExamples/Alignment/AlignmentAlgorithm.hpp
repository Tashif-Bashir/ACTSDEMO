// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"

namespace ActsExamples {

class AlignmentAlgorithm final : public BareAlgorithm {
 public:
  using AlignResult = Acts::Result<ActsAlignment::AlignmentResult>;
  /// Alignment function that takes sets of input measurements, initial
  /// trackstate and alignment options and returns some alignment-specific
  /// result.
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<MeasurementCalibrator, Acts::VoidOutlierFinder>;
  using AlignmentFunction = std::function<AlignResult(
      const std::vector<std::vector<IndexSourceLink>>&,
      const TrackParametersContainer&,
      const ActsAlignment::AlignmentOptions<TrackFitterOptions>&)>;

  /// Create the fitter function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyways.
  static AlignmentFunction makeAlignmentFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      Options::BFieldVariant magneticField, Acts::Logging::Level lvl);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input source links collection.
    std::string inputSourceLinks;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for for each proto track.
    std::string inputInitialTrackParameters;
    /// Output fitted trajectories collection.
    std::string outputTrajectories;
    /// Type erased fitter function.
    AlignmentFunction align;
    /// The alignd transform updater
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater;
    /// The surfaces (or detector elements?) to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// The alignment mask at each iteration
    std::map<unsigned int, std::bitset<6>> iterationState;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff = 0.10;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<size_t, double> deltaChi2ONdfCutOff = {10, 0.00001};
    /// Maximum number of iterations
    size_t maxNumIterations = 100;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorihtm
  /// @param level is the logging level
  AlignmentAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the alignment algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algporithm flow
  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext& ctx) const final override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
