// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "ActsExamples/EventData/Index.hpp"

namespace ActsExamples {

struct AlreadyCalibratedCalibrator {
  using Proxy =
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::TrackStateProxy;
  using ConstProxy = Acts::MultiTrajectory<
      Acts::ConstVectorMultiTrajectory>::ConstTrackStateProxy;

  std::unordered_map<ActsExamples::Index, ConstProxy> callibratedStates;

  void calibrate(const Acts::GeometryContext& /*gctx*/, Proxy trackState) const;
};

}  // namespace ActsExamples
