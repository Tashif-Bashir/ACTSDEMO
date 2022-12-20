// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Local include(s).
#include "TestHostCuts.hpp"

float TestHostCuts::seedWeight(
    const Acts::InternalSpacePoint<TestSpacePoint>& bottom,
    const Acts::InternalSpacePoint<TestSpacePoint>&,
    const Acts::InternalSpacePoint<TestSpacePoint>& top) const {
  float weight = 0;
  if (bottom.radius() > 150) {
    weight = 400;
  }
  if (top.radius() < 150) {
    weight = 200;
  }
  return weight;
}

bool TestHostCuts::singleSeedCut(
    float weight, const Acts::InternalSpacePoint<TestSpacePoint>& b,
    const Acts::InternalSpacePoint<TestSpacePoint>&,
    const Acts::InternalSpacePoint<TestSpacePoint>&) const {
  return !(b.radius() > 150. && weight < 380.);
}

std::vector<
  typename CandidatesForSpM<InternalSpacePoint<SpacePoint>>::output_type>
    ATLASCuts<SpacePoint>::cutPerMiddleSP(
        std::vector<
        typename CandidatesForSpM<InternalSpacePoint<SpacePoint>>::output_type>
        seedCandidates) const {
  std::vector<
      typename CandidatesForSpM<InternalSpacePoint<SpacePoint>>::output_type>
      newSeedsVector;
    if (seedCandidates.size() <= 1)
    return seedCandidates;

  newSeedsVector.push_back(std::move(seedCandidates[0]));
  std::size_t itLength = std::min(seedCandidates.size(), std::size_t(5));
  // don't cut first element
  for (std::size_t i(1); i < itLength; i++) {
    float weight = std::get<CandidatesForSpM<SpacePoint>::Components::WEIGHT>(
        seedCandidates[i]);
    const auto& bottom =
        std::get<CandidatesForSpM<SpacePoint>::Components::BSP>(
            seedCandidates[i]);
    if (weight > 200. or bottom->radius() > 43.) {
      newSeedsVector.push_back(std::move(seedCandidates[i]));
    }
  }
  return newSeedsVector;
}
}
