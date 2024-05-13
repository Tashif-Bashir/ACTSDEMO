// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include <cstdint>

namespace Acts {
// constructor
template <typename SpacePointContainer>
HashingTrainingAlgorithm<SpacePointContainer>::HashingTrainingAlgorithm(
    const Acts::HashingTrainingConfig& cfg)
    : m_cfg(cfg) {
  if (m_cfg.f <= 0) {
    throw std::invalid_argument("Invalid f, f must be positive");
  }
  if (m_cfg.AnnoySeed <= 0) {
    throw std::invalid_argument(
        "Invalid Annoy random seed, Annoy random seed must be positive");
  }
}

template <typename SpacePointContainer>
AnnoyModel HashingTrainingAlgorithm<SpacePointContainer>::execute(
    SpacePointContainer spacePoints) const {
  const unsigned int AnnoySeed = m_cfg.AnnoySeed;
  const std::int32_t f = m_cfg.f;

  AnnoyModel annoyModel = AnnoyModel(f);

  using Scalar = Acts::ActsScalar;

  annoyModel.set_seed(AnnoySeed);

  unsigned int spacePointIndex = 0;
  // Add spacePoints parameters to Annoy
  for (const auto& spacePoint : spacePoints) {
    Scalar x = spacePoint->x() / Acts::UnitConstants::mm;
    Scalar y = spacePoint->y() / Acts::UnitConstants::mm;

    // Helix transform
    Scalar phi = atan2(y, x);

    double* vec = (double*)malloc(f * sizeof(double));
    vec[0] = (double)phi;
    if (f >= 2) {
      Scalar z = spacePoint->z() / Acts::UnitConstants::mm;
      Scalar r2 = x * x + y * y;
      Scalar rho = sqrt(r2 + z * z);
      Scalar theta = acos(z / rho);
      Scalar eta = -log(tan(0.5 * theta));
      vec[1] = (double)eta;
    }

    annoyModel.add_item(spacePointIndex, vec);
    spacePointIndex++;
  }

  unsigned int n_trees = 2 * f;

  annoyModel.build(n_trees);

  return annoyModel;
}

}  // namespace Acts
