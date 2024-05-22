// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/detail/periodic.hpp"

namespace Acts {

inline std::vector<double> CylinderBounds::values() const {
  std::vector<double> valvector;
  valvector.insert(valvector.begin(), m_values.begin(), m_values.end());
  return valvector;
}

inline bool CylinderBounds::coversFullAzimuth() const {
  return m_closed;
}

inline void CylinderBounds::checkConsistency() noexcept(false) {
  if (get(eR) <= 0.) {
    throw std::invalid_argument("CylinderBounds: invalid radial setup.");
  }
  if (get(eHalfLengthZ) <= 0.) {
    throw std::invalid_argument("CylinderBounds: invalid length setup.");
  }
  if (get(eHalfPhiSector) <= 0. || get(eHalfPhiSector) > M_PI) {
    throw std::invalid_argument("CylinderBounds: invalid phi sector setup.");
  }
  if (get(eAveragePhi) != detail::radian_sym(get(eAveragePhi))) {
    throw std::invalid_argument("CylinderBounds: invalid phi positioning.");
  }
  if (get(eBevelMinZ) != detail::radian_sym(get(eBevelMinZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at min Z.");
  }
  if (get(eBevelMaxZ) != detail::radian_sym(get(eBevelMaxZ))) {
    throw std::invalid_argument("CylinderBounds: invalid bevel at max Z.");
  }
}

}  // namespace Acts
