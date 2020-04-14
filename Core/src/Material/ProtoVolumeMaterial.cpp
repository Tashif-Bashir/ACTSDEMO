// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ProtoSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/ProtoVolumeMaterial.hpp"

Acts::ProtoVolumeMaterial::ProtoVolumeMaterial(const BinUtility& binUtility)
    : m_binUtility(binUtility), m_material() {}

std::ostream& Acts::ProtoVolumeMaterial::toStream(std::ostream& sl) const {
  sl << "Acts::ProtoVolumeMaterial : " << std::endl;
  if (m_binUtility.bins(0) * m_binUtility.bins(1) * m_binUtility.bins(2) > 1) {
    sl << "   - Number of Material bins [0,1] : " << m_binUtility.bins(0)
       << " / " << m_binUtility.bins(1) << " / " << m_binUtility.bins(2)
       << std::endl;
  } else {
    sl << "   - Homogeneous Material" << std::endl;
  }
  return sl;
}
