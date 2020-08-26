// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialProperties.hpp"

#include "Acts/Material/detail/AverageMaterials.hpp"

#include <climits>
#include <limits>
#include <numeric>
#include <ostream>

static constexpr auto eps = 2 * std::numeric_limits<float>::epsilon();

Acts::MaterialProperties::MaterialProperties(float thickness)
    : m_thickness(thickness) {}

Acts::MaterialProperties::MaterialProperties(const Material& material,
                                             float thickness)
    : m_material(material),
      m_thickness(thickness),
      m_thicknessInX0((eps < material.X0()) ? (thickness / material.X0()) : 0),
      m_thicknessInL0((eps < material.L0()) ? (thickness / material.L0()) : 0) {
}

Acts::MaterialProperties::MaterialProperties(
    const std::vector<MaterialProperties>& layers)
    : MaterialProperties(std::reduce(layers.begin(), layers.end(),
                                     MaterialProperties(),
                                     detail::averageMaterials)) {
  // NOTE 2020-08-26 msmk
  // the reduce work best (in the numerical stability sense) if the input
  // layers are sorted by by thickness/mass density. then, the late terms
  // of the averaging are only small corrections to the large average of
  // the initial layers. this could be enforce by sorting the layers first, but
  // i am not sure if this is actually a problem.
}

void Acts::MaterialProperties::scaleThickness(float scale) {
  m_thickness *= scale;
  m_thicknessInX0 *= scale;
  m_thicknessInL0 *= scale;
}

std::ostream& Acts::operator<<(std::ostream& os,
                               const MaterialProperties& materialProperties) {
  os << materialProperties.material()
     << "|t=" << materialProperties.thickness();
  return os;
}
