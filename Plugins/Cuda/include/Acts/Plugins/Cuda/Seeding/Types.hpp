// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {
namespace Cuda {
namespace details {

/// Structure used in the CUDA-based triplet finding
struct Triplet {
  int bIndex;
  int tIndex;
  float topRadius;
  float impactParameter;
  float invHelixDiameter;
  float weight;
};  // struct Triplet

}  // namespace details
}  // namespace Cuda
}  // namespace Acts
