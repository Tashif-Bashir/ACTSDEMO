// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <iomanip>
#include <ostream>

std::ostream& Acts::operator<<(std::ostream& os, Acts::GeometryIdentifier id) {
  // zero represents an invalid/undefined identifier
  if (not id.value()) {
    os << "(undefined)";
    return os;
  }

  static const char* const names[] = {
      "volume=", "boundary=", "layer=", "approach=", "sensitive=",
  };
  const GeometryIdentifier::Value levels[] = {
      id.volume(), id.boundary(), id.layer(), id.approach(), id.sensitive(),
  };

  os << '(';
  bool writeSeparator = false;
  for (auto i = 0u; i < 5u; ++i) {
    if (levels[i]) {
      if (writeSeparator) {
        os << '|';
      }
      os << names[i] << levels[i];
      writeSeparator = true;
    }
  }
  os << ')';
  return os;
}
