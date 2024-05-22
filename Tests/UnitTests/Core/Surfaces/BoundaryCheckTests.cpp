// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <vector>

#include "BoundaryCheckTestsRefs.hpp"

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Surfaces)

// See: https://en.wikipedia.org/wiki/Bounding_volume
//
// Aligned box w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxSimple) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  AlignedBoxBoundaryCheck check(ll, ur, BoundaryTolerance::None());
  BOOST_CHECK(check.inside({0, 0}, std::nullopt));
  BOOST_CHECK(!check.inside({2, 2}, std::nullopt));
  BOOST_CHECK(!check.inside({0, 2}, std::nullopt));
  BOOST_CHECK(!check.inside({2, 0}, std::nullopt));
}
// Aligned box w/ tolerance check along first axis
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxToleranceLoc0) {
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  AlignedBoxBoundaryCheck check(ll, ur,
                                BoundaryTolerance::AbsoluteBound(1.5, 100.0));
  BOOST_CHECK(check.inside({0, 0}, std::nullopt));
  BOOST_CHECK(check.inside({2, 2}, std::nullopt));
  BOOST_CHECK(!check.inside({4, 4}, std::nullopt));
  BOOST_CHECK(check.inside({0, 2}, std::nullopt));
  BOOST_CHECK(check.inside({2, 0}, std::nullopt));
}

// Aligned box w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckBoxCovariance) {
  SquareMatrix2 cov;
  cov << 1, 0.5, 0.5, 2;
  Vector2 ll(-1, -1);
  Vector2 ur(1, 1);
  AlignedBoxBoundaryCheck check(
      ll, ur, BoundaryTolerance::Chi2Bound(cov.inverse(), 3.0));
  BOOST_CHECK(check.inside({0, 0}, std::nullopt));
  BOOST_CHECK(check.inside({2, 2}, std::nullopt));
  BOOST_CHECK(!check.inside({4, 4}, std::nullopt));
  BOOST_CHECK(check.inside({0, 3}, std::nullopt));
  BOOST_CHECK(check.inside({3, 0}, std::nullopt));
}

// Triangle w/ simple check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleSimple) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  PolygonBoundaryCheck check(vertices, BoundaryTolerance::None());
  BOOST_CHECK(check.inside({0, 0}, std::nullopt));
  BOOST_CHECK(check.inside({0, 1}, std::nullopt));
  BOOST_CHECK(!check.inside({2, 2}, std::nullopt));
  BOOST_CHECK(!check.inside({0, -1}, std::nullopt));
}
// Triangle w/ covariance check
BOOST_AUTO_TEST_CASE(BoundaryCheckTriangleCovariance) {
  Vector2 vertices[] = {{-2, 0}, {2, 0}, {0, 2}};
  SquareMatrix2 cov;
  cov << 0.5, 0, 0, 0.5;
  PolygonBoundaryCheck check(vertices,
                             BoundaryTolerance::Chi2Bound(cov.inverse(), 4.1));
  BOOST_CHECK(check.inside({0, 0}, std::nullopt));
  BOOST_CHECK(check.inside({0, 1}, std::nullopt));
  BOOST_CHECK(check.inside({0, 2}, std::nullopt));
  BOOST_CHECK(check.inside({0, 3}, std::nullopt));
  BOOST_CHECK(check.inside({0, 4}, std::nullopt));
  BOOST_CHECK(!check.inside({0, 5}, std::nullopt));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
