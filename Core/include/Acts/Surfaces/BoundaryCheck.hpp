// This file is part of the Acts project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>
#include <iterator>
#include <optional>
#include <variant>
#include <vector>

namespace Acts {

class BoundaryTolerance {
 public:
  struct Infinite {};

  struct None {};

  struct AbsoluteBound {
    double tolerance0{};
    double tolerance1{};

    AbsoluteBound() = default;
    AbsoluteBound(double tolerance0_, double tolerance1_)
        : tolerance0(tolerance0_), tolerance1(tolerance1_) {}
  };

  struct AbsoluteCartesian {
    double tolerance0{};
    double tolerance1{};

    AbsoluteCartesian() = default;
    AbsoluteCartesian(double tolerance0_, double tolerance1_)
        : tolerance0(tolerance0_), tolerance1(tolerance1_) {}
  };

  struct AbsoluteEuclidean {
    double tolerance{};

    AbsoluteEuclidean() = default;
    AbsoluteEuclidean(double tolerance_) : tolerance(tolerance_) {}
  };

  struct Chi2Bound {
    SquareMatrix2 weight = SquareMatrix2::Identity();
    double maxChi2{};

    Chi2Bound() = default;
    Chi2Bound(const SquareMatrix2& weight_, double maxChi2_)
        : weight(weight_), maxChi2(maxChi2_) {}
  };

  using Variant = std::variant<Infinite, None, AbsoluteBound, AbsoluteCartesian,
                               AbsoluteEuclidean, Chi2Bound>;

  BoundaryTolerance(const Infinite& infinite);
  BoundaryTolerance(const None& none);
  BoundaryTolerance(const AbsoluteBound& AbsoluteBound);
  BoundaryTolerance(const AbsoluteCartesian& absoluteCartesian);
  BoundaryTolerance(const AbsoluteEuclidean& absoluteEuclidean);
  BoundaryTolerance(const Chi2Bound& Chi2Bound);

  /// Construct from variant
  BoundaryTolerance(Variant variant);

  bool isInfinite() const;
  bool isNone() const;
  bool hasAbsoluteBound(bool isCartesian = false) const;
  bool hasAbsoluteCartesian() const;
  bool hasAbsoluteEuclidean() const;
  bool hasChi2Bound() const;

  bool hasTolerance() const;

  AbsoluteBound asAbsoluteBound(bool isCartesian = false) const;
  const AbsoluteCartesian& asAbsoluteCartesian() const;
  const AbsoluteEuclidean& asAbsoluteEuclidean() const;
  const Chi2Bound& asChi2Bound() const;

  std::optional<AbsoluteBound> asAbsoluteBoundOpt(
      bool isCartesian = false) const;

  bool isTolerated(const Vector2& distance,
                   const std::optional<SquareMatrix2>& jacobianOpt) const;

  bool hasMetric(bool hasJacobian) const;

  SquareMatrix2 getMetric(const std::optional<SquareMatrix2>& jacobian) const;

 private:
  Variant m_variant;

  /// Check if the boundary check is of a specific type.
  template <typename T>
  bool holdsVariant() const;

  /// Get the specific underlying type.
  template <typename T>
  const T& getVariant() const;

  template <typename T>
  const T* getVariantPtr() const;
};

class AlignedBoxBoundaryCheck {
 public:
  AlignedBoxBoundaryCheck(const Vector2& lowerLeft, const Vector2& upperRight,
                          BoundaryTolerance tolerance);

  /// Get the lower left corner of the box.
  const Vector2& lowerLeft() const;

  /// Get the upper right corner of the box.
  const Vector2& upperRight() const;

  std::array<Vector2, 4> vertices() const;

  /// Get the tolerance.
  const BoundaryTolerance& tolerance() const;

  /// Check if the point is inside the box.
  bool inside(const Vector2& point,
              const std::optional<SquareMatrix2>& jacobian) const;

 private:
  Vector2 m_lowerLeft;
  Vector2 m_upperRight;
  BoundaryTolerance m_tolerance;
};

template <typename Vector2Container>
class PolygonBoundaryCheck {
 public:
  PolygonBoundaryCheck(const Vector2Container& vertices,
                       BoundaryTolerance tolerance);

  /// Get the vertices of the polygon.
  const Vector2Container& vertices() const;

  /// Get the tolerance.
  const BoundaryTolerance& tolerance() const;

  /// Check if the point is inside the polygon.
  bool inside(const Vector2& point,
              const std::optional<SquareMatrix2>& jacobian) const;

 private:
  const Vector2Container& m_vertices;
  BoundaryTolerance m_tolerance;
};

}  // namespace Acts

#include "Acts/Surfaces/BoundaryCheck.ipp"
