// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

namespace Acts {

using namespace UnitLiterals;

namespace Experimental {
namespace detail {
/// @brief  This file contains helper methods to build common support structures
/// such as support cylinders or discs.
///
/// It allows to model those as Disc/CylinderSurface objects, but also - if
/// configured such - as approximations built from palanr surfaces
namespace SupportSurfacesHelper {

using SupportSurfaceComponents =
    std::tuple<Surface::SurfaceType, std::vector<ActsScalar>, Acts::Transform3>;

/// @brief  A support creator turns an extend into a vector of bound values
using SurfaceComponentsCreator =
    std::function<SupportSurfaceComponents(const Extent&)>;

/// @brief function descriptor for cylindrical support
struct CylindricalSupport {
  /// Offset in R
  /// - negative indicates inner support
  /// - zero is centered (not recommended)
  /// - positive indicates outer support
  ActsScalar rOffset = 0.;

  /// Clearance in z in order to make the support surfaces
  /// not touch the volume boundaries
  std::array<ActsScalar, 2u> zClearance = {1_mm, 1_mm};

  /// Clearance in phi if a sectoral support is chosen
  /// not to touch the volume boundaries
  std::array<ActsScalar, 2u> phiClearance = {0.0001_rad, 0.0001_rad};

  // Type is obviously a cylinder
  static constexpr Surface::SurfaceType type = Surface::SurfaceType::Cylinder;

  /// The support creator function
  ///
  /// @param lExtent the layer and/or volume extent
  ///
  /// @return the support surface components
  SupportSurfaceComponents operator()(const Extent& lExtent) const;
};

/// @brief function descriptor for disc-like support
struct DiscSupport {
  /// Offset in z
  /// - negative indicates support an z min
  /// - zero is centered
  /// - positive indicates support at z max
  ActsScalar zOffset = 0.;

  /// Clearance in r in order to make the support surfaces
  /// not touch the volume boundaries
  std::array<ActsScalar, 2u> rClearance = {1_mm, 1_mm};

  /// Clearance in phi if a sectoral support is chosen
  /// not to touch the volume boundaries
  std::array<ActsScalar, 2u> phiClearance = {0.0001_rad, 0.0001_rad};

  // Type is obviously a disc
  static constexpr Surface::SurfaceType type = Surface::SurfaceType::Disc;

  /// The support creator function
  ///
  /// @param lExtent the layer and/or volume extent
  ///
  /// @return the support surface components
  SupportSurfaceComponents operator()(const Extent& lExtent) const;
};

/// @brief Helper method to build planar support structure
struct RectangularSupport {
  /// Placement - the remaining loc0, loc1 are then cyclic
  BinningValue pPlacement = BinningValue::binZ;

  /// Offset in position placement
  ActsScalar pOffset = 0.;

  /// Clearance in first local direction - cyclic order
  std::array<ActsScalar, 2u> loc0Clearance = {1_mm, 1_mm};

  /// Clearance in phi if a sectoral support is chosen
  /// not to touch the volume boundaries
  std::array<ActsScalar, 2u> loc1Clearance = {1_mm, 1_mm};

  // Type is obviously a plane
  static constexpr Surface::SurfaceType type = Surface::SurfaceType::Plane;

  /// The support creator function
  ///
  /// @param lExtent the layer and/or volume extent
  ///
  /// @return the support surface components
  SupportSurfaceComponents operator()(const Extent& lExtent) const;
};

/// @brief Helper method to build cylindrical support structure
///
/// @param components are the components generated by the SurfaceComponentsCreator function
/// @param splits the number of surfaces through which the surface is approximated (1u ... cylinder)
///
/// @return a vector of surfaces that represent this support
std::vector<std::shared_ptr<Surface>> cylindricalSupport(
    const SupportSurfaceComponents& components, unsigned int splits = 1u);

/// @brief Helper method to build disc support structure
///
/// @param components are the components generated by the SurfaceComponentsCreator function
/// @param splits the number of surfaces through which the surface is approximated (1u ... disc)
///
/// @return a vector of surfaces that represent this support
std::vector<std::shared_ptr<Surface>> discSupport(
    const SupportSurfaceComponents& components, unsigned int splits = 1u);

/// @brief Helper method to build planar support structure
///
/// @param components are the components generated by the SurfaceComponentsCreator function
///
/// @return a vector of surfaces that represent this support
std::vector<std::shared_ptr<Surface>> rectangularSupport(
    const SupportSurfaceComponents& components);

/// Add support to already existing surfaces
///
/// @param layerSurfaces [in, out] the surfaces to which those are added
/// @param assignToAll [in, out] indices that are assigned to all bins in the indexing
/// @param layerExtent the externally provided layer Extent
/// @param componentCreator a function the support component creator
/// @param supportSplits the number of splits if splitting is configured
///
/// @note this modifies the layerSurfaces and toAllIndices
void addSupport(std::vector<std::shared_ptr<Surface>>& layerSurfaces,
                std::vector<std::size_t>& assignToAll,
                const Extent& layerExtent,
                const SurfaceComponentsCreator& componentCreator,
                unsigned int supportSplits = 1u);

}  // namespace SupportSurfacesHelper
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
