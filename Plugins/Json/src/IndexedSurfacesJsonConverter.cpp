// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/IndexedSurfacesJsonConverter.hpp"

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Navigation/NavigationStateUpdators.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/IndexedGridJsonHelper.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace {

/// @brief  The generator struct
struct IndexedSurfacesGenerator {
  using value_type = std::vector<std::size_t>;

  /// @brief  Helper function to create and connect the IndexedSurfacesImpl
  ///
  /// @tparam grid_type the type of the grid, indicates also the dimension
  ///
  /// @param grid the grid object
  /// @param bv the bin value array
  /// @param transform the transform for the indexed surfaces inmplementaiton
  ///
  /// @return a connected SurfaceCandidatesUpdator object
  template <typename grid_type>
  Acts::Experimental::SurfaceCandidatesUpdator createUpdator(
      grid_type&& grid,
      const std::array<Acts::BinningValue, grid_type::DIM>& bv,
      const Acts::Transform3& transform) {
    Acts::Experimental::IndexedSurfacesImpl<grid_type> indexedSurfaces(
        std::move(grid), bv, transform);

    // The portal delegate
    Acts::Experimental::AllPortalsImpl allPortals;

    // The chained delegate: indexed surfaces and all portals
    using DelegateType = Acts::Experimental::IndexedSurfacesAllPortalsImpl<
        grid_type, Acts::Experimental::IndexedSurfacesImpl>;
    auto indexedSurfacesAllPortals = std::make_unique<const DelegateType>(
        std::tie(allPortals, indexedSurfaces));

    // Create the delegate and connect it
    Acts::Experimental::SurfaceCandidatesUpdator nStateUpdator;
    nStateUpdator.connect<&DelegateType::update>(
        std::move(indexedSurfacesAllPortals));

    return nStateUpdator;
  }
};

}  // namespace

Acts::Experimental::SurfaceCandidatesUpdator
Acts::IndexedSurfacesJsonConverter::fromJson(
    const nlohmann::json& jSurfaceNavigation) {
  if (!jSurfaceNavigation.is_null()) {
    // The return object
    auto sfCandidates = IndexedGridJsonHelper::generateFromJson<
        Experimental::SurfaceCandidatesUpdator, IndexedSurfacesGenerator>(
        jSurfaceNavigation, "IndexedSurfaces");
    if (sfCandidates.connected()) {
      return sfCandidates;
    }
  }
  // Return the object
  return Experimental::tryAllPortalsAndSurfaces();
}
