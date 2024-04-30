// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/GridAccessHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>
#include <vector>

class G4VPhysicalVolume;

namespace Acts {
class Surface;
}

namespace ActsExamples {

// Helper struct to find the sensitive surface candidates
struct SensitiveCandidates {
  std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
  /// Find the sensitive surfaces for a given position
  ///
  /// This fullfills the concept of a SensitiveCandidates
  ///
  /// @param gctx the geometry context
  /// @param position the position to look for sensitive surfaces
  ///
  /// @return a vector of sensitive surfaces
  std::vector<const Acts::Surface*> operator()(
      const Acts::GeometryContext& gctx, const Acts::Vector3& position) const {
    std::vector<const Acts::Surface*> surfaces;

    if (trackingGeometry != nullptr) {
      auto layer = trackingGeometry->associatedLayer(gctx, position);

      if (layer->surfaceArray() != nullptr) {
        for (const auto& surface : layer->surfaceArray()->surfaces()) {
          if (surface->associatedDetectorElement() != nullptr) {
            surfaces.push_back(surface);
          }
        }
      }
    }
    return surfaces;
  }
};

/// This Mapper takes a (non-const) Geant4 geometry and maps
/// it such that name will be containing the mapping prefix
/// and the Acts::GeometryIdentifier of the surface.
///
/// The mapping is done by matching the position of the G4 physical volume with
/// the center position of an Acts::Surface.
///
/// This allows to directly associate Geant4 hits to the sensitive
/// elements of the Acts::TrackingGeoemtry w/o map lookup.
class SensitiveSurfaceMapper {
 public:
  using SensitiveCandidates = std::function<std::vector<const Acts::Surface*>(
      const Acts::GeometryContext&, const Acts::Vector3&)>;

  constexpr static std::string_view mappingPrefix = "ActsGeoID#";

  /// Configuration struct for the surface mapper
  struct Config {
    /// For which G4 material names we try to find a mapping
    std::vector<std::string> materialMappings;

    /// For which G4 volume names we try to find a mapping
    std::vector<std::string> volumeMappings;

    /// The sensitive surfaces that are being mapped to
    SensitiveCandidates candidateSurfaces;
  };

  /// Constructor with:
  ///
  /// @param cfg the configuration struct
  /// @param logger the logging instance
  SensitiveSurfaceMapper(const Config& cfg,
                         std::unique_ptr<const Acts::Logger> logger =
                             Acts::getDefaultLogger("SensitiveSurfaceMapper",
                                                    Acts::Logging::INFO));
  ~SensitiveSurfaceMapper() = default;

  /// Recursive mapping function that walks through the Geant4
  /// hierarchy and applies name remapping to the Physical volumes
  /// of the Geant4 geometry.
  ///
  /// @param gctx the geometry context
  /// @param g4PhysicalVolume the current physical volume in process
  /// @param motherPosition the absolute position of the mother
  /// @param sCounter  a counter of how many volumes have been remapped
  void remapSensitiveNames(const Acts::GeometryContext& gctx,
                           G4VPhysicalVolume* g4PhysicalVolume,
                           const Acts::Transform3& motherTransform,
                           int& sCounter) const;

 protected:
  /// Configuration object
  Config m_cfg;

 private:
  /// Private access method to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The looging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples
