// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorVolumeBuilder.hpp"

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"

#include <stdexcept>

Acts::Experimental::DetectorVolumeBuilder::DetectorVolumeBuilder(
    const Acts::Experimental::DetectorVolumeBuilder::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : IDetectorComponentBuilder(), m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.externalsBuilder == nullptr) {
    throw std::invalid_argument(
        "DetectorVolumeBuilder: no external structure builder defined.");
  }
}

Acts::Experimental::DetectorComponent
Acts::Experimental::DetectorVolumeBuilder::construct(
    const GeometryContext& gctx) const {
  // The outgoing root volumes
  std::vector<std::shared_ptr<DetectorVolume>> rootVolumes;

  // Screen printout of the auxilliary information
  if (not m_cfg.auxilliary.empty()) {
    ACTS_DEBUG(m_cfg.auxilliary);
  }
  ACTS_DEBUG("Building a volume with name " << m_cfg.name);

  // Get transform and bounds from the volume
  auto [transform, bounds, portalGenerator] =
      m_cfg.externalsBuilder->construct(gctx);

  // Although only a single volume, describe it as a
  // portal container for further processing
  DetectorComponent::PortalContainer portalContainer;
  // The detector volume to be constructed
  std::shared_ptr<DetectorVolume> dVolume = nullptr;
  // If there are no internals, the volume is fully defined
  if (m_cfg.internalsBuilder == nullptr) {
    ACTS_VERBOSE("No internal structure present.")
    // Contruct the DetectorVolume
    dVolume = DetectorVolumeFactory::construct(
        portalGenerator, gctx, m_cfg.name, transform, std::move(bounds),
        tryAllPortals());
  } else {
    // Internal structure is present
    ACTS_VERBOSE("Internal structure is being built.")
    auto [surfaces, volumes, surfacesUpdator, volumeUpdator] =
        m_cfg.internalsBuilder->construct(gctx);

    // Add the internally created volumes as root volumes
    if (m_cfg.addInternalsToRoot) {
      for (const auto& v : volumes) {
        rootVolumes.push_back(v);
      }
    }
    // Contruct the DetectorVolume
    dVolume = DetectorVolumeFactory::construct(
        portalGenerator, gctx, m_cfg.name, transform, std::move(bounds),
        surfaces, volumes, std::move(surfacesUpdator),
        std::move(volumeUpdator));
  }
  // All portals are defined and build the current shell
  for (auto [ip, p] : enumerate(dVolume->portalPtrs())) {
    portalContainer[ip] = {p};
  }
  // Add to the root volume collection if configured
  if (m_cfg.addToRoot) {
    ACTS_VERBOSE("DetectorVolume is being attached to the root volumes.");
    rootVolumes.push_back(dVolume);
  }
  // The newly built volume is the only produced volume
  return Acts::Experimental::DetectorComponent{
      {dVolume},
      portalContainer,
      RootDetectorVolumes{rootVolumes, tryRootVolumes()}};
}
