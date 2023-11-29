// This file is part of the Acts project.
//
// Copyright (C) 2022-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Portal.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/ThrowAssert.hpp"

#include <cstddef>
#include <stdexcept>
#include <utility>

namespace Acts {
namespace Experimental {
class DetectorVolume;
}  // namespace Experimental
}  // namespace Acts

Acts::Experimental::Portal::Portal(std::shared_ptr<RegularSurface> surface)
    : m_surface(std::move(surface)) {
  throw_assert(m_surface, "Portal surface is nullptr");
}

std::shared_ptr<Acts::Experimental::Portal>
Acts::Experimental::Portal::makeShared(
    std::shared_ptr<RegularSurface> surface) {
  return std::shared_ptr<Portal>(new Portal(std::move(surface)));
}

const Acts::RegularSurface& Acts::Experimental::Portal::surface() const {
  return *m_surface.get();
}

Acts::RegularSurface& Acts::Experimental::Portal::surface() {
  return *m_surface.get();
}

const Acts::Experimental::Portal::DetectorVolumeUpdaters&
Acts::Experimental::Portal::detectorVolumeUpdaters() const {
  return m_volumeUpdaters;
}

Acts::Experimental::Portal::AttachedDetectorVolumes&
Acts::Experimental::Portal::attachedDetectorVolumes() {
  return m_attachedVolumes;
}

std::shared_ptr<Acts::Experimental::Portal>
Acts::Experimental::Portal::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::Portal>
Acts::Experimental::Portal::getSharedPtr() const {
  return shared_from_this();
}

void Acts::Experimental::Portal::assignGeometryId(
    const GeometryIdentifier& geometryId) {
  m_surface->assignGeometryId(geometryId);
}

void Acts::Experimental::Portal::fuse(std::shared_ptr<Portal>& other) {
  Direction bDir = Direction::Backward;

  // Determine this directioon
  Direction tDir = (!m_volumeUpdaters[bDir.index()].connected())
                       ? Direction::Forward
                       : Direction::Backward;

  if (!m_volumeUpdaters[tDir.index()].connected()) {
    throw std::invalid_argument(
        "Portal: trying to fuse portal (keep) with no links.");
  }
  // And now check other direction
  Direction oDir = tDir.invert();
  if (!other->m_volumeUpdaters[oDir.index()].connected()) {
    throw std::runtime_error(
        "Portal: trying to fuse portal (waste) with no links.");
  }

  if (m_surface->surfaceMaterial() != nullptr &&
      other->surface().surfaceMaterial() != nullptr) {
    throw std::runtime_error(
        "Portal: both surfaces have surface material, fusing will lead to "
        "information loss.");
  } else if (other->surface().surfaceMaterial() != nullptr) {
    m_surface->assignSurfaceMaterial(
        other->surface().surfaceMaterialSharedPtr());
  }

  auto odx = oDir.index();
  m_volumeUpdaters[odx] = std::move(other->m_volumeUpdaters[odx]);
  m_attachedVolumes[odx] = other->m_attachedVolumes[odx];
  // And finally overwrite the original portal
  other = getSharedPtr();
}

void Acts::Experimental::Portal::assignDetectorVolumeUpdater(
    Direction dir, DetectorVolumeUpdater dVolumeUpdater,
    std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes) {
  auto idx = dir.index();
  m_volumeUpdaters[idx] = std::move(dVolumeUpdater);
  m_attachedVolumes[idx] = std::move(attachedVolumes);
}

void Acts::Experimental::Portal::assignDetectorVolumeUpdater(
    DetectorVolumeUpdater dVolumeUpdater,
    std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes) {
  // Check and throw exceptions
  if (!m_volumeUpdaters[0u].connected() && !m_volumeUpdaters[1u].connected()) {
    throw std::runtime_error("Portal: portal has no link on either side.");
  }
  if (m_volumeUpdaters[0u].connected() && m_volumeUpdaters[1u].connected()) {
    throw std::runtime_error("Portal: portal already has links on both sides.");
  }
  std::size_t idx = m_volumeUpdaters[0u].connected() ? 1u : 0u;
  m_volumeUpdaters[idx] = std::move(dVolumeUpdater);
  m_attachedVolumes[idx] = std::move(attachedVolumes);
}

void Acts::Experimental::Portal::updateDetectorVolume(
    const GeometryContext& gctx, NavigationState& nState) const {
  const auto& position = nState.position;
  const auto& direction = nState.direction;
  const Vector3 normal = surface().normal(gctx, position);
  Direction dir = Direction::fromScalar(normal.dot(direction));
  const auto& vUpdater = m_volumeUpdaters[dir.index()];
  if (vUpdater.connected()) {
    vUpdater(gctx, nState);
  } else {
    nState.currentVolume = nullptr;
  }
}
