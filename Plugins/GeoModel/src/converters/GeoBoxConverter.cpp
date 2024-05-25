// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/converters/GeoBoxConverter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/Units.h>

Acts::Result<Acts::GeoModelSensitiveSurface>
Acts::GeoBoxConverter::toSensitiveSurface(const GeoFullPhysVol& geoFPV) const {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();
    auto geoBox = dynamic_cast<const GeoBox*>(geoShape);
    if (geoBox != nullptr) {
      return Result<GeoModelSensitiveSurface>::success(
          toSurface(geoFPV, *geoBox, transform, true));
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }
  return Result<GeoModelSensitiveSurface>::failure(
      GeoModelConversionError::MissingLogicalVolume);
}

Acts::Result<std::shared_ptr<Acts::Surface>>
Acts::GeoBoxConverter::toPassiveSurface(const GeoFullPhysVol& geoFPV) const {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();

    auto geoBox = dynamic_cast<const GeoBox*>(geoShape);
    if (geoBox != nullptr) {
      // Conversion function call with sensitive = false
      auto [element, surface] = toSurface(geoFPV, *geoBox, transform, false);
      return Result<std::shared_ptr<Surface>>::success(surface);
    }
    return Result<std::shared_ptr<Surface>>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }
  return Result<std::shared_ptr<Surface>>::failure(
      GeoModelConversionError::MissingLogicalVolume);
}

std::tuple<std::shared_ptr<Acts::GeoModelDetectorElement>,
           std::shared_ptr<Acts::Surface>>
Acts::GeoBoxConverter::toSurface(const GeoFullPhysVol& geoFPV,
                                 const GeoBox& geoBox,
                                 const Transform3& absTransform,
                                 bool sensitive) const {
  /// auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / GeoModelKernelUnits::millimeter;

  // Create the surface transform
  Transform3 transform = Transform3::Identity();
  transform.translation() = unitLength * absTransform.translation();
  auto rotation = absTransform.rotation();
  // Get the half lengths
  std::vector<ActsScalar> halfLengths = {geoBox.getXHalfLength(),
                                         geoBox.getYHalfLength(),
                                         geoBox.getZHalfLength()};
  // Create the surface
  auto minElement = std::min_element(halfLengths.begin(), halfLengths.end());
  auto zIndex = std::distance(halfLengths.begin(), minElement);
  std::size_t yIndex = zIndex > 0u ? zIndex - 1u : 2u;
  std::size_t xIndex = yIndex > 0u ? yIndex - 1u : 2u;

  Vector3 colX = rotation.col(xIndex);
  Vector3 colY = rotation.col(yIndex);
  Vector3 colZ = rotation.col(zIndex);
  rotation.col(0) = colX;
  rotation.col(1) = colY;
  rotation.col(2) = colZ;
  transform.linear() = rotation;

  // Create the surface bounds
  ActsScalar halfX = unitLength * halfLengths[xIndex];
  ActsScalar halfY = unitLength * halfLengths[yIndex];
  auto rectangleBounds = std::make_shared<Acts::RectangleBounds>(halfX, halfY);
  if (!sensitive) {
    auto surface =
        Surface::makeShared<PlaneSurface>(transform, rectangleBounds);
    return std::make_tuple(nullptr, surface);
  }
  // Create the element and the surface
  auto detectorElement =
      GeoModelDetectorElement::createDetectorElement<PlaneSurface>(
          geoFPV, rectangleBounds, transform,
          2 * unitLength * halfLengths[zIndex]);
  auto surface = detectorElement->surface().getSharedPtr();
  // Return the detector element and surface
  return std::make_tuple(detectorElement, surface);
}
