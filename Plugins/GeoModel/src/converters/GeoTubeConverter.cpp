// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/converters/GeoTubeConverter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/Units.h>

Acts::Result<Acts::GeoModelSensitiveSurface>
Acts::GeoTubeConverter::toSensitiveSurface(const GeoFullPhysVol& geoFPV) const {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();
    auto geoTube = dynamic_cast<const GeoTube*>(geoShape);
    if (geoTube != nullptr) {
      return Result<GeoModelSensitiveSurface>::success(
          toSurface(geoFPV, *geoTube, transform, true));
    }
    return Result<GeoModelSensitiveSurface>::failure(
        GeoModelConversionError::WrongShapeForConverter);
  }
  return Result<GeoModelSensitiveSurface>::failure(
      GeoModelConversionError::MissingLogicalVolume);
}

Acts::Result<std::shared_ptr<Acts::Surface>>
Acts::GeoTubeConverter::toPassiveSurface(const GeoFullPhysVol& geoFPV) const {
  // Retrieve logcal volume and absolute transform
  const GeoLogVol* logVol = geoFPV.getLogVol();
  const Transform3& transform = geoFPV.getAbsoluteTransform(nullptr);
  if (logVol != nullptr) {
    const GeoShape* geoShape = logVol->getShape();

    auto geoTube = dynamic_cast<const GeoTube*>(geoShape);
    if (geoTube != nullptr) {
      // Conversion function call with sensitive = false
      auto [element, surface] = toSurface(geoFPV, *geoTube, transform, false);
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
Acts::GeoTubeConverter::toSurface(const GeoFullPhysVol& geoFPV,
                                  const GeoTube& geoTube,
                                  const Transform3& absTransform,
                                  bool sensitive) const {
  /// auto-calculate the unit length conversion
  static constexpr ActsScalar unitLength =
      Acts::UnitConstants::mm / GeoModelKernelUnits::millimeter;

  // Create the surface transform
  Transform3 transform = Transform3::Identity();
  transform.translation() = unitLength * absTransform.translation();
  transform.linear() = absTransform.rotation();

  // Create the surface
  ActsScalar innerRadius = unitLength * geoTube.getRMin();
  ActsScalar outerRadius = unitLength * geoTube.getRMax();
  ActsScalar halfZ = unitLength * geoTube.getZHalfLength();

  if (targetShape == Surface::SurfaceType::Straw) {
    // Create the element and the surface
    auto lineBounds = std::make_shared<LineBounds>(outerRadius, halfZ);
    if (!sensitive) {
      auto surface = Surface::makeShared<StrawSurface>(transform, lineBounds);
      return std::make_tuple(nullptr, surface);
    }

    auto detectorElement =
        GeoModelDetectorElement::createDetectorElement<StrawSurface>(
            geoFPV, lineBounds, transform, 2 * outerRadius);
    auto surface = detectorElement->surface().getSharedPtr();
    return std::make_tuple(detectorElement, surface);
    // Next option is tranlsation to disc
  } else if (targetShape == Surface::SurfaceType::Disc) {
    auto radialBounds =
        std::make_shared<RadialBounds>(innerRadius, outerRadius);
    if (!sensitive) {
      auto surface = Surface::makeShared<DiscSurface>(transform, radialBounds);
      return std::make_tuple(nullptr, surface);
    }

    // Create the element and the surface
    auto detectorElement =
        GeoModelDetectorElement::createDetectorElement<DiscSurface>(
            geoFPV, radialBounds, transform, 2 * halfZ);
    auto surface = detectorElement->surface().getSharedPtr();
    return std::make_tuple(detectorElement, surface);
  }
  // Finally cylinder to cylinder
  auto cylinderBounds = std::make_shared<CylinderBounds>(outerRadius, halfZ);
  if (!sensitive) {
    auto surface =
        Surface::makeShared<CylinderSurface>(transform, cylinderBounds);
    return std::make_tuple(nullptr, surface);
  }
  // Create the element and the surface
  auto detectorElement =
      GeoModelDetectorElement::createDetectorElement<CylinderSurface>(
          geoFPV, cylinderBounds, transform, outerRadius - innerRadius);
  auto surface = detectorElement->surface().getSharedPtr();
  return std::make_tuple(detectorElement, surface);
}
