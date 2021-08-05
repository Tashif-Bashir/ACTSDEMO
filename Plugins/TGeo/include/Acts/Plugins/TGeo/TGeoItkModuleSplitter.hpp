// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/TGeo/ITGeoDetectorElementSplitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

class TGeoNode;

namespace Acts {

class TGeoDetectorElement;

/// @brief TGeoItkModuleSplitter
///
/// Split Itk modules into submodules, depending on the sensor type
class TGeoItkModuleSplitter : public ITGeoDetectorElementSplitter {
  public:
  /// Nested configuration struct
  struct Config {
    // Map the nodes name to the splitting parameters
    std::map<std::string, unsigned int> paramMap = {};
  };

  /// Constructor
  ///
  /// @param cfg the configuration struct
  /// @param logger the logging object
  TGeoItkModuleSplitter(const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "TGeoItkModuleSplitter", Acts::Logging::INFO));

  virtual ~TGeoItkModuleSplitter() = default;

  /// Take a geometry context and TGeoElement and find the correct splitting 
  /// method for the module type.
  ///
  /// @param gctx is a geometry context object
  /// @param detElement is a TGeoDetectorElement that is eventually split
  ///
  /// @note If no split is performed the unsplit detector element is returned
  ///
  /// @return a vector of TGeoDetectorElement objects
  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> split(
      const GeometryContext& gctx,
      std::shared_ptr<const Acts::TGeoDetectorElement> detElement) const;

 private:

  /// Take a geometry context and TGeoElement in the Itk barrel region 
  /// and split it into sub elements.
  ///
  /// @param gctx is a geometry context object
  /// @param detElement is a TGeoDetectorElement that is eventually split
  ///
  /// @note If no split is performed the unsplit detector element is returned
  ///
  /// @return a vector of TGeoDetectorElement objects
  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> splitBarrelModule(
      const GeometryContext& gctx,
      std::shared_ptr<const Acts::TGeoDetectorElement> detElement,
      unsigned int nSegments) const;

  /// Take a geometry context and TGeoElement in the Itk disks and split it 
  /// into sub elements.
  ///
  /// @param gctx is a geometry context object
  /// @param detElement is a TGeoDetectorElement that is eventually split
  ///
  /// @note If no split is performed the unsplit detector element is returned
  ///
  /// @return a vector of TGeoDetectorElement objects
  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> 
  splitDiskModule(
      const GeometryContext& gctx,
      std::shared_ptr<const Acts::TGeoDetectorElement> detElement,
      unsigned int nSegments) const;

  /// Contains the splitting parameters, sorted by sensor type
  Config m_cfg;

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts