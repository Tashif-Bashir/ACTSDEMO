// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include "Acts/Geometry/HierarchicalGeometryContainer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <nlohmann/json.hpp>

#include <map>

namespace Acts {

/// @class JsonGeometryConverter
///
/// @brief read the material from Json
template <typename object_t>
class JsonHierarchicalObjectConverter {
 public:
  /// @class Config
  /// Configuration of the Reader
  class Config {
   public:
    /// The geometry version
    std::string geoversion = "undefined";
    /// The detector tag
    std::string detkey = "detector";
    /// The volume identification string
    std::string volkey = "volumes";
    /// The boundary surface string
    std::string boukey = "boundaries";
    /// The layer identification string
    std::string laykey = "layers";
    /// The approach identification string
    std::string appkey = "approach";
    /// The sensitive identification string
    std::string senkey = "sensitive";
    /// The representing idntification string
    std::string repkey = "representing";
    /// The name identification
    std::string namekey = "Name";
    /// The object key
    std::string datakey = "Object";

    /// The default logger
    std::shared_ptr<const Logger> logger;
    /// The name of the writer
    std::string name = "";

    /// Steering to handle sensitive data
    bool processSensitives = true;
    /// Steering to handle approach data
    bool processApproaches = true;
    /// Steering to handle representing data
    bool processRepresenting = true;
    /// Steering to handle boundary data
    bool processBoundaries = true;
    /// Steering to handle volume data
    bool processVolumes = true;

    /// Constructor
    ///
    /// @param lname Name of the writer tool
    /// @param lvl The output logging level
    Config(const std::string& lname = "JsonHierarchicalObjectConverter",
           Logging::Level lvl = Logging::INFO)
        : logger(getDefaultLogger(lname, lvl)), name(lname) {}
  };

  using geo_id_value = uint64_t;
  using Representation = std::map<geo_id_value, object_t>;

  /// @brief Layer representation for Json writing
  struct LayerRep {
    // the layer id
    GeometryID layerID;

    Representation sensitives;
    Representation approaches;
    object_t representing;
  };

  /// @brief Volume representation for Json writing
  struct VolumeRep {
    // The geometry id
    GeometryID volumeID;

    /// The namne
    std::string volumeName;

    std::map<geo_id_value, LayerRep> layers;
    Representation boundaries;
    object_t volume;
  };

  /// @brief Detector representation for Json writing
  struct DetectorRep {
    std::map<geo_id_value, VolumeRep> volumes;
  };

  /// Constructor
  ///
  /// @param cfg configuration struct for the reader
  JsonHierarchicalObjectConverter(const Config& cfg);

  /// Destructor
  ~JsonHierarchicalObjectConverter() = default;

  /// Convert method
  ///
  /// @param surfaceMaterialMap The indexed material map collection
  HierarchicalGeometryContainer<object_t> jsonToHierarchicalContainer(
      const nlohmann::json& map,
      std::function<object_t(const nlohmann::json&)> fromJson) const;

  /// Convert method
  ///
  /// @param surfaceMaterialMap The indexed material map collection
  nlohmann::json hierarchicalObjectToJson(
      const HierarchicalGeometryContainer<object_t>& hObject,
      std::function<nlohmann::json(const object_t&)> toJson) const;

  /// Write method
  ///
  /// @param tGeometry is the tracking geometry which contains the material
  nlohmann::json trackingGeometryToJson(
      const TrackingGeometry& tGeometry,
      std::function<nlohmann::json(const object_t&)> toJson,
      std::function<object_t(const GeometryID&)> initialise) const;

 private:
  /// Convert to internal representation method, recursive call
  ///
  /// @param tGeometry is the tracking geometry which contains the material
  void convertToRep(
      DetectorRep& detRep, const TrackingVolume& tVolume,
      std::function<object_t(const GeometryID&)> initialise) const;

  /// Convert to internal representation method
  ///
  /// @param tGeometry is the tracking geometry which contains the material
  LayerRep convertToRep(
      const Layer& tLayer,
      std::function<object_t(const GeometryID&)> initialise) const;

  /// Convert to internal representation method
  ///
  /// @param tGeometry is the tracking geometry which contains the material
  void convertToRep(
      DetectorRep& detRep,
      const Acts::HierarchicalGeometryContainer<object_t>& hObject) const;

  /// Create Json from a detector represenation
  nlohmann::json detectorRepToJson(
      const DetectorRep& detRep,
      std::function<nlohmann::json(const object_t&)> toJson) const;

  /// The config class
  Config m_cfg;

  /// Private access to the logging instance
  const Logger& logger() const { return *m_cfg.logger; }
};

}  // namespace Acts
#include "Acts/Plugins/Json/JsonHierarchicalObjectConverter.ipp"
