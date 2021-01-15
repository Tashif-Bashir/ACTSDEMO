// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/JsonUtilitiesConverter.hpp"

// Custom Json encoder/decoders. Maming is mandated by nlohman::json and thus
// can not match our naming guidelines.

void to_json(nlohmann::json& j, const Acts::BinningData& bd) {
  // Common to all bin utilities
  j["min"] = bd.min;
  j["max"] = bd.max;
  j["option"] = (bd.option == Acts::open ? "open" : "closed");
  j["value"] = Acts::binningValueNames[bd.binvalue];
  int bins = bd.bins();
  // Write sub bin data if present
  if (bd.subBinningData != nullptr) {
    nlohmann::json subjson;
    to_json(subjson, *bd.subBinningData);
    j["subdata"] = subjson;
    j["subadditive"] = bd.subBinningAdditive;
    // this modifies the bins as bins() returns total number in general
    if (bd.subBinningAdditive){
      bins -= static_cast<int>(subjson["bins"]) + 1;
    } else {
      bins /= static_cast<int>(subjson["bins"]);
    }
  }
  // Now distinguish between equidistant / arbitrary
  if (bd.type == Acts::equidistant) {
    j["type"] = "equidistant";
  } else if (bd.type == Acts::arbitrary) {
    j["type"] = "arbitrary";
    j["boundaries"] = bd.boundaries();
  }
  j["bins"] = bins;
}

void from_json(const nlohmann::json& j, Acts::BinningData& bd) {
  // Common to all bin utilities
  float min = j["min"];
  float max = j["max"];
  int bins = j["bins"];
  std::string valueName = j["value"];
  auto valueIter = std::find(Acts::binningValueNames.begin(),
                             Acts::binningValueNames.end(), valueName);
  Acts::BinningValue bValue = static_cast<Acts::BinningValue>(
      valueIter - Acts::binningValueNames.begin());
  if (bins == 1) {
    bd = Acts::BinningData(bValue, min, max);
    return;
  }
  Acts::BinningOption bOption =
      (j["option"] == "open") ? Acts::open : Acts::closed;
  Acts::BinningType bType =
      (j["type"] == "equidistant") ? Acts::equidistant : Acts::arbitrary;

  std::unique_ptr<Acts::BinningData> subBinning = nullptr;
  bool subBinningAdditive = false;
  if (j.find("subdata") != j.end()) {
    subBinningAdditive = j["subadditive"];
  }

  if (bType == Acts::equidistant) {
    bd = Acts::BinningData(bOption, bValue, bins, min, max,
                           std::move(subBinning), subBinningAdditive);
  } else {
    std::vector<float> boundaries = j["boundaries"];
    bd = Acts::BinningData(bOption, bValue, boundaries, std::move(subBinning));
  }
}

void to_json(nlohmann::json& j, const Acts::BinUtility& bu) {
  for (const auto& bdata : bu.binningData()) {
    nlohmann::json jdata;
    to_json(jdata, bdata);
    j.push_back(jdata);
  }
}

void from_json(const nlohmann::json& j, Acts::BinUtility& bu) {
  for (const auto& jdata : j) {
    Acts::BinningData bd;
    from_json(jdata, bd);
    bu += Acts::BinUtility(bd);
  }
}
