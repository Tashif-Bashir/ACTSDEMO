// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>

#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <mutex>
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

namespace ActsExamples {

class WhiteBoard;
namespace DD4hepG4 {
class DD4hepToG4Svc;
}

class EventRecording final : public ActsExamples::BareAlgorithm {
 public:
  /// @class Config
  struct Config {
    std::string eventInput = "";
    std::string eventOutput = "geant-outcome-tracks";

    std::unique_ptr<G4VUserDetectorConstruction> detectorConstruction = nullptr;

    /// random number seed 1
    int seed1 = 12345;
    /// random number seed 2
    int seed2 = 45678;
  };

  /// Constructor
  EventRecording(Config&& cnf,
                 Acts::Logging::Level level = Acts::Logging::INFO);
  ~EventRecording() { m_runManager = nullptr; };

  ActsExamples::ProcessCode execute(
      const AlgorithmContext& context) const final override;

 private:
  /// The config object
  Config m_cfg;
  /// G4 run manager
  std::unique_ptr<G4RunManager> m_runManager;

  // has to be mutable; algorithm interface enforces object constness
  mutable std::mutex m_runManagerLock;
};
}  // namespace ActsExamples
