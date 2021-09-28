// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Geant4/G4DetectorConstructionFactory.hpp"
#include "ActsExamples/Geant4/Geant4Simulation.hpp"
#include "ActsExamples/Geant4/MaterialPhysicsList.hpp"
#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"
#include "ActsExamples/Geant4/SimParticleTranslation.hpp"

#include <memory>

#include <G4MagneticField.hh>
#include <G4RunManager.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4UserTrackingAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {
void addGeant4HepMC3(Context& ctx);
}

PYBIND11_MODULE(ActsPythonBindingsGeant4, mod) {
  py::class_<G4RunManager>(mod, "G4RunManager");
  py::class_<G4VUserPrimaryGeneratorAction>(mod,
                                            "G4VUserPrimaryGeneratorAction");
  py::class_<G4UserRunAction>(mod, "G4UserRunAction");
  py::class_<G4UserTrackingAction>(mod, "G4UserTrackingAction");
  py::class_<G4UserSteppingAction>(mod, "G4UserSteppingAction");
  py::class_<G4MagneticField>(mod, "G4MagneticField");
  py::class_<G4VUserDetectorConstruction>(mod, "G4VUserDetectorConstruction");

  py::class_<SensitiveSurfaceMapper, std::shared_ptr<SensitiveSurfaceMapper>>(
      mod, "SensitiveSurfaceMapper");

  {
    using Alg = Geant4Simulation;

    auto alg =
        py::class_<Alg, ActsExamples::BareAlgorithm, std::shared_ptr<Alg>>(
            mod, "GeantinoRecording")
            .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(outputSimHits);
    ACTS_PYTHON_MEMBER(outputParticlesInitial);
    ACTS_PYTHON_MEMBER(outputParticlesFinal);
    ACTS_PYTHON_MEMBER(outputMaterialTracks);
    ACTS_PYTHON_MEMBER(physicsList);
    ACTS_PYTHON_MEMBER(primaryGeneratorAction);
    ACTS_PYTHON_MEMBER(runActions);
    ACTS_PYTHON_MEMBER(eventActions);
    ACTS_PYTHON_MEMBER(trackingActions);
    ACTS_PYTHON_MEMBER(steppingActions);
    ACTS_PYTHON_MEMBER(detectorConstruction);
    ACTS_PYTHON_MEMBER(magneticField);
    ACTS_PYTHON_MEMBER(sensitiveSurfaceMapper);
    ACTS_PYTHON_STRUCT_END();
  }

  mod.def("materialRecordingConfig",
          [](Acts::Logging::Level level, G4VUserDetectorConstruction* detector,
             const std::string& inputParticles,
             const std::string& outputMaterialTracks) {
            auto* physicsList = new MaterialPhysicsList(
                Acts::getDefaultLogger("MaterialPhysicsList", level));

            // The Geant4 actions needed
            std::vector<G4UserRunAction*> runActions = {};
            std::vector<G4UserEventAction*> eventActions = {};
            std::vector<G4UserTrackingAction*> trackingActions = {};

            MaterialSteppingAction::Config mStepCfg;
            mStepCfg.excludeMaterials = {"Air", "Vacuum"};
            std::vector<G4UserSteppingAction*> steppingActions = {
                new MaterialSteppingAction(
                    mStepCfg,
                    Acts::getDefaultLogger("MaterialSteppingAction", level))};

            // Set up the Geant4 Simulation

            // Set the main Geant4 algorithm, primary generation, detector
            // construction
            Geant4Simulation::Config g4Cfg;

            // Read the particle from the generator
            SimParticleTranslation::Config g4PrCfg;
            g4PrCfg.inputParticles = inputParticles;
            g4PrCfg.forceParticle = true;
            g4PrCfg.forcedMass = 0.;
            g4PrCfg.forcedPdgCode = 999;
            // Set the material tracks at output
            g4Cfg.outputMaterialTracks = outputMaterialTracks;

            // Set the primarty generator
            g4Cfg.primaryGeneratorAction = new SimParticleTranslation(
                g4PrCfg,
                Acts::getDefaultLogger("SimParticleTranslation", level));
            g4Cfg.detectorConstruction = detector;

            g4Cfg.physicsList = physicsList;

            // Set the user actions
            g4Cfg.runActions = runActions;
            g4Cfg.eventActions = eventActions;
            g4Cfg.trackingActions = trackingActions;
            g4Cfg.steppingActions = steppingActions;

            return g4Cfg;
          });

  Acts::Python::Context ctx;
  ctx.modules["geant4"] = &mod;

  addGeant4HepMC3(ctx);
}
