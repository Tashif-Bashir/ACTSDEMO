// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Seeding/SeedingAlgorithm.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Utilities/Units.hpp>

#include <memory>

#include <boost/program_options.hpp>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // --------------------------------------------------------------------------------
  GenericDetector detector;

  // Setup and parse options

  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addOutputOptions(desc);
  Options::addInputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Read particles (initial states) and clusters from CSV files
  auto particleReader = Options::readCsvParticleReaderConfig(vm);
  particleReader.inputStem = "particles_initial";
  particleReader.outputParticles = "particles_initial";
  sequencer.addReader(
      std::make_shared<CsvParticleReader>(particleReader, logLevel));

  // Read truth hits from CSV files
  auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vm);
  simHitReaderCfg.inputStem = "simhits";
  simHitReaderCfg.outputSimHits = "simhits";
  sequencer.addReader(
      std::make_shared<CsvSimHitReader>(simHitReaderCfg, logLevel));

  // Create smeared measurements
  HitSmearing::Config hitSmearingCfg;
  hitSmearingCfg.inputSimHits = simHitReaderCfg.outputSimHits;
  hitSmearingCfg.outputMeasurements = "measurements";
  hitSmearingCfg.outputSourceLinks = "sourcelinks";
  hitSmearingCfg.outputMeasurementParticlesMap = "measurement_particles_map";
  hitSmearingCfg.outputMeasurementSimHitsMap = "measurement_simhits_map";
  hitSmearingCfg.sigmaLoc0 = 25_um;
  hitSmearingCfg.sigmaLoc1 = 100_um;
  hitSmearingCfg.randomNumbers = rnd;
  hitSmearingCfg.trackingGeometry = tGeometry;
  sequencer.addAlgorithm(
      std::make_shared<HitSmearing>(hitSmearingCfg, logLevel));

  // Seeding algorithm
  SeedingAlgorithm::Config seeding;
  seeding.outputSeeds = "seeds";
  seeding.inputMeasurements = hitSmearingCfg.outputMeasurements;
  sequencer.addAlgorithm(std::make_shared<SeedingAlgorithm>(seeding, logLevel));

  // Performance Writer
  SeedingPerformanceWriter::Config seedPerfCfg;
  seedPerfCfg.inputSeeds = seeding.outputSeeds;
  seedPerfCfg.inputParticles = particleReader.outputParticles;
  seedPerfCfg.inputMeasurementParticlesMap = hitSmearingCfg.outputMeasurementParticlesMap;
  seedPerfCfg.outputFilename = "performance.root";
  sequencer.addWriter(
      std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));

  return sequencer.run();
}
