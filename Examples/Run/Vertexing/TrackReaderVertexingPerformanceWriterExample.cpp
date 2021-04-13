// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file Find vertices using truth particle information as input
///
/// Reads truth particles from TrackMl files and use the truth information
/// to generate smeared track parameters. Use this pseudo-reconstructed
/// tracks as the input to the vertex finder.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Printers/TrackParametersPrinter.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryParametersReader.hpp"
#include "ActsExamples/Io/Root/RootVertexPerformanceWriter.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsExamples/TruthTracking/TrackSelector.hpp"

#include <memory>

using namespace Acts::UnitLiterals;
using namespace ActsExamples;

int main(int argc, char* argv[]) {
  // setup and parse options
  auto desc = Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addVertexingOptions(desc);
  Options::addInputOptions(desc);
  Options::addMagneticFieldOptions(desc);
  Options::addOutputOptions(desc, OutputFormat::DirectoryOnly);
  auto vars = Options::parse(desc, argc, argv);
  if (vars.empty()) {
    return EXIT_FAILURE;
  }

  // basic setup
  auto logLevel = Options::readLogLevel(vars);
  auto rnd =
      std::make_shared<RandomNumbers>(Options::readRandomNumbersConfig(vars));
  Sequencer sequencer(Options::readSequencerConfig(vars));

  auto outputDir = ensureWritableDirectory(vars["output-dir"].as<std::string>());

  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vars);

  RootTrajectoryParametersReader::Config trackParamsReader;
  trackParamsReader.trackCollection = "fittedTrackParameters";
  trackParamsReader.particleCollection = "truthParticles";
  trackParamsReader.inputFile = "./data/reco_generic_new/ttbar_mu200/trackparams_fitter.root";
  sequencer.addReader(std::make_shared<RootTrajectoryParametersReader>(
      trackParamsReader));

  // Apply some primary vertexing selection cuts
  TrackSelector::Config trackSelectorConfig;
  trackSelectorConfig.inputTrackParameters = trackParamsReader.trackCollection;
  trackSelectorConfig.outputTrackParameters = "selectedTracks";
  trackSelectorConfig.outputTrackIndices = "outputTrackIndices";
  trackSelectorConfig.removeNeutral = true;
  trackSelectorConfig.absEtaMax = vars["vertexing-eta-max"].as<double>();
  trackSelectorConfig.loc0Max = vars["vertexing-rho-max"].as<double>() * 1_mm;
  trackSelectorConfig.ptMin = vars["vertexing-pt-min"].as<double>() * 1_MeV;
  sequencer.addAlgorithm(
      std::make_shared<TrackSelector>(trackSelectorConfig, logLevel));

  
  // find vertices
  AdaptiveMultiVertexFinderAlgorithm::Config findVertices(magneticField);
  findVertices.inputTrackParameters = trackSelectorConfig.outputTrackParameters;
  findVertices.outputProtoVertices = "fittedProtoVertices";
  findVertices.outputVertices = "fittedVertices";
  sequencer.addAlgorithm(std::make_shared<AdaptiveMultiVertexFinderAlgorithm>(
      findVertices, logLevel));

  // write track parameters from fitting
  RootVertexPerformanceWriter::Config vertexWriterConfig;
  vertexWriterConfig.inputTruthParticles = trackParamsReader.particleCollection;
  vertexWriterConfig.inputVertices = findVertices.outputVertices;
  vertexWriterConfig.outputFilename = "vertexperformance_AMVF.root";
  vertexWriterConfig.outputTreename = "amvf";
  vertexWriterConfig.outputDir = outputDir;
  sequencer.addWriter(std::make_shared<RootVertexPerformanceWriter>(
      vertexWriterConfig, logLevel));


  return sequencer.run();
}
