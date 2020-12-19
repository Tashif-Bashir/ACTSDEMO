// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"

std::vector<SpacePoint> readFile(std::string filename) {
  std::string line;
  int layer;
  std::vector<SpacePoint> sps;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, varianceR, varianceZ;
      if (linetype == "lxyz") {
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        r = std::sqrt(x * x + y * y);
        float f22 = varianceR;
        float wid = varianceZ;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }
        sps.push_back(SpacePoint{x, y, z, r, layer, varianceR, varianceZ});
      }
    }
  }
  return sps;
}

int main(int argc, char** argv) {
  std::string file{"sp.txt"};
  bool help(false);
  bool quiet(false);
  bool estimateTrackParams(false);

  int opt;
  while ((opt = getopt(argc, argv, "hf:qe")) != -1) {
    switch (opt) {
      case 'f':
        file = optarg;
        break;
      case 'q':
        quiet = true;
        break;
      case 'e':
        estimateTrackParams = true;
        break;
      case 'h':
        help = true;
        [[fallthrough]];
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-hq] [-f FILENAME]\n";
        if (help) {
          std::cout << "      -h : this help" << std::endl;
          std::cout
              << "      -f FILE : read spacepoints from FILE. Default is \""
              << file << "\"" << std::endl;
          std::cout << "      -q : don't print out all found seeds"
                    << std::endl;
          std::cout << "      -e : estimate track parameters from found seeds"
                    << std::endl;
        }

        exit(EXIT_FAILURE);
    }
  }

  std::ifstream f(file);
  if (!f.good()) {
    std::cerr << "input file \"" << file << "\" does not exist\n";
    exit(EXIT_FAILURE);
  }

  auto start_read = std::chrono::system_clock::now();
  std::vector<SpacePoint> sps = readFile(file);
  std::vector<const SpacePoint*> spVec;
  spVec.reserve(sps.size());
  // Store the pointers in the vector
  std::transform(sps.begin(), sps.end(), std::back_inserter(spVec),
                 [](const SpacePoint& sp) { return &sp; });
  auto end_read = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_read = end_read - start_read;

  std::cout << "read " << spVec.size() << " SP from file " << file << " in "
            << elapsed_read.count() << "s" << std::endl;

  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));
  Acts::Seedfinder<SpacePoint> a(config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2 {
    return {sp.varianceR, sp.varianceZ};
  };

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  auto start = std::chrono::system_clock::now();
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(a.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  std::cout << "Number of regions: " << seedVector.size() << std::endl;
  int numSeeds = 0;
  for (auto& regionVec : seedVector) {
    numSeeds += regionVec.size();
  }
  std::cout << "Number of seeds generated: " << numSeeds << std::endl;
  if (!quiet) {
    for (size_t ir = 0; ir < seedVector.size(); ++ir) {
      for (size_t is = 0; is < seedVector[ir].size(); ++is) {
        const auto& seed = seedVector[ir][is];
        const SpacePoint* sp = seed.sp()[0];
        std::cout << "Print out info of found seed: " << is
                  << " in region: " << ir << std::endl;
        std::cout << sp->layer << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed.sp()[1];
        std::cout << sp->layer << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed.sp()[2];
        std::cout << sp->layer << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        std::cout << std::endl;
        if (estimateTrackParams) {
          auto transParamsRes = Acts::estimateTrackParamsFromSeed(seed.sp());
          if (transParamsRes.has_value()) {
            auto transParams = transParamsRes.value();
            std::cout << "Estimated (d0, curvature, phi): \n " << transParams[0]
                      << "  " << transParams[1] << "  " << transParams[2]
                      << std::endl;
          } else {
            std::cout
                << "WARNING: estimation of (d0, curvature, phi) for found seed "
                << is << " in region " << ir << " failed." << std::endl;
          }
          auto boundParamsRes = Acts::estimateTrackParamsFromSeed(
              seed.sp(), Acts::Transform3::Identity(), config.bFieldInZ * 1000,
              config.minPt * 0.001);
          if (boundParamsRes.has_value()) {
            auto boundParams = boundParamsRes.value();
            std::cout << "Estimated (loc0, loc1, phi, theta, q/p, t) on the "
                         "surface with identify transform: \n"
                      << boundParams.transpose() << std::endl;
          } else {
            std::cout
                << "WARNING: estimation of bound parameters for found seed "
                << is << " in region " << ir << " failed." << std::endl;
          }
        }
      }
    }
  }
  return 0;
}
