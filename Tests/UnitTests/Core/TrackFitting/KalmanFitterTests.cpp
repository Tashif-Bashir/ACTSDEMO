// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"

#include <algorithm>
#include <memory>
#include <random>

#include "FitterTestsCommon.hpp"

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using KalmanFitter = Acts::KalmanFitter<ConstantFieldPropagator>;

KalmanUpdater kfUpdater;
KalmanSmoother kfSmoother;

KalmanFitterExtensions getExtensions() {
  KalmanFitterExtensions extensions;
  extensions.calibrator.connect<&testSourceLinkCalibrator>();
  extensions.updater.connect<&KalmanUpdater::operator()>(&kfUpdater);
  extensions.smoother.connect<&KalmanSmoother::operator()>(&kfSmoother);
  return extensions;
}

// Construct initial track parameters.
Acts::CurvilinearTrackParameters makeParameters() {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  Acts::BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // define a track in the transverse plane along x
  Acts::Vector4 mPos4(-3_m, 0., 0., 42_ns);
  return Acts::CurvilinearTrackParameters(mPos4, 0_degree, 90_degree, 1_GeV,
                                          1_e, cov);
}

// Instantiate the tester
const FitterTester tester;

// reconstruction propagator and fitter
const auto kfLogger = getDefaultLogger("KalmanFilter", Logging::INFO);
const auto kfZeroPropagator =
    makeConstantFieldPropagator<ConstantFieldStepper>(tester.geometry, 0_T);
const auto kfZero = KalmanFitter(kfZeroPropagator);

std::default_random_engine rng(42);

auto makeDefaultKalmanFitterOptions() {
  return KalmanFitterOptions(tester.geoCtx, tester.magCtx, tester.calCtx,
                             getExtensions(), LoggerWrapper{*kfLogger},
                             PropagatorPlainOptions());
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingKalmanFitter)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldNoSurfaceForward(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceForward(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceBackward(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithSurfaceAtExit(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldShuffled(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_ZeroFieldWithHole(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  auto start = makeParameters();

  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  auto extensions = getExtensions();
  TestOutlierFinder tof{5_mm};
  extensions.outlierFinder.connect<&TestOutlierFinder::operator()>(&tof);

  KalmanFitterOptions kfOptions(tester.geoCtx, tester.magCtx, tester.calCtx,
                                extensions, LoggerWrapper{*kfLogger},
                                PropagatorPlainOptions());

  tester.test_ZeroFieldWithOutliers(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithReverseFiltering) {
  auto start = makeParameters();

  auto test = [&](double threshold, bool reverse, bool expected_reversed,
                  bool expected_smoothed) {
    // Reverse filtering threshold set at 0.5 GeV
    auto extensions = getExtensions();
    TestReverseFilteringLogic trfl{threshold};
    extensions.reverseFilteringLogic
        .connect<&TestReverseFilteringLogic::operator()>(&trfl);

    KalmanFitterOptions kfOptions(tester.geoCtx, tester.magCtx, tester.calCtx,
                                  extensions, LoggerWrapper{*kfLogger},
                                  PropagatorPlainOptions());
    kfOptions.reversedFiltering = reverse;
    tester.test_ZeroFieldWithReverseFiltering(
        kfZero, kfOptions, start, rng, expected_reversed, expected_smoothed);
  };

  // Track of 1 GeV with a threshold set at 0.1 GeV, reversed filtering should
  // not be used
  test(0.1_GeV, false, false, true);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(10._GeV, false, true, false);

  // Track of 1 GeV with a threshold set at 10 GeV, reversed filtering should
  // be used
  test(0.1_GeV, true, true, false);
}

// TODO this is not really Kalman fitter specific. is probably better tested
// with a synthetic trajectory.
BOOST_AUTO_TEST_CASE(GlobalCovariance) {
  auto start = makeParameters();
  auto kfOptions = makeDefaultKalmanFitterOptions();

  tester.test_GlobalCovariance(kfZero, kfOptions, start, rng);
}

BOOST_AUTO_TEST_SUITE_END()
