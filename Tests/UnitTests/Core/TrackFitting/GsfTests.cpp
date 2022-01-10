// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/MultiEigenStepperLoop.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GaussianSumFitter.hpp"
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

using Stepper = Acts::MultiEigenStepperLoop<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;

using KalmanUpdater = Acts::GainMatrixUpdater;
using KalmanSmoother = Acts::GainMatrixSmoother;
using Gsf = Acts::KalmanFitter<Propagator>;

KalmanUpdater kfUpdater;
KalmanSmoother kfSmoother;

KalmanFitterExtensions getExtensions() {
  KalmanFitterExtensions extensions;
  extensions.calibrator.connect<&testSourceLinkCalibrator>();
  extensions.updater.connect<&KalmanUpdater::operator()>(&kfUpdater);
  extensions.smoother.connect<&KalmanSmoother::operator()>(&kfSmoother);
  return extensions;
}

// Instatiate the tester
FitterTester tester;

// reconstruction propagator and fitter
const auto logger = getDefaultLogger("GSF", Logging::VERBOSE);
const auto gsfZeroPropagator =
    makeConstantFieldPropagator<Stepper>(tester.geometry, 0_T);
const auto gsfZero = GaussianSumFitter(std::move(gsfZeroPropagator));

std::default_random_engine rng(42);

auto makeDefaultGsfOptions() {
  return GsfOptions{tester.geoCtx,          tester.magCtx,
                    tester.calCtx,          getExtensions(),
                    LoggerWrapper{*logger}, PropagatorPlainOptions()};
}

// A Helper type to allow us to put the MultiComponentBoundTrackParameters into
// the function so that it can also be used as SingleBoundTrackParameters for
// the MeasurementsCreator
template <typename charge_t>
struct MultiCmpsParsInterface : public SingleBoundTrackParameters<charge_t> {
  using MultiPars = MultiComponentBoundTrackParameters<charge_t>;

  MultiCmpsParsInterface(const MultiPars &p)
      : SingleBoundTrackParameters<charge_t>(
            p.referenceSurface().getSharedPtr(), p.parameters(),
            p.covariance()),
        pars(p) {}

  MultiPars pars;

  operator MultiComponentBoundTrackParameters<charge_t>() { return pars; }
};

auto makeParameters() {
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
  Acts::CurvilinearTrackParameters cp(mPos4, 0_degree, 90_degree, 1_GeV, 1_e,
                                      cov);

  // Construct bound multi component parameters from curvilinear ones
  Acts::BoundVector deltaLOC0 = Acts::BoundVector::Zero();
  deltaLOC0[eBoundLoc0] = 0.5_mm;

  Acts::BoundVector deltaLOC1 = Acts::BoundVector::Zero();
  deltaLOC1[eBoundLoc1] = 0.5_mm;

  Acts::BoundVector deltaQOP = Acts::BoundVector::Zero();
  deltaQOP[eBoundQOverP] = 0.01_GeV;

  std::vector<std::tuple<double, BoundVector, BoundSymMatrix>> cmps = {
      {1.0, cp.parameters(), cov}/*,
      {0.2, cp.parameters() + deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() + deltaLOC0 - deltaLOC1 - deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 + deltaLOC1 + deltaQOP, cov},
      {0.2, cp.parameters() - deltaLOC0 - deltaLOC1 - deltaQOP, cov}*/};

  return MultiCmpsParsInterface<SinglyCharged>(
      Acts::MultiComponentBoundTrackParameters<SinglyCharged>(
          cp.referenceSurface().getSharedPtr(), cmps));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingKalmanFitter)

// BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
//   auto multi_pars = makeParameters();
//   auto options = makeDefaultGsfOptions();
// 
//   tester.test_ZeroFieldNoSurfaceForward(gsfZero, options, multi_pars, rng, true,
//                                         true);
// }


BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceForward) {
  auto multi_pars = makeParameters();
  auto options = makeDefaultGsfOptions();

  tester.test_ZeroFieldWithSurfaceForward(gsfZero, options, multi_pars, rng, true,
                                        true);
}

// BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceBackward) {
//   auto multi_pars = makeParameters();
//   auto options = makeDefaultGsfOptions();
// 
//   tester.test_ZeroFieldWithSurfaceBackward(gsfZero, options, multi_pars, rng, true,
//                                         true);
// }
/*
BOOST_AUTO_TEST_CASE(ZeroFieldWithSurfaceAtExit) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithSurfaceAtExit(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldShuffled) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldShuffled(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithHole) {
  auto options = makeDefaultGsfOptions();

  test_ZeroFieldWithHole(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithOutliers) {
  // fitter options w/o target surface. outlier distance is set to be below the
  // default outlier distance in the `MeasurementsCreator`
  auto extensions = getExtensions();
  TestOutlierFinder tof{5_mm};
  extensions.outlierFinder.connect<&TestOutlierFinder::operator()>(&tof);

  KalmanFitterOptions options(geoCtx, magCtx, calCtx, extensions,
                                LoggerWrapper{*kfLogger},
                                PropagatorPlainOptions());

  test_ZeroFieldWithOutliers(kfZero, options, rng);
}

BOOST_AUTO_TEST_CASE(ZeroFieldWithReverseFiltering) {
  auto test = [](double threshold, bool reverse, bool expected_reversed,
                 bool expected_smoothed) {
    // Reverse filtering threshold set at 0.5 GeV
    auto extensions = getExtensions();
    TestReverseFilteringLogic trfl{threshold};
    extensions.reverseFilteringLogic
        .connect<&TestReverseFilteringLogic::operator()>(&trfl);

    KalmanFitterOptions options(geoCtx, magCtx, calCtx, extensions,
                                  LoggerWrapper{*kfLogger},
                                  PropagatorPlainOptions());
    options.reversedFiltering = reverse;
    test_ZeroFieldWithReverseFiltering(kfZero, options, rng,
                                       expected_reversed, expected_smoothed);
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
  auto options = makeDefaultGsfOptions();

  test_GlobalCovariance(kfZero, options, rng);
}
*/

BOOST_AUTO_TEST_SUITE_END()
