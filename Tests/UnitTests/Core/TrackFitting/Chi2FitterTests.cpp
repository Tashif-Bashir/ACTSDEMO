// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/TrackFitting/Chi2Fitter.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/TrackFitting/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <algorithm>
#include <memory>
#include <random>

namespace {

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;

// using KalmanUpdater = Acts::GainMatrixUpdater;
// using KalmanSmoother = Acts::GainMatrixSmoother;
// using KalmanFitter =  Acts::KalmanFitter<ConstantFieldPropagator,
// KalmanUpdater, KalmanSmoother>;
using Chi2Fitter = Acts::Chi2Fitter<ConstantFieldPropagator>;

Chi2FitterExtensions getExtensions() {
  Chi2FitterExtensions extensions;
  extensions.calibrator.connect<&testSourceLinkCalibrator>();
  return extensions;
}

/// Find outliers using plain distance for testing purposes.
///
/// In a real setup, the outlier classification can be much more involved, e.g.
/// by computing the weighted distance/ local chi2 value. Here, the purpose is
/// to test that the basic principle works using simplified, synthetic data.
/// Thus, the simplest possible implementation should do.
struct TestOutlierFinder {
  double distanceMax = std::numeric_limits<double>::max();

  /// Classify a measurement as a valid one or an outlier.
  ///
  /// @tparam track_state_t Type of the track state
  /// @param state The track state to classify
  /// @retval False if the measurement is not an outlier
  /// @retval True if the measurement is an outlier
  bool operator()(MultiTrajectory::ConstTrackStateProxy state) const {
    // can't determine an outlier w/o a measurement or predicted parameters
    if (not state.hasCalibrated() or not state.hasPredicted()) {
      return false;
    }
    auto residuals = state.calibrated() - state.projector() * state.predicted();
    auto distance = residuals.norm();
    return (distanceMax <= distance);
  }
};

// Construct a straight-line propagator.
StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Acts::StraightLineStepper stepper;
  return StraightPropagator(std::move(stepper), std::move(navigator));
}

// Construct a propagator using a constant magnetic field along z.
ConstantFieldPropagator makeConstantFieldPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo, double bz) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, bz));
  ConstantFieldStepper stepper(std::move(field));
  return ConstantFieldPropagator(std::move(stepper), std::move(navigator));
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

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();
// expected number of measurements for the given detector
constexpr size_t nMeasurements = 6u;

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
const MeasurementResolution resStrip1 = {MeasurementType::eLoc1, {150_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(2), resPixel},
    {GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
    {GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
};

// simulation propagator
const auto simPropagator = makeStraightPropagator(geometry);

// reconstruction propagator and fitter
const auto chi2Logger =
    getDefaultLogger("Chi2Fitter", Logging::VERBOSE);  // TODO: INFO for KF
const auto chi2ZeroPropagator = makeConstantFieldPropagator(geometry, 0_T);
const auto chi2Zero = Chi2Fitter(chi2ZeroPropagator);

std::default_random_engine rng(42);

}  // namespace

BOOST_AUTO_TEST_SUITE(TrackFittingChi2Fitter)

BOOST_AUTO_TEST_CASE(ZeroFieldNoSurfaceForward) {
  auto start = makeParameters();
  auto measurements = createMeasurements(simPropagator, geoCtx, magCtx, start,
                                         resolutions, rng);
  const auto& sourceLinks = measurements.sourceLinks;
  BOOST_REQUIRE_EQUAL(sourceLinks.size(), nMeasurements);

  Chi2FitterOptions chi2Options(geoCtx, magCtx, calCtx, getExtensions(),
                                LoggerWrapper{*chi2Logger},
                                PropagatorPlainOptions());

  // chi2Options.nUpdates = 2; //  χ² = 17.9695 -> 11.0035 -> 11.0035 ...

  // BOOST_TEST_INFO("Test Case ZeroFieldNoSurfaceForward: running .fit()...");

  auto res =
      chi2Zero.fit(sourceLinks.begin(), sourceLinks.end(), start, chi2Options);
  BOOST_REQUIRE(res.ok());

  const auto& val = res.value();
  //   BOOST_CHECK_NE(val.trackTip, SIZE_MAX);
  //   BOOST_CHECK(not val.fittedParameters);
  //   BOOST_CHECK_EQUAL(val.measurementStates, sourceLinks.size());
  // check the output status flags
  //   BOOST_CHECK(val.smoothed);
  //   BOOST_CHECK(not val.reversed);
  //   BOOST_CHECK(not val.reset);
  BOOST_CHECK(val.finished);
  //   BOOST_CHECK_EQUAL(val.missedActiveSurfaces.size(), 0u);
}

// TODO: add more test cases, for holes, outliers, ...

BOOST_AUTO_TEST_SUITE_END()
