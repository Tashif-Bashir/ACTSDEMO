// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Propagator/TryAllNavigator.hpp"
#include "Acts/Propagator/TryAllOverstepNavigator.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicBVHTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>

namespace Acts {
class Layer;
struct FreeToBoundCorrection;
}  // namespace Acts

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;
using Acts::VectorHelpers::perp;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

/// This is a simple cache struct to mimic the
/// Propagator cache
struct PropagatorState {
  /// This is a simple cache struct to mimic a Stepper
  struct Stepper {
    // comply with concept
    using Jacobian = BoundMatrix;
    using Covariance = BoundSquareMatrix;
    using BoundState = std::tuple<BoundTrackParameters, Jacobian, double>;
    using CurvilinearState =
        std::tuple<CurvilinearTrackParameters, Jacobian, double>;

    /// This is a simple cache struct to mimic the
    /// Stepper cache in the propagation
    struct State {
      /// Position
      Vector4 pos4 = Vector4(0., 0., 0., 0.);

      /// Direction
      Vector3 dir = Vector3(1., 0., 0.);

      /// Momentum
      double p = 0;

      /// Charge
      double q = 0;

      /// Particle hypothesis
      ParticleHypothesis particleHypothesis = ParticleHypothesis::pion();

      // accummulated path length cache
      double pathAccumulated = 0.;

      // adaptive sep size of the runge-kutta integration
      ConstrainedStep stepSize = ConstrainedStep(100_cm);

      // Previous step size for overstep estimation (ignored here)
      double previousStepSize = 0.;

      GeometryContext geoContext = GeometryContext();
    };

    /// State resetter
    void resetState(State& /*state*/, const BoundVector& /*boundParams*/,
                    const BoundSquareMatrix& /*cov*/,
                    const Surface& /*surface*/,
                    const double /*stepSize*/) const {}

    /// Global particle position accessor
    Vector3 position(const State& state) const {
      return state.pos4.segment<3>(Acts::ePos0);
    }

    /// Time access
    double time(const State& state) const { return state.pos4[Acts::eTime]; }

    /// Momentum direction accessor
    Vector3 direction(const State& state) const { return state.dir; }

    /// QoP accessor
    double qOverP(const State& state) const {
      return (state.q == 0 ? 1 : state.q) / state.p;
    }

    /// Absolute momentum accessor
    double absoluteMomentum(const State& state) const { return state.p; }

    /// Momentum accessor
    Vector3 momentum(const State& state) const { return state.p * state.dir; }

    /// Charge access
    double charge(const State& state) const { return state.q; }

    Intersection3D::Status updateSurfaceStatus(
        State& state, const Surface& surface, std::uint8_t index,
        Direction navDir, const BoundaryCheck& bcheck,
        ActsScalar surfaceTolerance, const Logger& logger) const {
      return detail::updateSingleSurfaceStatus<Stepper>(
          *this, state, surface, index, navDir, bcheck, surfaceTolerance,
          logger);
    }

    template <typename object_intersection_t>
    void updateStepSize(State& state,
                        const object_intersection_t& oIntersection,
                        Direction /*direction*/, bool release = true) const {
      detail::updateSingleStepSize<Stepper>(state, oIntersection, release);
    }

    void setStepSize(State& state, double stepSize,
                     ConstrainedStep::Type stype = ConstrainedStep::actor,
                     bool release = true) const {
      state.previousStepSize = state.stepSize.value();
      state.stepSize.update(stepSize, stype, release);
    }

    double getStepSize(const State& state, ConstrainedStep::Type stype) const {
      return state.stepSize.value(stype);
    }

    void releaseStepSize(State& state) const {
      state.stepSize.release(ConstrainedStep::actor);
    }

    std::string outputStepSize(const State& state) const {
      return state.stepSize.toString();
    }

    Result<BoundState> boundState(
        State& state, const Surface& surface, bool /*transportCov*/,
        const FreeToBoundCorrection& /*freeToBoundCorrection*/
    ) const {
      auto bound = BoundTrackParameters::create(
          surface.getSharedPtr(), tgContext, state.pos4, state.dir,
          state.q / state.p, std::nullopt, state.particleHypothesis);
      if (!bound.ok()) {
        return bound.error();
      }
      BoundState bState{std::move(*bound), Jacobian::Identity(),
                        state.pathAccumulated};
      return bState;
    }

    CurvilinearState curvilinearState(State& state, bool /*transportCov*/
    ) const {
      CurvilinearTrackParameters parameters(state.pos4, state.dir,
                                            state.q / state.p, std::nullopt,
                                            state.particleHypothesis);
      // Create the bound state
      CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                                 state.pathAccumulated};
      return curvState;
    }

    void update(State& /*state*/, const FreeVector& /*freePars*/,
                const BoundVector& /*boundPars*/, const Covariance& /*cov*/,
                const Surface& /*surface*/) const {}

    void update(State& /*state*/, const Vector3& /*uposition*/,
                const Vector3& /*udirection*/, double /*up*/,
                double /*time*/) const {}

    void transportCovarianceToCurvilinear(State& /*state*/) const {}

    void transportCovarianceToBound(
        State& /*state*/, const Surface& /*surface*/,
        const FreeToBoundCorrection& /*freeToBoundCorrection*/) const {}

    Result<Vector3> getField(State& /*state*/, const Vector3& /*pos*/) const {
      // get the field from the cell
      return Result<Vector3>::success({0., 0., 0.});
    }
  };

  static_assert(StepperConcept<Stepper>,
                "Dummy stepper does not fulfill concept");

  /// emulate the options template
  struct Options {
    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    std::size_t debugPfxWidth = 30;
    std::size_t debugMsgWidth = 50;

    Direction direction = Direction::Forward;

    const Acts::Logger& logger = Acts::getDummyLogger();

    ActsScalar surfaceTolerance = s_onSurfaceTolerance;
  };

  /// Navigation cache: the start surface
  const Surface* startSurface = nullptr;

  /// Navigation cache: the current surface
  const Surface* currentSurface = nullptr;

  /// Navigation cache: the target surface
  const Surface* targetSurface = nullptr;
  bool targetReached = false;

  /// Give some options
  Options options;

  /// The Stepper state - internal statew of the Stepper
  Stepper::State stepping;

  /// Navigation state - internal state of the Navigator
  Navigator::State navigation;

  // The context cache for this propagation
  GeometryContext geoContext = GeometryContext();
};

template <typename stepper_state_t>
void step(stepper_state_t& sstate) {
  // update the cache position
  sstate.pos4[Acts::ePos0] += sstate.stepSize.value() * sstate.dir[Acts::eMom0];
  sstate.pos4[Acts::ePos1] += sstate.stepSize.value() * sstate.dir[Acts::eMom1];
  sstate.pos4[Acts::ePos2] += sstate.stepSize.value() * sstate.dir[Acts::eMom2];
  // create navigation parameters
  return;
}

#define testNavigatorStateVectors(state, cand, extSurf)        \
  do {                                                         \
    BOOST_CHECK_EQUAL(state.candidates.size(), cand);          \
    BOOST_CHECK_EQUAL(state.externalSurfaces.size(), extSurf); \
  } while (0)

#define testNavigatorStatePointers(state, startVol, startLay, startSurf,    \
                                   currSurf, currVol, targetVol, targetLay, \
                                   targetSurf)                              \
  do {                                                                      \
    BOOST_CHECK_EQUAL(state.startVolume, startVol);                         \
    BOOST_CHECK_EQUAL(state.startLayer, startLay);                          \
    BOOST_CHECK_EQUAL(state.startSurface, startSurf);                       \
    BOOST_CHECK_EQUAL(state.currentSurface, currSurf);                      \
    BOOST_CHECK_EQUAL(state.currentVolume, currVol);                        \
    BOOST_CHECK_EQUAL(state.targetVolume, targetVol);                       \
    BOOST_CHECK_EQUAL(state.targetLayer, targetLay);                        \
    BOOST_CHECK_EQUAL(state.targetSurface, targetSurf);                     \
  } while (0)

// the surface cache & the creation of the geometry

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

const double Bz = 2_T;
auto bField = std::make_shared<ConstantBField>(Vector3{0, 0, Bz});

// the debug boolean
bool debug = false;

BOOST_AUTO_TEST_CASE(Navigator_preStep) {
  // position and direction vector
  Vector4 position4(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  // the propagator cache
  PropagatorState state;
  state.options.debug = debug;

  // the stepper cache
  state.stepping.pos4 = position4;
  state.stepping.dir = momentum.normalized();

  // Stepper
  PropagatorState::Stepper stepper;

  //
  // (1) Test for inactivity
  //
  // Run without anything present
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    Navigator navigator{navCfg};

    navigator.postStep(state, stepper);
    testNavigatorStateVectors(state.navigation, 0u, 0u);
    testNavigatorStatePointers(state.navigation, nullptr, nullptr, nullptr,
                               nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  // Run with geometry but without resolving
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = false;
    navCfg.resolveMaterial = false;
    navCfg.resolvePassive = false;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    navigator.postStep(state, stepper);
    testNavigatorStateVectors(state.navigation, 0u, 0u);
    testNavigatorStatePointers(state.navigation, nullptr, nullptr, nullptr,
                               nullptr, nullptr, nullptr, nullptr, nullptr);
  }

  // Run with geometry and resolving but broken navigation for various reasons
  {
    Navigator::Config navCfg;
    navCfg.resolveSensitive = true;
    navCfg.resolveMaterial = true;
    navCfg.resolvePassive = true;
    navCfg.trackingGeometry = tGeometry;
    Navigator navigator{navCfg};

    state.navigation.navigationBreak = true;
    // a) Because target is reached
    state.navigation.targetReached = true;
    navigator.postStep(state, stepper);
    testNavigatorStateVectors(state.navigation, 0u, 0u);
    testNavigatorStatePointers(state.navigation, nullptr, nullptr, nullptr,
                               nullptr, nullptr, nullptr, nullptr, nullptr);

    // b) Because of no target surface
    state.navigation.targetReached = false;
    state.navigation.targetSurface = nullptr;
    navigator.postStep(state, stepper);
    testNavigatorStateVectors(state.navigation, 0u, 0u);
    testNavigatorStatePointers(state.navigation, nullptr, nullptr, nullptr,
                               nullptr, nullptr, nullptr, nullptr, nullptr);

    //
    // (2) Test the initialisation
    //
    // a) Initialise without additional information
    state.navigation = Navigator::State();
    state.stepping.pos4 << 0., 0., 0., 0.;
    const TrackingVolume* startVol = tGeometry->lowestTrackingVolume(
        state.geoContext, stepper.position(state.stepping));
    const Layer* startLay = startVol->associatedLayer(
        state.geoContext, stepper.position(state.stepping));
    navigator.initialize(state, stepper);
    testNavigatorStateVectors(state.navigation, 2u, 0u);
    testNavigatorStatePointers(state.navigation, startVol, nullptr, nullptr,
                               nullptr, startVol, nullptr, nullptr, nullptr);

    // b) Initialise having a start surface
    const Surface* startSurf = tGeometry->getBeamline();
    state.navigation = Navigator::State();
    state.navigation.startSurface = startSurf;
    navigator.initialize(state, stepper);
    testNavigatorStateVectors(state.navigation, 2u, 0u);
    testNavigatorStatePointers(state.navigation, startVol, nullptr, startSurf,
                               startSurf, startVol, nullptr, nullptr, nullptr);

    // c) Initialise having a start volume
    state.navigation = Navigator::State();
    state.navigation.startVolume = startVol;
    navigator.initialize(state, stepper);
    testNavigatorStateVectors(state.navigation, 2u, 0u);
    testNavigatorStatePointers(state.navigation, startVol, startLay, nullptr,
                               nullptr, startVol, nullptr, nullptr, nullptr);
  }
}

BOOST_AUTO_TEST_CASE(Navigator_postStep) {
  // create a navigator
  Navigator::Config navCfg;
  navCfg.trackingGeometry = tGeometry;
  navCfg.resolveSensitive = true;
  navCfg.resolveMaterial = true;
  navCfg.resolvePassive = false;
  Navigator navigator{navCfg};

  // create a navigator for the Bounding Volume Hierarchy test
  CubicBVHTrackingGeometry grid(20, 1000, 5);
  Navigator::Config bvhNavCfg;
  bvhNavCfg.trackingGeometry = grid.trackingGeometry;
  bvhNavCfg.resolveSensitive = true;
  bvhNavCfg.resolveMaterial = true;
  bvhNavCfg.resolvePassive = false;
  Navigator BVHNavigator{bvhNavCfg};

  // position and direction vector
  Vector4 position4(0., 0., 0, 0);
  Vector3 momentum(1., 1., 0);

  // the propagator cache
  PropagatorState state;
  state.options.debug = debug;

  // the stepper cache
  state.stepping.pos4 = position4;
  state.stepping.dir = momentum.normalized();

  // forward navigation ----------------------------------------------
  if (debug) {
    std::cout << "<<<<<<<<<<<<<<<<<<<<< FORWARD NAVIGATION >>>>>>>>>>>>>>>>>>"
              << std::endl;
  }

  // Stepper
  PropagatorState::Stepper stepper;

  // (1) Initialization navigation from start point
  // - this will call resolveLayers() as well
  // - and thus should call a return to the stepper
  navigator.initialize(state, stepper);
  // Check that the currentVolume is set
  BOOST_CHECK_NE(state.navigation.currentVolume, nullptr);
  // Check that the currentVolume is the startVolume
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  // Check that the currentSurface is reset to:
  BOOST_CHECK_EQUAL(state.navigation.currentSurface, nullptr);
  // No layer has been found
  BOOST_CHECK_EQUAL(state.navigation.candidates.size(), 2u);
  // ACTORS-ABORTERS-TARGET
  navigator.preStep(state, stepper);
  // A layer has been found
  BOOST_CHECK_EQUAL(state.navigation.candidates.size(), 2u);
  // The index should points to the begin
  BOOST_CHECK_EQUAL(state.navigation.candidateIndex, 0);
  // Cache the beam pipe radius
  double beamPipeR = perp(state.navigation.candidate().intersection.position());
  // step size has been updated
  CHECK_CLOSE_ABS(state.stepping.stepSize.value(), beamPipeR,
                  s_onSurfaceTolerance);
  if (debug) {
    std::cout << "<<< Test 1a >>> initialize at "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    // Clear the debug string for the next test
    state.options.debugString = "";
  }

  // Do the step towards the beam pipe
  step(state.stepping);

  // (2) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // Check that the currentVolume is the still startVolume
  BOOST_CHECK_EQUAL(state.navigation.currentVolume,
                    state.navigation.startVolume);
  // The layer number has not changed
  BOOST_CHECK_EQUAL(state.navigation.candidates.size(), 1u);
  // The index still points to the begin
  BOOST_CHECK_EQUAL(state.navigation.candidateIndex, 0);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1b >>> step to the BeamPipe at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Do the step towards the boundary
  step(state.stepping);

  // (3) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1c >>> step to the Boundary at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // positive return: do the step
  step(state.stepping);
  // (4) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1d >>> step to 1st layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on first layer
  for (std::size_t isf = 0; isf < 5; ++isf) {
    step(state.stepping);
    // (5-9) re-entering navigator:
    // POST STEP
    navigator.postStep(state, stepper);
    // ACTORS - ABORTERS - PRE STEP
    navigator.preStep(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1e-1i >>> step within 1st layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (10) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1j >>> step to 2nd layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on second layer
  for (std::size_t isf = 0; isf < 5; ++isf) {
    step(state.stepping);
    // (11-15) re-entering navigator:
    // POST STEP
    navigator.postStep(state, stepper);
    // ACTORS - ABORTERS - PRE STEP
    navigator.preStep(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1k-1o >>> step within 2nd layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (16) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1p >>> step to 3rd layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on third layer
  for (std::size_t isf = 0; isf < 3; ++isf) {
    step(state.stepping);
    // (17-19) re-entering navigator:
    // POST STEP
    navigator.postStep(state, stepper);
    // ACTORS - ABORTERS - PRE STEP
    navigator.preStep(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1q-1s >>> step within 3rd layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (20) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1t >>> step to 4th layer at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }

  // Step through the surfaces on second layer
  for (std::size_t isf = 0; isf < 3; ++isf) {
    step(state.stepping);
    // (21-23) re-entering navigator:
    // POST STEP
    navigator.postStep(state, stepper);
    // ACTORS - ABORTERS - PRE STEP
    navigator.preStep(state, stepper);

    if (debug) {
      std::cout << "<<< Test 1t-1v >>> step within 4th layer at  "
                << toString(state.stepping.pos4) << std::endl;
      std::cout << state.options.debugString << std::endl;
      state.options.debugString = "";
    }
  }

  // positive return: do the step
  step(state.stepping);
  // (24) re-entering navigator:
  // POST STEP
  navigator.postStep(state, stepper);
  // ACTORS - ABORTERS - PRE STEP
  navigator.preStep(state, stepper);

  if (debug) {
    std::cout << "<<< Test 1w >>> step to boundary at  "
              << toString(state.stepping.pos4) << std::endl;
    std::cout << state.options.debugString << std::endl;
    state.options.debugString = "";
  }
}

using SurfaceCollector = SurfaceCollector<SurfaceSelector>;

std::vector<GeometryIdentifier> collectGeoIds(
    const SurfaceCollector::result_type& surfaceHits) {
  std::vector<GeometryIdentifier> geoIds;
  for (const auto& surfaceHit : surfaceHits.collected) {
    geoIds.push_back(surfaceHit.surface->geometryId());
  }
  return geoIds;
}

/// the actual test nethod that runs the test can be used with several
/// propagator types
///
/// @tparam propagator_t is the actual propagator type
///
/// @param prop is the propagator instance
/// @param start start parameters for propagation
/// @param debugMode toggle debug mode
template <typename propagator_t>
void runSelfConsistencyTest(const propagator_t& prop,
                            const CurvilinearTrackParameters& start,
                            bool debugMode) {
  // Action list and abort list
  using ActionListType = ActionList<SurfaceCollector>;
  using AbortListType = AbortList<>;
  using Options = PropagatorOptions<ActionListType, AbortListType>;

  // forward surface test
  Options fwdOptions(tgContext, mfContext);
  fwdOptions.pathLimit = 25_cm;
  fwdOptions.maxStepSize = 1_cm;

  // get the surface collector and configure it
  auto& fwdSurfaceCollector =
      fwdOptions.actionList.template get<SurfaceCollector>();
  fwdSurfaceCollector.selector.selectSensitive = true;
  fwdSurfaceCollector.selector.selectMaterial = true;
  fwdSurfaceCollector.selector.selectPassive = true;

  if (debugMode) {
    std::cout << ">>> Forward Propagation : start." << std::endl;
  }
  auto fwdResult = prop.propagate(start, fwdOptions).value();
  auto fwdSurfaceHits =
      fwdResult.template get<SurfaceCollector::result_type>().collected;
  auto fwdSurfaces =
      collectGeoIds(fwdResult.template get<SurfaceCollector::result_type>());

  // get the forward output to the screen
  if (debugMode) {
    // check if the surfaces are free
    std::cout << ">>> Surface hits found on ..." << std::endl;
    for (const auto& fwdSteps : fwdSurfaces) {
      std::cout << "--> Surface with " << fwdSteps << std::endl;
    }
    std::cout << ">>> Forward Propagation : end." << std::endl;
  }

  // backward surface test
  Options bwdOptions(tgContext, mfContext);
  bwdOptions.pathLimit = 25_cm;
  bwdOptions.maxStepSize = 1_cm;
  bwdOptions.direction = Direction::Backward;

  // get the surface collector and configure it
  auto& bwdMSurfaceCollector =
      bwdOptions.actionList.template get<SurfaceCollector>();
  bwdMSurfaceCollector.selector.selectSensitive = true;
  bwdMSurfaceCollector.selector.selectMaterial = true;
  bwdMSurfaceCollector.selector.selectPassive = true;

  const auto& startSurface = start.referenceSurface();

  if (debugMode) {
    std::cout << ">>> Backward Propagation : start." << std::endl;
  }
  auto bwdResult =
      prop.propagate(*fwdResult.endParameters, startSurface, bwdOptions)
          .value();
  auto bwdSurfaceHits =
      bwdResult.template get<SurfaceCollector::result_type>().collected;
  auto bwdSurfaces =
      collectGeoIds(bwdResult.template get<SurfaceCollector::result_type>());

  // get the backward output to the screen
  if (debugMode) {
    // check if the surfaces are free
    std::cout << ">>> Surface hits found on ..." << std::endl;
    for (auto& bwdSteps : bwdSurfaces) {
      std::cout << "--> Surface with " << bwdSteps << std::endl;
    }
    std::cout << ">>> Backward Propagation : end." << std::endl;
  }

  // forward-backward compatibility test
  {
    // remove the undefined surface and reverse to make comparable
    auto fwdTmp = std::vector(fwdSurfaces.begin() + 1, fwdSurfaces.end());
    auto bwdTmp = std::vector(bwdSurfaces.begin() + 1, bwdSurfaces.end());
    std::reverse(bwdTmp.begin(), bwdTmp.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(bwdTmp.begin(), bwdTmp.end(), fwdTmp.begin(),
                                  fwdTmp.end());
  }

  // stepping from one surface to the next
  // now go from surface to surface and check
  Options fwdStepOptions(tgContext, mfContext);
  fwdStepOptions.maxStepSize = 1_cm;

  // get the surface collector and configure it
  auto& fwdStepSurfaceCollector =
      fwdOptions.actionList.template get<SurfaceCollector>();
  fwdStepSurfaceCollector.selector.selectSensitive = true;
  fwdStepSurfaceCollector.selector.selectMaterial = true;
  fwdStepSurfaceCollector.selector.selectPassive = true;

  std::vector<GeometryIdentifier> fwdStepSurfaces;

  // move forward step by step through the surfaces
  BoundTrackParameters sParameters = start;
  std::vector<BoundTrackParameters> stepParameters;
  for (auto& fwdSteps : fwdSurfaceHits) {
    if (debugMode) {
      std::cout << ">>> Forward step : "
                << sParameters.referenceSurface().geometryId() << " --> "
                << fwdSteps.surface->geometryId() << std::endl;
    }

    // make a forward step
    auto fwdStep =
        prop.propagate(sParameters, *fwdSteps.surface, fwdStepOptions).value();

    auto fwdStepSurfacesTmp =
        collectGeoIds(fwdStep.template get<SurfaceCollector::result_type>());
    fwdStepSurfaces.insert(fwdStepSurfaces.end(), fwdStepSurfacesTmp.begin(),
                           fwdStepSurfacesTmp.end());

    if (fwdStep.endParameters.has_value()) {
      // make sure the parameters do not run out of scope
      stepParameters.push_back(*fwdStep.endParameters);
      sParameters = stepParameters.back();
    }
  }
  // final destination surface
  const Surface& dSurface = fwdResult.endParameters->referenceSurface();
  if (debugMode) {
    std::cout << ">>> Forward step : "
              << sParameters.referenceSurface().geometryId() << " --> "
              << dSurface.geometryId() << std::endl;
  }
  auto fwdStepFinal =
      prop.propagate(sParameters, dSurface, fwdStepOptions).value();
  auto fwdStepSurfacesTmp =
      collectGeoIds(fwdStepFinal.template get<SurfaceCollector::result_type>());
  fwdStepSurfaces.insert(fwdStepSurfaces.end(), fwdStepSurfacesTmp.begin(),
                         fwdStepSurfacesTmp.end());

  // TODO forward-forward step compatibility test

  // stepping from one surface to the next : backwards
  // now go from surface to surface and check
  Options bwdStepOptions(tgContext, mfContext);
  bwdStepOptions.maxStepSize = 1_cm;
  bwdStepOptions.direction = Direction::Backward;

  // get the surface collector and configure it
  auto& bwdStepSurfaceCollector =
      bwdOptions.actionList.template get<SurfaceCollector>();
  bwdStepSurfaceCollector.selector.selectSensitive = true;
  bwdStepSurfaceCollector.selector.selectMaterial = true;
  bwdStepSurfaceCollector.selector.selectPassive = true;

  std::vector<GeometryIdentifier> bwdStepSurfaces;

  // move forward step by step through the surfaces
  sParameters = *fwdResult.endParameters;
  for (auto& bwdSteps : bwdSurfaceHits) {
    if (debugMode) {
      std::cout << ">>> Backward step : "
                << sParameters.referenceSurface().geometryId() << " --> "
                << bwdSteps.surface->geometryId() << std::endl;
    }

    // make a forward step
    auto bwdStep =
        prop.propagate(sParameters, *bwdSteps.surface, bwdStepOptions).value();

    auto bwdStepSurfacesTmp =
        collectGeoIds(bwdStep.template get<SurfaceCollector::result_type>());
    bwdStepSurfaces.insert(bwdStepSurfaces.end(), bwdStepSurfacesTmp.begin(),
                           bwdStepSurfacesTmp.end());

    if (bwdStep.endParameters.has_value()) {
      // make sure the parameters do not run out of scope
      stepParameters.push_back(*bwdStep.endParameters);
      sParameters = stepParameters.back();
    }
  }
  // final destination surface
  const Surface& dbSurface = start.referenceSurface();
  if (debugMode) {
    std::cout << ">>> Backward step : "
              << sParameters.referenceSurface().geometryId() << " --> "
              << dSurface.geometryId() << std::endl;
  }
  auto bwdStepFinal =
      prop.propagate(sParameters, dbSurface, bwdStepOptions).value();
  auto bwdStepSurfacesTmp =
      collectGeoIds(bwdStepFinal.template get<SurfaceCollector::result_type>());
  bwdStepSurfaces.insert(bwdStepSurfaces.end(), bwdStepSurfacesTmp.begin(),
                         bwdStepSurfacesTmp.end());

  // TODO backward-backward step compatibility test

  std::reverse(bwdStepSurfaces.begin(), bwdStepSurfaces.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(bwdStepSurfaces.begin(), bwdStepSurfaces.end(),
                                fwdStepSurfaces.begin(), fwdStepSurfaces.end());
}

/// the actual test nethod that runs the test can be used with several
/// propagator types
///
/// @tparam propagator_probe_t is the probe propagator type
/// @tparam propagator_ref_t is the reference propagator type
///
/// @param propProbe is the probe propagator instance
/// @param propRef is the reference propagator instance
/// @param start start parameters for propagation
/// @param debugMode toggle debug mode
template <typename propagator_probe_t, typename propagator_ref_t>
void runConsistencyTest(const propagator_probe_t& propProbe,
                        const propagator_ref_t& propRef,
                        const CurvilinearTrackParameters& start,
                        bool debugMode) {
  // Action list and abort list
  using ActionListType = ActionList<SurfaceCollector>;
  using AbortListType = AbortList<>;
  using Options = PropagatorOptions<ActionListType, AbortListType>;

  auto run = [&](const auto& prop) {
    // forward surface test
    Options fwdOptions(tgContext, mfContext);
    fwdOptions.maxStepSize = 1_cm;

    // get the surface collector and configure it
    auto& fwdSurfaceCollector =
        fwdOptions.actionList.template get<SurfaceCollector>();
    fwdSurfaceCollector.selector.selectSensitive = true;
    fwdSurfaceCollector.selector.selectMaterial = true;
    fwdSurfaceCollector.selector.selectPassive = true;

    auto fwdResult = prop.propagate(start, fwdOptions).value();
    auto fwdSurfaces =
        collectGeoIds(fwdResult.template get<SurfaceCollector::result_type>());

    // get the forward output to the screen
    if (debugMode) {
      // check if the surfaces are free
      std::cout << ">>> Surface hits found on ..." << std::endl;
      for (const auto& fwdSteps : fwdSurfaces) {
        std::cout << "--> Surface with " << fwdSteps << std::endl;
      }
    }

    return fwdSurfaces;
  };

  if (debugMode) {
    std::cout << ">>> Probe Propagation : start." << std::endl;
  }
  const auto& probeSurfaces = run(propProbe);
  if (debugMode) {
    std::cout << ">>> Probe Propagation : end." << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Reference Propagation : start." << std::endl;
  }
  const auto& refSurfaces = run(propRef);
  if (debugMode) {
    std::cout << ">>> Reference Propagation : end." << std::endl;
  }

  // probe-ref compatibility test
  BOOST_CHECK_EQUAL_COLLECTIONS(probeSurfaces.begin(), probeSurfaces.end(),
                                refSurfaces.begin(), refSurfaces.end());
}

int ntests = 80;
int skip = 0;
bool debugMode = false;

using EigenStepper = Acts::EigenStepper<>;
using EigenPropagator = Propagator<EigenStepper, Navigator>;
using StraightLinePropagator = Propagator<StraightLineStepper, Navigator>;
using Reference1EigenPropagator = Propagator<EigenStepper, TryAllNavigator>;
using Reference1StraightLinePropagator =
    Propagator<StraightLineStepper, TryAllNavigator>;
using Reference2EigenPropagator =
    Propagator<EigenStepper, TryAllOverstepNavigator>;
using Reference2StraightLinePropagator =
    Propagator<StraightLineStepper, TryAllOverstepNavigator>;

EigenStepper estepper(bField);
StraightLineStepper slstepper;

EigenPropagator epropagator(estepper,
                            Navigator({tGeometry, true, true, true,
                                       BoundaryCheck(false)},
                                      getDefaultLogger("nav", Logging::INFO)),
                            getDefaultLogger("prop", Logging::INFO));
StraightLinePropagator slpropagator(
    slstepper,
    Navigator({tGeometry, true, true, true, BoundaryCheck(false)},
              getDefaultLogger("nav", Logging::INFO)),
    getDefaultLogger("prop", Logging::INFO));

Reference1EigenPropagator refepropagator1(
    estepper,
    TryAllNavigator({tGeometry, true, true, true,
                     BoundaryCheck(true, true, 1, 1)},
                    getDefaultLogger("nav", Logging::INFO)),
    getDefaultLogger("prop", Logging::INFO));
Reference1StraightLinePropagator refslpropagator1(
    slstepper,
    TryAllNavigator({tGeometry, true, true, true,
                     BoundaryCheck(true, true, 1, 1)},
                    getDefaultLogger("nav", Logging::INFO)),
    getDefaultLogger("prop", Logging::INFO));

Reference2EigenPropagator refepropagator2(
    estepper,
    TryAllOverstepNavigator({tGeometry, true, true, true,
                             BoundaryCheck(true, true, 1, 1)},
                            getDefaultLogger("nav", Logging::INFO)),
    getDefaultLogger("prop", Logging::INFO));
Reference2StraightLinePropagator refslpropagator2(
    slstepper,
    TryAllOverstepNavigator({tGeometry, true, true, true,
                             BoundaryCheck(true, true, 1, 1)},
                            getDefaultLogger("nav", Logging::INFO)),
    getDefaultLogger("prop", Logging::INFO));

BOOST_DATA_TEST_CASE(
    Navigator_random,
    bdata::random((bdata::engine = std::mt19937(), bdata::seed = 20,
                   bdata::distribution = std::uniform_real_distribution<double>(
                       0.5_GeV, 10_GeV))) ^
        bdata::random((bdata::engine = std::mt19937(), bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<double>(-M_PI,
                                                                  M_PI))) ^
        bdata::random(
            (bdata::engine = std::mt19937(), bdata::seed = 22,
             bdata::distribution =
                 std::uniform_real_distribution<double>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::engine = std::mt19937(), bdata::seed = 23,
             bdata::distribution = std::uniform_int_distribution<int>(0, 1))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, index) {
  if (index < skip) {
    return;
  }

  double p = pT / sin(theta);
  double q = -1 + 2 * charge;
  CurvilinearTrackParameters start(Vector4(0, 0, 0, 0), phi, theta, q / p,
                                   std::nullopt, ParticleHypothesis::pion());

  if (debugMode) {
    std::cout << ">>> Run navigation tests with pT = " << pT
              << "; phi = " << phi << "; theta = " << theta
              << "; charge = " << charge << "; index = " << index << ";"
              << std::endl;
  }

  if (debugMode) {
    std::cout << ">>> Test self consistency epropagator" << std::endl;
  }
  runSelfConsistencyTest(epropagator, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test self consistency slpropagator" << std::endl;
  }
  runSelfConsistencyTest(slpropagator, start, debugMode);

  if (debugMode) {
    std::cout << ">>> Test reference 1 consistency epropagator" << std::endl;
  }
  runConsistencyTest(epropagator, refepropagator1, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test reference 1 consistency slpropagator" << std::endl;
  }
  runConsistencyTest(slpropagator, refslpropagator1, start, debugMode);

  if (debugMode) {
    std::cout << ">>> Test reference 2 consistency epropagator" << std::endl;
  }
  runConsistencyTest(epropagator, refepropagator2, start, debugMode);
  if (debugMode) {
    std::cout << ">>> Test reference 2 consistency slpropagator" << std::endl;
  }
  runConsistencyTest(slpropagator, refslpropagator2, start, debugMode);
}

}  // namespace Test
}  // namespace Acts
