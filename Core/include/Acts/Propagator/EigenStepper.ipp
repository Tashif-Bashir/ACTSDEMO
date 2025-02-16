// This file is part of the Acts project.
//
// Copyright (C) 2019-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/CovarianceEngine.hpp"
#include "Acts/Utilities/QuickMath.hpp"

#include <limits>

template <typename E, typename A>
Acts::EigenStepper<E, A>::EigenStepper(
    std::shared_ptr<const MagneticFieldProvider> bField,
    double /*overstepLimit*/)
    : m_bField(std::move(bField)) {}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::makeState(
    std::reference_wrapper<const GeometryContext> gctx,
    std::reference_wrapper<const MagneticFieldContext> mctx,
    const BoundTrackParameters& par, double ssize) const -> State {
  return State{gctx, m_bField->makeCache(mctx), par, ssize};
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::resetState(State& state,
                                          const BoundVector& boundParams,
                                          const BoundSquareMatrix& cov,
                                          const Surface& surface,
                                          const double stepSize) const {
  FreeVector freeParams =
      transformBoundToFreeParameters(surface, state.geoContext, boundParams);

  // Update the stepping state
  state.pars = freeParams;
  state.cov = cov;
  state.stepSize = ConstrainedStep(stepSize);
  state.pathAccumulated = 0.;

  // Reinitialize the stepping jacobian
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
  state.jacobian = BoundMatrix::Identity();
  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const
    -> Result<BoundState> {
  return detail::boundState(
      state.geoContext, surface, state.cov, state.jacobian, state.jacTransport,
      state.derivative, state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated,
      freeToBoundCorrection);
}

template <typename E, typename A>
template <typename propagator_state_t, typename navigator_t>
bool Acts::EigenStepper<E, A>::prepareCurvilinearState(
    propagator_state_t& prop_state, const navigator_t& navigator) const {
  // test whether the accumulated path has still its initial value.
  if (prop_state.stepping.pathAccumulated == 0.) {
    // if no step was executed the path length derivates have not been
    // computed but are needed to compute the curvilinear covariance. The
    // derivates are given by k1 for a zero step width.
    if (prop_state.stepping.extension.validExtensionForStep(prop_state, *this,
                                                            navigator)) {
      // First Runge-Kutta point (at current position)
      auto& sd = prop_state.stepping.stepData;
      auto pos = position(prop_state.stepping);
      auto fieldRes = getField(prop_state.stepping, pos);
      if (fieldRes.ok()) {
        sd.B_first = *fieldRes;
        if (prop_state.stepping.extension.k1(prop_state, *this, navigator,
                                             sd.k1, sd.B_first, sd.kQoP)) {
          // dr/ds :
          prop_state.stepping.derivative.template head<3>() =
              prop_state.stepping.pars.template segment<3>(eFreeDir0);
          // d (dr/ds) / ds :
          prop_state.stepping.derivative.template segment<3>(4) = sd.k1;
          // to set dt/ds :
          prop_state.stepping.extension.finalize(
              prop_state, *this, navigator,
              prop_state.stepping.pathAccumulated);
          return true;
        }
      }
    }
    return false;
  }
  return true;
}

template <typename E, typename A>
auto Acts::EigenStepper<E, A>::curvilinearState(State& state,
                                                bool transportCov) const
    -> CurvilinearState {
  return detail::curvilinearState(
      state.cov, state.jacobian, state.jacTransport, state.derivative,
      state.jacToGlobal, state.pars, state.particleHypothesis,
      state.covTransport && transportCov, state.pathAccumulated);
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state,
                                      const FreeVector& freeParams,
                                      const BoundVector& /*boundParams*/,
                                      const Covariance& covariance,
                                      const Surface& surface) const {
  state.pars = freeParams;
  state.cov = covariance;
  state.jacToGlobal = surface.boundToFreeJacobian(
      state.geoContext, freeParams.template segment<3>(eFreePos0),
      freeParams.template segment<3>(eFreeDir0));
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::update(State& state, const Vector3& uposition,
                                      const Vector3& udirection, double qOverP,
                                      double time) const {
  state.pars.template segment<3>(eFreePos0) = uposition;
  state.pars.template segment<3>(eFreeDir0) = udirection;
  state.pars[eFreeTime] = time;
  state.pars[eFreeQOverP] = qOverP;
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToCurvilinear(
    State& state) const {
  detail::transportCovarianceToCurvilinear(state.cov, state.jacobian,
                                           state.jacTransport, state.derivative,
                                           state.jacToGlobal, direction(state));
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::transportCovarianceToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const {
  detail::transportCovarianceToBound(state.geoContext.get(), surface, state.cov,
                                     state.jacobian, state.jacTransport,
                                     state.derivative, state.jacToGlobal,
                                     state.pars, freeToBoundCorrection);
}

template <typename E, typename A>
template <typename propagator_state_t, typename navigator_t>
Acts::Result<double> Acts::EigenStepper<E, A>::step(
    propagator_state_t& state, const navigator_t& navigator) const {
  // Runge-Kutta integrator state
  auto& sd = state.stepping.stepData;

  double errorEstimate = 0.;
  double h2 = 0, half_h = 0;

  auto pos = position(state.stepping);
  auto dir = direction(state.stepping);

  // First Runge-Kutta point (at current position)
  auto fieldRes = getField(state.stepping, pos);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  sd.B_first = *fieldRes;
  if (!state.stepping.extension.validExtensionForStep(state, *this,
                                                      navigator) ||
      !state.stepping.extension.k1(state, *this, navigator, sd.k1, sd.B_first,
                                   sd.kQoP)) {
    return 0.;
  }

  const auto calcStepSizeScaling = [&](const double errorEstimate_) -> double {
    // For details about these values see ATL-SOFT-PUB-2009-001 for details
    constexpr double lower = 0.25;
    constexpr double upper = 4.0;
    // This is given by the order of the Runge-Kutta method
    constexpr double exponent = 0.25;

    // Whether to use fast power function if available
    constexpr bool tryUseFastPow{false};

    double x = state.options.stepTolerance / errorEstimate_;

    if constexpr (exponent == 0.25 && !tryUseFastPow) {
      // This is 3x faster than std::pow
      x = std::sqrt(std::sqrt(x));
    } else if constexpr (std::numeric_limits<double>::is_iec559 &&
                         tryUseFastPow) {
      x = fastPow(x, exponent);
    } else {
      x = std::pow(x, exponent);
    }

    return std::clamp(x, lower, upper);
  };

  // The following functor starts to perform a Runge-Kutta step of a certain
  // size, going up to the point where it can return an estimate of the local
  // integration error. The results are stated in the local variables above,
  // allowing integration to continue once the error is deemed satisfactory
  const auto tryRungeKuttaStep = [&](const double h) -> Result<bool> {
    // helpers because bool and std::error_code are ambiguous
    constexpr auto success = &Result<bool>::success;
    constexpr auto failure = &Result<bool>::failure;

    // State the square and half of the step size
    h2 = h * h;
    half_h = h * 0.5;

    // Second Runge-Kutta point
    const Vector3 pos1 = pos + half_h * dir + h2 * 0.125 * sd.k1;
    auto field = getField(state.stepping, pos1);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_middle = *field;

    if (!state.stepping.extension.k2(state, *this, navigator, sd.k2,
                                     sd.B_middle, sd.kQoP, half_h, sd.k1)) {
      return success(false);
    }

    // Third Runge-Kutta point
    if (!state.stepping.extension.k3(state, *this, navigator, sd.k3,
                                     sd.B_middle, sd.kQoP, half_h, sd.k2)) {
      return success(false);
    }

    // Last Runge-Kutta point
    const Vector3 pos2 = pos + h * dir + h2 * 0.5 * sd.k3;
    field = getField(state.stepping, pos2);
    if (!field.ok()) {
      return failure(field.error());
    }
    sd.B_last = *field;
    if (!state.stepping.extension.k4(state, *this, navigator, sd.k4, sd.B_last,
                                     sd.kQoP, h, sd.k3)) {
      return success(false);
    }

    // Compute and check the local integration error estimate
    errorEstimate =
        h2 * ((sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>() +
              std::abs(sd.kQoP[0] - sd.kQoP[1] - sd.kQoP[2] + sd.kQoP[3]));
    // Protect against division by zero
    errorEstimate = std::max(1e-20, errorEstimate);

    return success(errorEstimate <= state.options.stepTolerance);
  };

  const double initialH =
      state.stepping.stepSize.value() * state.options.direction;
  double h = initialH;
  std::size_t nStepTrials = 0;
  // Select and adjust the appropriate Runge-Kutta step size as given
  // ATL-SOFT-PUB-2009-001
  while (true) {
    nStepTrials++;
    auto res = tryRungeKuttaStep(h);
    if (!res.ok()) {
      return res.error();
    }
    if (!!res.value()) {
      break;
    }

    h *= calcStepSizeScaling(2 * errorEstimate);

    // If step size becomes too small the particle remains at the initial
    // place
    if (std::abs(h) < std::abs(state.options.stepSizeCutOff)) {
      // Not moving due to too low momentum needs an aborter
      return EigenStepperError::StepSizeStalled;
    }

    // If the parameter is off track too much or given stepSize is not
    // appropriate
    if (nStepTrials > state.options.maxRungeKuttaStepTrials) {
      // Too many trials, have to abort
      return EigenStepperError::StepSizeAdjustmentFailed;
    }
  }

  // When doing error propagation, update the associated Jacobian matrix
  if (state.stepping.covTransport) {
    // The step transport matrix in global coordinates
    FreeMatrix D;
    if (!state.stepping.extension.finalize(state, *this, navigator, h, D)) {
      return EigenStepperError::StepInvalid;
    }

    // See the documentation of Acts::blockedMult for a description of blocked
    // matrix multiplication. However, we can go one further. Let's assume that
    // some of these sub-matrices are zero matrices 0₈ and identity matrices
    // I₈, namely:
    //
    // D₁₁ = I₈, J₁₁ = I₈, D₂₁ = 0₈, J₂₁ = 0₈
    //
    // Which gives:
    //
    // K₁₁ = I₈  * I₈  + D₁₂ * 0₈  = I₈
    // K₁₂ = I₈  * J₁₂ + D₁₂ * J₂₂ = J₁₂ + D₁₂ * J₂₂
    // K₂₁ = 0₈  * I₈  + D₂₂ * 0₈  = 0₈
    // K₂₂ = 0₈  * J₁₂ + D₂₂ * J₂₂ = D₂₂ * J₂₂
    //
    // Furthermore, we're constructing K in place of J, and since
    // K₁₁ = I₈ = J₁₁ and K₂₁ = 0₈ = D₂₁, we don't actually need to touch those
    // sub-matrices at all!
    assert((D.topLeftCorner<4, 4>().isIdentity()));
    assert((D.bottomLeftCorner<4, 4>().isZero()));
    assert((state.stepping.jacTransport.template topLeftCorner<4, 4>()
                .isIdentity()));
    assert((state.stepping.jacTransport.template bottomLeftCorner<4, 4>()
                .isZero()));

    state.stepping.jacTransport.template topRightCorner<4, 4>() +=
        D.topRightCorner<4, 4>() *
        state.stepping.jacTransport.template bottomRightCorner<4, 4>();
    state.stepping.jacTransport.template bottomRightCorner<4, 4>() =
        (D.bottomRightCorner<4, 4>() *
         state.stepping.jacTransport.template bottomRightCorner<4, 4>())
            .eval();
  } else {
    if (!state.stepping.extension.finalize(state, *this, navigator, h)) {
      return EigenStepperError::StepInvalid;
    }
  }

  // Update the track parameters according to the equations of motion
  state.stepping.pars.template segment<3>(eFreePos0) +=
      h * dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
  state.stepping.pars.template segment<3>(eFreeDir0) +=
      h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
  (state.stepping.pars.template segment<3>(eFreeDir0)).normalize();

  if (state.stepping.covTransport) {
    state.stepping.derivative.template head<3>() =
        state.stepping.pars.template segment<3>(eFreeDir0);
    state.stepping.derivative.template segment<3>(4) = sd.k4;
  }

  state.stepping.pathAccumulated += h;
  ++state.stepping.nSteps;
  state.stepping.nStepTrials += nStepTrials;

  const double stepSizeScaling = calcStepSizeScaling(errorEstimate);
  const double nextAccuracy = std::abs(h * stepSizeScaling);
  const double previousAccuracy = std::abs(state.stepping.stepSize.accuracy());
  const double initialStepLength = std::abs(initialH);
  if (nextAccuracy < initialStepLength || nextAccuracy > previousAccuracy) {
    state.stepping.stepSize.setAccuracy(nextAccuracy);
  }

  return h;
}

template <typename E, typename A>
void Acts::EigenStepper<E, A>::setIdentityJacobian(State& state) const {
  state.jacobian = BoundMatrix::Identity();
}
