// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <utility>

namespace Acts {

enum class TrackExtrapolationStrategy {
  /// Use the first track state to reach target surface
  first,
  /// Use the last track state to reach target surface
  last,
  /// Use the first or last track state to reach target surface depending on the
  /// distance
  firstOrLast,
};

enum class TrackExtrapolationError {
  CompatibleTrackStateNotFound = 1,
  ReferenceSurfaceUnreachable = 2,
};

std::error_code make_error_code(TrackExtrapolationError e);

/// @brief Smooth a track using the gain matrix smoother
///
/// @tparam track_container_t The track container type
///
/// @param geoContext The geometry context
/// @param trackContainer The track container
/// @param track The track to smooth
/// @param logger The logger
///
/// @return The result of the smoothing
template <typename track_container_t>
Result<void> smoothTrack(
    const GeometryContext &geoContext, track_container_t &trackContainer,
    typename track_container_t::TrackProxy &track,
    const Logger &logger = *getDefaultLogger("TrackSmoother", Logging::INFO)) {
  Acts::GainMatrixSmoother smoother;

  auto smoothingResult =
      smoother(geoContext, trackContainer.trackStateContainer(),
               track.tipIndex(), logger);

  if (!smoothingResult.ok()) {
    ACTS_ERROR("Smoothing track " << track.index() << " failed with error "
                                  << smoothingResult.error());
    return smoothingResult.error();
  }

  return Result<void>::success();
}

/// @brief Smooth tracks using the gain matrix smoother
///
/// @tparam track_container_t The track container type
///
/// @param geoContext The geometry context
/// @param trackContainer The track container
/// @param logger The logger
///
/// @return The result of the smoothing
template <typename track_container_t>
Result<void> smoothTracks(
    const GeometryContext &geoContext, const track_container_t &trackContainer,
    const Logger &logger = *getDefaultLogger("TrackSmoother", Logging::INFO)) {
  Result<void> result = Result<void>::success();

  for (const auto &track : trackContainer) {
    auto smoothingResult =
        smoothTrack(geoContext, trackContainer, track, logger);

    if (!smoothingResult.ok()) {
      result = smoothingResult.error();
      continue;
    }
  }

  return result;
}

/// @brief Find a track state for extrapolation
///
/// @tparam track_proxy_t The track proxy type
///
/// @param geoContext The geometry context
/// @param track The track
/// @param referenceSurface The reference surface
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the search
template <typename track_proxy_t>
Result<std::pair<typename track_proxy_t::ConstTrackStateProxy, double>>
findTrackStateForExtrapolation(
    const GeometryContext &geoContext, track_proxy_t &track,
    const Surface &referenceSurface, TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  using TrackStateProxy = typename track_proxy_t::ConstTrackStateProxy;

  auto findFirst = [&]() -> Result<TrackStateProxy> {
    // TODO specialize if track is forward linked

    auto result = Result<TrackStateProxy>::failure(
        TrackExtrapolationError::CompatibleTrackStateNotFound);

    for (const auto &trackState : track.trackStatesReversed()) {
      bool isMeasurement =
          trackState.typeFlags().test(TrackStateFlag::MeasurementFlag);
      bool isOutlier = trackState.typeFlags().test(TrackStateFlag::OutlierFlag);

      if (isMeasurement && !isOutlier) {
        result = trackState;
      }
    }

    return result;
  };

  auto findLast = [&]() -> Result<TrackStateProxy> {
    for (const auto &trackState : track.trackStatesReversed()) {
      bool isMeasurement =
          trackState.typeFlags().test(TrackStateFlag::MeasurementFlag);
      bool isOutlier = trackState.typeFlags().test(TrackStateFlag::OutlierFlag);

      if (isMeasurement && !isOutlier) {
        return trackState;
      }
    }

    return Result<TrackStateProxy>::failure(
        TrackExtrapolationError::CompatibleTrackStateNotFound);
  };

  auto intersect = [&](const TrackStateProxy &state) -> SurfaceIntersection {
    auto freeVector = MultiTrajectoryHelpers::freeSmoothed(geoContext, state);

    return referenceSurface
        .intersect(geoContext, freeVector.template segment<3>(eFreePos0),
                   freeVector.template segment<3>(eFreeDir0),
                   BoundaryCheck(true), s_onSurfaceTolerance)
        .closest();
  };

  if (strategy == TrackExtrapolationStrategy::first) {
    ACTS_VERBOSE("looking for first track state");

    auto first = findFirst();
    if (!first.ok()) {
      ACTS_ERROR("no first track state found");
      return first.error();
    }

    SurfaceIntersection intersection = intersect(*first);
    if (!intersection) {
      ACTS_ERROR("no intersection found");
      return Result<std::pair<TrackStateProxy, double>>::failure(
          TrackExtrapolationError::ReferenceSurfaceUnreachable);
    }

    ACTS_VERBOSE("found intersection at " << intersection.pathLength());
    return std::make_pair(*first, intersection.pathLength());
  }

  if (strategy == TrackExtrapolationStrategy::last) {
    ACTS_VERBOSE("looking for last track state");

    auto last = findLast();
    if (!last.ok()) {
      ACTS_ERROR("no last track state found");
      return last.error();
    }

    SurfaceIntersection intersection = intersect(*last);
    if (!intersection) {
      ACTS_ERROR("no intersection found");
      return Result<std::pair<TrackStateProxy, double>>::failure(
          TrackExtrapolationError::ReferenceSurfaceUnreachable);
    }

    ACTS_VERBOSE("found intersection at " << intersection.pathLength());
    return std::make_pair(*last, intersection.pathLength());
  }

  if (strategy == TrackExtrapolationStrategy::firstOrLast) {
    ACTS_VERBOSE("looking for first or last track state");

    auto first = findFirst();
    if (!first.ok()) {
      ACTS_ERROR("no first track state found");
      return first.error();
    }

    auto last = findLast();
    if (!last.ok()) {
      ACTS_ERROR("no last track state found");
      return last.error();
    }

    SurfaceIntersection intersectionFirst = intersect(*first);
    SurfaceIntersection intersectionLast = intersect(*last);

    double absDistanceFirst = std::abs(intersectionFirst.pathLength());
    double absDistanceLast = std::abs(intersectionLast.pathLength());

    if (intersectionFirst && absDistanceFirst <= absDistanceLast) {
      ACTS_VERBOSE("using first track state with intersection at "
                   << intersectionFirst.pathLength());
      return std::make_pair(*first, intersectionFirst.pathLength());
    }

    if (intersectionLast && absDistanceLast <= absDistanceFirst) {
      ACTS_VERBOSE("using last track state with intersection at "
                   << intersectionLast.pathLength());
      return std::make_pair(*last, intersectionLast.pathLength());
    }

    ACTS_ERROR("no intersection found");
    return Result<std::pair<TrackStateProxy, double>>::failure(
        TrackExtrapolationError::ReferenceSurfaceUnreachable);
  }

  throw std::runtime_error("invalid extrapolation strategy");
}

/// @brief Extrapolate a track to a reference surface
///
/// @tparam track_proxy_t The track proxy type
/// @tparam propagator_t The propagator type
/// @tparam propagator_options_t The propagator options type
///
/// @param track The track
/// @param referenceSurface The reference surface
/// @param propagator The propagator
/// @param options The propagator options
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the extrapolation
template <typename track_proxy_t, typename propagator_t,
          typename propagator_options_t>
Result<void> extrapolateTrackToReferenceSurface(
    track_proxy_t &track, const Surface &referenceSurface,
    const propagator_t &propagator, propagator_options_t options,
    TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  auto findResult = findTrackStateForExtrapolation(
      options.geoContext, track, referenceSurface, strategy, logger);

  if (!findResult.ok()) {
    ACTS_ERROR("failed to find track state for extrapolation");
    return findResult.error();
  }

  auto &[trackState, distance] = *findResult;

  options.direction = Direction::fromScalarZeroAsPositive(distance);
  auto propagateResult = propagator.propagate(
      track.createParametersFromState(trackState), referenceSurface, options);

  if (!propagateResult.ok()) {
    ACTS_ERROR("failed to extrapolate track: " << propagateResult.error());
    return propagateResult.error();
  }

  track.setReferenceSurface(referenceSurface.getSharedPtr());
  track.parameters() = propagateResult->endParameters.value().parameters();
  track.covariance() =
      propagateResult->endParameters.value().covariance().value();

  return Result<void>::success();
}

/// @brief Extrapolate tracks to a reference surface
///
/// @tparam track_container_t The track container type
/// @tparam propagator_t The propagator type
/// @tparam propagator_options_t The propagator options type
///
/// @param trackContainer The track container
/// @param referenceSurface The reference surface
/// @param propagator The propagator
/// @param options The propagator options
/// @param strategy The extrapolation strategy
/// @param logger The logger
///
/// @return The result of the extrapolation
template <typename track_container_t, typename propagator_t,
          typename propagator_options_t>
Result<void> extrapolateTracksToReferenceSurface(
    const track_container_t &trackContainer, const Surface &referenceSurface,
    const propagator_t &propagator, propagator_options_t options,
    TrackExtrapolationStrategy strategy,
    const Logger &logger = *getDefaultLogger("TrackExtrapolation",
                                             Logging::INFO)) {
  Result<void> result = Result<void>::success();

  for (const auto &track : trackContainer) {
    auto extrapolateResult = extrapolateTrackToReferenceSurface(
        track, referenceSurface, propagator, options, strategy, logger);
    if (!extrapolateResult.ok()) {
      result = extrapolateResult.error();
      continue;
    }
  }

  return result;
}

}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::TrackExtrapolationError> : std::true_type {};
}  // namespace std
