// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/AmbiguityResolution/GreedyAmbiguityResolution.hpp"

#include <unordered_map>

namespace Acts {

template <typename track_container_t, typename traj_t,
          template <typename> class holder_t, typename source_link_hash_t,
          typename source_link_equality_t>
void GreedyAmbiguityResolution::computeInitialState(
    const TrackContainer<track_container_t, traj_t, holder_t>& tracks,
    State& state, source_link_hash_t&& sourceLinkHash,
    source_link_equality_t&& sourceLinkEquality) const {
  auto measurementIndexMap =
      std::unordered_map<SourceLink, std::size_t, source_link_hash_t,
                         source_link_equality_t>(0, sourceLinkHash,
                                                 sourceLinkEquality);

  for (const auto& track : tracks) {
    if (track.nMeasurements() < m_cfg.nMeasurementsMin) {
      continue;
    }
    std::vector<std::size_t> measurements;
    tracks.trackStateContainer().visitBackwards(
        track.tipIndex(), [&](const auto& hit) {
          if (hit.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
            SourceLink sourceLink = hit.getUncalibratedSourceLink();
            // assign a new measurement index if the source link was not seen
            // yet
            auto emplace = measurementIndexMap.try_emplace(
                sourceLink, measurementIndexMap.size());
            measurements.push_back(emplace.first->second);
          }
          return true;
        });

    state.trackTips.push_back(track.index());
    state.trackChi2.push_back(track.chi2());
    state.measurementsPerTrack.push_back(std::move(measurements));
    state.selectedTracks.insert(state.numberOfTracks);

    ++state.numberOfTracks;
  }

  for (std::size_t iTrack = 0; iTrack < state.numberOfTracks; ++iTrack) {
    for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
      state.tracksPerMeasurement[iMeasurement].insert(iTrack);
    }
  }
  state.sharedMeasurementsPerTrack =
      std::vector<std::size_t>(state.trackTips.size(), 0);

  for (std::size_t iTrack = 0; iTrack < state.numberOfTracks; ++iTrack) {
    for (auto iMeasurement : state.measurementsPerTrack[iTrack]) {
      if (state.tracksPerMeasurement[iMeasurement].size() > 1) {
        ++state.sharedMeasurementsPerTrack[iTrack];
      }
    }
  }
}

}  // namespace Acts
