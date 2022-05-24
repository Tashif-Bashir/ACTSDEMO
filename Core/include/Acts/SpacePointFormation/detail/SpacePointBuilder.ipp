// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {
namespace detail {

/// @brief Storage container for variables related to the calculation of space
/// points
struct SpacePointParameters {
  /// Vector pointing from bottom to top end of first SDE
  Vector3 q;
  /// Vector pointing from bottom to top end of second SDE
  Vector3 r;
  /// Twice the vector pointing from vertex to to midpoint of first SDE
  Vector3 s;
  /// Twice the vector pointing from vertex to to midpoint of second SDE
  Vector3 t;
  /// Cross product between SpacePointParameters::q and
  /// SpacePointParameters::s
  Vector3 qs;
  /// Cross product between SpacePointParameters::r and
  /// SpacePointParameters::t
  Vector3 rt;
  /// Magnitude of SpacePointParameters::q
  double qmag = 0.;
  /// Parameter that determines the hit position on the first SDE
  double m = 0.;
  /// Parameter that determines the hit position on the second SDE
  double n = 0.;
  /// Regular limit of the absolut values of SpacePointParameters::m and
  /// SpacePointParameters::n
  double limit = 1.;
  /// Limit of SpacePointParameters::m and SpacePointParameters::n in case of
  /// variable vertex
  double limitExtended = 0.;
};

/// @brief Calculates (Delta theta)^2 + (Delta phi)^2 between two measurements
///
/// @param [in] pos1 position of the first measurement
/// @param [in] pos2 position the second measurement
/// @param [in] maxDistance Maximum distance between two measurements
/// @param [in] maxAngleTheta2 Maximum squared theta angle between two
/// measurements
/// @param [in] maxAnglePhi2 Maximum squared phi angle between two measurements
///
/// @return The squared sum within configuration parameters, otherwise -1
inline double differenceOfMeasurementsChecked(const Vector3& pos1,
                                              const Vector3& pos2,
                                              const Vector3& posVertex,
                                              const double maxDistance,
                                              const double maxAngleTheta2,
                                              const double maxAnglePhi2) {
  // Check if measurements are close enough to each other
  if ((pos1 - pos2).norm() > maxDistance) {
    return -1.;
  }

  // Calculate the angles of the vectors
  double phi1 = VectorHelpers::phi(pos1 - posVertex);
  double theta1 = VectorHelpers::theta(pos1 - posVertex);
  double phi2 = VectorHelpers::phi(pos2 - posVertex);
  double theta2 = VectorHelpers::theta(pos2 - posVertex);

  // Calculate the squared difference between the theta angles
  double diffTheta2 = (theta1 - theta2) * (theta1 - theta2);
  if (diffTheta2 > maxAngleTheta2) {
    return -1.;
  }
  // Calculate the squared difference between the phi angles
  double diffPhi2 = (phi1 - phi2) * (phi1 - phi2);
  if (diffPhi2 > maxAnglePhi2) {
    return -1.;
  }
  // Return the squared distance between both vector
  return diffTheta2 + diffPhi2;
}

/// @brief Calculates a space point whithout using the vertex
/// @note This is mostly to resolve space points from cosmic data
/// @param a vector to the top end of the first SDE
/// @param c vector to the top end of the second SDE
/// @param q vector from the bottom to the top end of the first SDE
/// @param r vector from the bottom to the top end of the second SDE
/// @return parameter that indicates the location of the space point; returns
/// 1. if it failed
/// @note The meaning of the parameter is explained in more detail in the
/// function body
inline double calcPerpendicularProjection(const Vector3& a, const Vector3& c,
                                          const Vector3& q, const Vector3& r) {
  /// This approach assumes that no vertex is available. This option aims to
  /// approximate the space points from cosmic data.
  /// The underlying assumption is that the best point is given by the  closest
  /// distance between both lines describing the SDEs.
  /// The point x on the first SDE is parametrized as a + lambda0 * q with  the
  /// top end a of the strip and the vector q = a - b(ottom end of the  strip).
  /// An analogous parametrization is performed of the second SDE with y = c  +
  /// lambda1 * r.
  /// x get resolved by resolving lambda0 from the condition that |x-y| is  the
  /// shortest distance between two skew lines.

  Vector3 ac = c - a;
  double qr = q.dot(r);
  double denom = q.dot(q) - qr * qr;

  // Check for numerical stability
  if (fabs(denom) > 10e-7) {
    // Return lambda0
    return (ac.dot(r) * qr - ac.dot(q) * r.dot(r)) / denom;
  }
  // lambda0 is in the interval [-1,0]. This return serves as error check.
  return 1.;
}

/// @brief This function tests if a space point can be estimated by a more
/// tolerant treatment of construction. In fact, this function indirectly
/// allows shifts of the vertex.
///
/// @param [in] spaPoPa container that stores geometric parameters and rules of
/// the space point formation
/// @param [in] stripLengthGapTolerance Tolerance scaling factor of the gap
/// between strip detector elements
///
/// @return indicator if the test was successful
inline bool recoverSpacePoint(SpacePointParameters& spaPoPa,
                              double stripLengthGapTolerance) {
  /// Consider some cases that would allow an easy exit
  // Check if the limits are allowed to be increased
  if (stripLengthGapTolerance <= 0.) {
    return false;
  }
  spaPoPa.qmag = spaPoPa.q.norm();
  // Increase the limits. This allows a check if the point is just slightly
  // outside the SDE
  spaPoPa.limitExtended =
      spaPoPa.limit + stripLengthGapTolerance / spaPoPa.qmag;
  // Check if m is just slightly outside
  if (fabs(spaPoPa.m) > spaPoPa.limitExtended) {
    return false;
  }
  // Calculate n if not performed previously
  if (spaPoPa.n == 0.) {
    spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) / spaPoPa.r.dot(spaPoPa.qs);
  }
  // Check if n is just slightly outside
  if (fabs(spaPoPa.n) > spaPoPa.limitExtended) {
    return false;
  }

  /// The following code considers an overshoot of m and n in the same direction
  /// of their SDE. The term "overshoot" represents the amount of m or n outside
  /// its regular interval (-1, 1).
  /// It calculates which overshoot is worse. In order to compare both, the
  /// overshoot in n is projected onto the first surface by considering the
  /// normalized projection of r onto q.
  /// This allows a rescaling of the overshoot. The worse overshoot will be set
  /// to +/-1, the parameter with less overshoot will be moved towards 0 by the
  /// worse overshoot.
  /// In order to treat both SDEs equally, the rescaling eventually needs to be
  /// performed several times. If these shifts allows m and n to be in the
  /// limits, the space point can be stored.
  /// @note This shift can be understood as a shift of the particle's
  /// trajectory. This is leads to a shift of the vertex. Since these two points
  /// are treated independently from other measurement, it is also possible to
  /// consider this as a change in the slope of the particle's trajectory.
  ///  The would also move the vertex position.

  // Calculate the scaling factor to project lengths of the second SDE on the
  // first SDE
  double secOnFirstScale =
      spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.qmag * spaPoPa.qmag);
  // Check if both overshoots are in the same direction
  if (spaPoPa.m > 1. && spaPoPa.n > 1.) {
    // Calculate the overshoots
    double mOvershoot = spaPoPa.m - 1.;
    double nOvershoot =
        (spaPoPa.n - 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m -= biggerOvershoot;
    spaPoPa.n -= (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // Check if both overshoots are in the same direction
  if (spaPoPa.m < -1. && spaPoPa.n < -1.) {
    // Calculate the overshoots
    double mOvershoot = -(spaPoPa.m + 1.);
    double nOvershoot =
        -(spaPoPa.n + 1.) * secOnFirstScale;  // Perform projection
    // Resolve worse overshoot
    double biggerOvershoot = std::max(mOvershoot, nOvershoot);
    // Move m and n towards 0
    spaPoPa.m += biggerOvershoot;
    spaPoPa.n += (biggerOvershoot / secOnFirstScale);
    // Check if this recovered the space point
    return fabs(spaPoPa.m) < spaPoPa.limit && fabs(spaPoPa.n) < spaPoPa.limit;
  }
  // No solution could be found
  return false;
}

/// @brief This function performs a straight forward calculation of a space
/// point and returns whether it was succesful or not.
///
/// @param [in] stripEnds1 Top and bottom end of the first strip detector
/// element
/// @param [in] stripEnds1 Top and bottom end of the second strip detector
/// element
/// @param [in] posVertex Position of the vertex
/// @param [in, out] spaPoPa Data container of the calculations
/// @param [in] stripLengthTolerance Tolerance scaling factor on the strip
/// detector element length
///
/// @return Boolean statement whether the space point calculation was succesful
inline bool calculateSpacePointPosition(
    const std::pair<Vector3, Vector3>& stripEnds1,
    const std::pair<Vector3, Vector3>& stripEnds2, const Vector3& posVertex,
    SpacePointParameters& spaPoPa, const double stripLengthTolerance) {
  /// The following algorithm is meant for finding the position on the first
  /// strip if there is a corresponding Measurement on the second strip. The
  /// resulting point is a point x on the first surfaces. This point is
  /// along a line between the points a (top end of the strip)
  /// and b (bottom end of the strip). The location can be parametrized as
  /// 	2 * x = (1 + m) a + (1 - m) b
  /// as function of the scalar m. m is a parameter in the interval
  /// -1 < m < 1 since the hit was on the strip. Furthermore, the vector
  /// from the vertex to the Measurement on the second strip y is needed to be a
  /// multiple k of the vector from vertex to the hit on the first strip x.
  /// As a consequence of this demand y = k * x needs to be on the
  /// connecting line between the top (c) and bottom (d) end of
  /// the second strip. If both measurements correspond to each other, the
  /// condition
  /// 	y * (c X d) = k * x (c X d) = 0 ("X" represents a cross product)
  /// needs to be fulfilled. Inserting the first equation into this
  /// equation leads to the condition for m as given in the following
  /// algorithm and therefore to the calculation of x.
  /// The same calculation can be repeated for y. Its corresponding
  /// parameter will be named n.

  spaPoPa.s = stripEnds1.first + stripEnds1.second - 2 * posVertex;
  spaPoPa.t = stripEnds2.first + stripEnds2.second - 2 * posVertex;
  spaPoPa.qs = spaPoPa.q.cross(spaPoPa.s);
  spaPoPa.rt = spaPoPa.r.cross(spaPoPa.t);
  spaPoPa.m = -spaPoPa.s.dot(spaPoPa.rt) / spaPoPa.q.dot(spaPoPa.rt);

  // Set the limit for the parameter
  if (spaPoPa.limit == 1. && stripLengthTolerance != 0.) {
    spaPoPa.limit = 1. + stripLengthTolerance;
  }

  // Check if m and n can be resolved in the interval (-1, 1)
  return (fabs(spaPoPa.m) <= spaPoPa.limit &&
          fabs(spaPoPa.n = -spaPoPa.t.dot(spaPoPa.qs) /
                           spaPoPa.r.dot(spaPoPa.qs)) <= spaPoPa.limit);
}
}  // namespace detail

template <typename spacepoint_t>
SpacePointBuilder<spacepoint_t>::SpacePointBuilder(
    SpacePointBuilderConfig cfg,
    std::function<
        spacepoint_t(Acts::Vector3, Acts::Vector2,
                     boost::container::static_vector<const SourceLink*, 2>)>
        func,
    std::unique_ptr<const Logger> logger)
    : m_config(cfg), m_spConstructor(func), m_logger(std::move(logger)) {
  m_spUtility = std::make_shared<SpacePointUtility>(cfg);
}

template <typename spacepoint_t>
template <template <typename...> typename container_t>
void SpacePointBuilder<spacepoint_t>::calculateSingleHitSpacePoints(
    const GeometryContext& gctx,
    const std::vector<const Measurement*>& measurements,
    std::back_insert_iterator<container_t<spacepoint_t>> spacePointIt) const {
  for (const auto& meas : measurements) {
    auto [gPos, gCov] = m_spUtility->globalCoords(gctx, *meas);
    const auto& slink =
        std::visit([](const auto& x) { return &x.sourceLink(); }, *meas);

    boost::container::static_vector<const SourceLink*, 2> slinks;
    slinks.emplace_back(slink);
    spacePointIt = m_spConstructor(gPos, gCov, std::move(slinks));
  }
}

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::makeMeasurementPairs(
    const GeometryContext& gctx,
    const std::vector<const Measurement*>& measurementsFront,
    const std::vector<const Measurement*>& measurementsBack,
    std::vector<std::pair<const Measurement*, const Measurement*>>&
        measurementPairs) const {
  // Return if no Measurements are given in a vector
  if (measurementsFront.empty() || measurementsBack.empty()) {
    return;
  }
  // Declare helper variables
  double currentDiff;
  double diffMin;
  unsigned int measurementMinDist;

  // Walk through all Measurements on both surfaces
  for (unsigned int iMeasurementsFront = 0;
       iMeasurementsFront < measurementsFront.size(); iMeasurementsFront++) {
    // Set the closest distance to the maximum of double
    diffMin = std::numeric_limits<double>::max();
    // Set the corresponding index to an element not in the list of Measurements
    measurementMinDist = measurementsBack.size();
    for (unsigned int iMeasurementsBack = 0;
         iMeasurementsBack < measurementsBack.size(); iMeasurementsBack++) {
      auto [gposFront, gcovFront] = m_spUtility->globalCoords(
          gctx, *(measurementsFront[iMeasurementsFront]));
      auto [gposBack, gcovBack] = m_spUtility->globalCoords(
          gctx, *(measurementsBack[iMeasurementsBack]));

      currentDiff = detail::differenceOfMeasurementsChecked(
          gposFront, gposBack, m_config.vertex, m_config.diffDist,
          m_config.diffPhi2, m_config.diffTheta2);

      // Store the closest Measurements (distance and index) calculated so far
      if (currentDiff < diffMin && currentDiff >= 0.) {
        diffMin = currentDiff;
        measurementMinDist = iMeasurementsBack;
      }
    }

    // Store the best (=closest) result
    if (measurementMinDist < measurementsBack.size()) {
      std::pair<const Measurement*, const Measurement*> measurementPair =
          std::make_pair(measurementsFront[iMeasurementsFront],
                         measurementsBack[measurementMinDist]);
      measurementPairs.push_back(measurementPair);
    }
  }
}

template <typename spacepoint_t>
void SpacePointBuilder<spacepoint_t>::calculateDoubleHitSpacePoint(
    const GeometryContext& gctx,
    const std::pair<const Measurement*, const Measurement*>& measurementPair,
    const std::pair<const std::pair<Vector3, Vector3>,
                    const std::pair<Vector3, Vector3>>& stripEndsPair,
    std::shared_ptr<const spacepoint_t> spacePoint) const {
  // Source of algorithm: Athena, SiSpacePointMakerTool::makeSCT_SpacePoint()

  detail::SpacePointParameters spaPoPa;

  const auto& ends1 = stripEndsPair.first;
  const auto& ends2 = stripEndsPair.second;

  spaPoPa.q = ends1.first - ends1.second;
  spaPoPa.r = ends2.first - ends2.second;

  const auto& slink_front =
      m_spUtility->getSourceLink(*(measurementPair.first));
  const auto& slink_back =
      m_spUtility->getSourceLink(*(measurementPair.second));

  boost::container::static_vector<const SourceLink*, 2> slinks;
  slinks.emplace_back(slink_front);
  slinks.emplace_back(slink_back);

  if (m_config.usePerpProj) {  // for cosmic without vertex constraint
    double resultPerpProj = detail::calcPerpendicularProjection(
        ends1.first, ends2.first, spaPoPa.q, spaPoPa.r);
    if (resultPerpProj <= 0.) {
      Vector3 gPos = ends1.first + resultPerpProj * spaPoPa.q;
      double theta = acos(spaPoPa.q.dot(spaPoPa.r) /
                          (spaPoPa.q.norm() * spaPoPa.r.norm()));
      const auto gCov = m_spUtility->calcGlobalVars(
          gctx, *(measurementPair.first), *(measurementPair.second), theta);

      spacePoint = std::make_shared<spacepoint_t>(
          m_spConstructor(gPos, gCov, std::move(slinks)));
    }
  }

  bool spFound = calculateSpacePointPosition(
      ends1, ends2, m_config.vertex, spaPoPa, m_config.stripLengthTolerance);
  if (!spFound)
    spFound =
        detail::recoverSpacePoint(spaPoPa, m_config.stripLengthGapTolerance);

  if (spFound) {
    Vector3 gPos = 0.5 * (ends1.first + ends1.second + spaPoPa.m * spaPoPa.q);

    double theta =
        acos(spaPoPa.q.dot(spaPoPa.r) / (spaPoPa.q.norm() * spaPoPa.r.norm()));

    const auto gCov = m_spUtility->calcGlobalVars(
        gctx, *(measurementPair.first), *(measurementPair.second), theta);

    spacePoint = std::make_shared<spacepoint_t>(
        m_spConstructor(gPos, gCov, std::move(slinks)));
  }
}

}  // namespace Acts
