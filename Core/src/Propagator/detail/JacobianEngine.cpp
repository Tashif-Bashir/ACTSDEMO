// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/JacobianEngine.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"

namespace Acts {

void detail::boundToBoundTransportJacobian(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters,
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the derivative of path length at the final surface or the
  // point-of-closest approach w.r.t. free parameters
  const FreeToPathMatrix freeToPath =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  // https://acts.readthedocs.io/en/latest/white_papers/correction-for-transport-jacobian.html
  // Calculate the full jacobian from the local/bound parameters at the start
  // surface to local/bound parameters at the final surface
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      freeToBoundJacobian *
      (FreeMatrix::Identity() + freeToPathDerivatives * freeToPath) *
      freeTransportJacobian * boundToFreeJacobian;
}

void detail::boundToCurvilinearTransportJacobian(
    const Vector3& direction, const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    BoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from global to local at the curvilinear surface
  FreeToBoundMatrix freeToBoundJacobian =
      CurvilinearSurface(direction).freeToBoundJacobian();

  // Update the jacobian to include the derivative of the path length at the
  // curvilinear surface w.r.t. the free parameters
  freeToBoundJacobian.topLeftCorner<6, 3>() +=
      (freeToBoundJacobian * freeToPathDerivatives) *
      (-1.0 * direction).transpose();

  // Calculate the full jocobian from the local parameters at the start surface
  // to curvilinear parameters
  // @note jac(locA->locB) = jac(gloB->locB)*(1+
  // pathCorrectionFactor(gloB))*jacTransport(gloA->gloB) *jac(locA->gloA)
  fullTransportJacobian =
      blockedMult(freeToBoundJacobian,
                  blockedMult(freeTransportJacobian, boundToFreeJacobian));
}

BoundToFreeMatrix detail::boundToFreeTransportJacobian(
    const BoundToFreeMatrix& boundToFreeJacobian,
    const FreeMatrix& freeTransportJacobian) {
  // Calculate the full jacobian, in this case simple a product of
  // jacobian(transport in free) * jacobian(bound to free)
  return freeTransportJacobian * boundToFreeJacobian;
}

void detail::freeToBoundTransportJacobian(
    const GeometryContext& geoContext, const Surface& surface,
    const FreeVector& freeParameters, const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives,
    FreeToBoundMatrix& fullTransportJacobian) {
  // Calculate the jacobian from free to bound at the final surface
  FreeToBoundMatrix freeToBoundJacobian =
      surface.freeToBoundJacobian(geoContext, freeParameters);
  FreeToPathMatrix sVec =
      surface.freeToPathDerivative(geoContext, freeParameters);
  // Return the jacobian to local
  fullTransportJacobian = freeToBoundJacobian * (freeTransportJacobian +
                                                 freeToPathDerivatives * sVec *
                                                     freeTransportJacobian);
}

FreeToBoundMatrix detail::freeToCurvilinearTransportJacobian(
    const Vector3& direction, const FreeMatrix& freeTransportJacobian,
    const FreeVector& freeToPathDerivatives) {
  auto sfactors = direction.transpose() *
                  freeTransportJacobian.template topLeftCorner<3, 8>();

  // Since the jacobian to local needs to calculated for the bound parameters
  // here, it is convenient to do the same here
  return CurvilinearSurface(direction).freeToBoundJacobian() *
         (freeTransportJacobian - freeToPathDerivatives * sfactors);
}

Result<void> detail::reinitializeJacobians(
    const GeometryContext& geoContext, const Surface& surface,
    FreeMatrix& freeTransportJacobian, FreeVector& freeToPathDerivatives,
    BoundToFreeMatrix& boundToFreeJacobian, const FreeVector& freeParameters) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();

  // Get the local position
  const Vector3 position = freeParameters.segment<3>(eFreePos0);
  const Vector3 direction = freeParameters.segment<3>(eFreeDir0);
  auto lpResult = surface.globalToLocal(geoContext, position, direction);
  if (!lpResult.ok()) {
    return lpResult.error();
  }
  // Transform from free to bound parameters
  Result<BoundVector> boundParameters = detail::transformFreeToBoundParameters(
      freeParameters, surface, geoContext);
  if (!boundParameters.ok()) {
    return boundParameters.error();
  }
  // Reset the jacobian from local to global
  boundToFreeJacobian =
      surface.boundToFreeJacobian(geoContext, *boundParameters);
  return Result<void>::success();
}

void detail::reinitializeJacobians(FreeMatrix& freeTransportJacobian,
                                   FreeVector& freeToPathDerivatives,
                                   BoundToFreeMatrix& boundToFreeJacobian,
                                   const Vector3& direction) {
  // Reset the jacobians
  freeTransportJacobian = FreeMatrix::Identity();
  freeToPathDerivatives = FreeVector::Zero();
  boundToFreeJacobian = CurvilinearSurface(direction).boundToFreeJacobian();
}

}  // namespace Acts
