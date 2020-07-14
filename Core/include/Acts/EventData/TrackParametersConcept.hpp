// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

class Surface;

namespace Concepts {
using namespace Acts::concept;

// nested types that must be available
template <typename T>
using TypeParametersVector = typename T::ParVector_t;
template <typename T>
using TypeCovarianceMatrix = typename T::CovMatrix_t;

template <typename T>
using ReturnTypeReferenceSurface =
    decltype(std::declval<T>().referenceSurface());
template <typename T>
using ReturnTypeParameters = decltype(std::declval<T>().parameters());
template <typename T>
using ReturnTypeCovariance = decltype(std::declval<T>().covariance());
template <typename T>
using ReturnTypePosition = decltype(std::declval<T>().position());
template <typename T>
using ReturnTypeTime = decltype(std::declval<T>().time());
template <typename T>
using ReturnTypeMomentum = decltype(std::declval<T>().momentum());
template <typename T>
using ReturnTypeCharge = decltype(std::declval<T>().charge());

template <typename T>
struct BoundTrackParametersConceptImpl {
  // check for required nested types
  constexpr static bool hasTypeParametersVector =
      exists<TypeParametersVector, const T>;
  constexpr static bool hasTypeCovarianceMatrix =
      exists<TypeCovarianceMatrix, const T>;

  // check for required methods
  constexpr static bool hasMethodReferenceSurface =
      identical_to<const Surface&, ReturnTypeReferenceSurface, const T>;
  constexpr static bool hasMethodParameters =
      identical_to<BoundVector, ReturnTypeParameters, const T>;
  constexpr static bool hasMethodCovariance =
      identical_to<const std::optional<BoundSymMatrix>&, ReturnTypeCovariance,
                   const T>;
  constexpr static bool hasMethodPosition =
      identical_to<Vector3D, ReturnTypePosition, const T>;
  constexpr static bool hasMethodTime =
      identical_to<double, ReturnTypeTime, const T>;
  constexpr static bool hasMethodMomentum =
      identical_to<Vector3D, ReturnTypeMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<double, ReturnTypeCharge, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeParametersVector, "Parameters vector type is missing");
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodReferenceSurface,
                "Missing/ invalid 'referenceSurface' method");
  static_assert(hasMethodParameters, "Missing/ invvalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing/ invvalid 'covariance' method");
  static_assert(hasMethodPosition, "Missing/ invvalid 'position' method");
  static_assert(hasMethodTime, "Missing/ invvalid 'time' method");
  static_assert(hasMethodMomentum, "Missing/ invvalid 'momentum' method");
  static_assert(hasMethodCharge, "Missing/ invvalid 'charge' method");

  constexpr static bool value =
      require<hasTypeParametersVector, hasTypeCovarianceMatrix,
              hasMethodReferenceSurface, hasMethodParameters,
              hasMethodCovariance, hasMethodPosition, hasMethodTime,
              hasMethodMomentum, hasMethodCharge>;
};

template <typename T>
struct FreeTrackParametersConceptImpl {
  // check for required nested types
  constexpr static bool hasTypeCovarianceMatrix =
      exists<TypeCovarianceMatrix, const T>;

  // check for required methods
  constexpr static bool hasMethodParameters =
      identical_to<FreeVector, ReturnTypeParameters, const T>;
  constexpr static bool hasMethodCovariance =
      identical_to<const std::optional<FreeSymMatrix>&, ReturnTypeCovariance,
                   const T>;
  constexpr static bool hasMethodPosition =
      identical_to<Vector3D, ReturnTypePosition, const T>;
  constexpr static bool hasMethodTime =
      identical_to<double, ReturnTypeTime, const T>;
  constexpr static bool hasMethodMomentum =
      identical_to<Vector3D, ReturnTypeMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<double, ReturnTypeCharge, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodParameters, "Missing/ invvalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing/ invvalid 'covariance' method");
  static_assert(hasMethodPosition, "Missing/ invvalid 'position' method");
  static_assert(hasMethodTime, "Missing/ invvalid 'time' method");
  static_assert(hasMethodMomentum, "Missing/ invvalid 'momentum' method");
  static_assert(hasMethodCharge, "Missing/ invvalid 'charge' method");

  constexpr static bool value =
      require<hasTypeCovarianceMatrix, hasMethodParameters, hasMethodCovariance,
              hasMethodPosition, hasMethodTime, hasMethodMomentum,
              hasMethodCharge>;
};

template <typename parameters_t>
constexpr bool BoundTrackParametersConcept =
    BoundTrackParametersConceptImpl<parameters_t>::value;

template <typename parameters_t>
constexpr bool FreeTrackParametersConcept =
    FreeTrackParametersConceptImpl<parameters_t>::value;

}  // namespace Concepts
}  // namespace Acts
