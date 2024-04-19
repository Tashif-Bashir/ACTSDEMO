// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/parameter_pack.hpp>

// acts includes
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

namespace Acts::CovfieConversion{

using builder_backend_t = covfie::backend::strided<
    covfie::vector::size3,
    covfie::backend::array<covfie::vector::float3>>;

using affine_linear_strided_field_t = covfie::field<covfie::backend::affine<covfie::backend::linear<builder_backend_t>>>;

using constant_field_t = covfie::field<covfie::backend::constant<covfie::vector::size3, covfie::vector::float3>>;

/// @brief Creates a covfie field from an interpolated magnetic field.
/// @param magneticField The acts interpolated magnetic field.
/// @return An affine linear strided covfie field.
affine_linear_strided_field_t covfieField(const Acts::InterpolatedMagneticField& magneticField);

/// @brief Creates a covfie field from a constant B field.
/// @param magneticField The acts constant magnetic field.
/// @return A constant covfie field.
constant_field_t covfieField(const Acts::ConstantBField& magneticField);

/// @brief Creates a covfie field from a magnetic field provider by sampling it.
/// @param magneticField The acts magnetic field provider.
/// @param cache The acts cache.
/// @param nBins 3D array of containing the number of bins for each axis.
/// @param min (min_x, min_y, min_z)
/// @param max (max_x, max_y, max_z)
/// @return An affine linear strided covfie field.
affine_linear_strided_field_t covfieField(const Acts::MagneticFieldProvider& magneticField, Acts::MagneticFieldProvider::Cache& cache, const std::vector<std::size_t>& nBins, const std::vector<double>& min, const std::vector<double>& max);

};
