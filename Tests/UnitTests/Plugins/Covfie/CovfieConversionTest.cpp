// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// system includes
#include <cmath>
#include <iostream>

// covfie core
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/parameter_pack.hpp>


// acts includes
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Definitions/Units.hpp"

// covfie plugin
#include "Acts/Plugins/Covfie/CovfieConversion.hpp"

// Boost.Test include(s)
#include <boost/test/unit_test.hpp>

using namespace Acts::UnitLiterals;

template <typename cache_t, typename view_t, typename iterator_t>
bool MagneticFieldEqual(const Acts::MagneticFieldProvider& fieldProvider, cache_t& cache, view_t view, iterator_t points, float error_margin_half_width){
    for (auto point : points){
        auto x = point[0], y = point[1], z = point[2];

        auto lookupResult = fieldProvider.getField(Acts::Vector3{x, y, z}, cache);
        if(!lookupResult.ok()) {
            throw std::runtime_error{"Field lookup failure"};
        }
        auto actsValueX = (*lookupResult)[0], actsValueY = (*lookupResult)[1], actsValueZ = (*lookupResult)[2];

        auto covfieValues = view.at(x, y, z);
        auto covfieValueX = covfieValues[0], covfieValueY = covfieValues[1], covfieValueZ = covfieValues[2];

        if (error_margin_half_width < std::abs(covfieValueX - actsValueX) || error_margin_half_width < std::abs(covfieValueY - actsValueY) || error_margin_half_width < std::abs(covfieValueZ - actsValueZ)){
            std::cout << "Failed at value (" << x << ", " << y << ", " << z << ") (acts: (" << actsValueX << ", " << actsValueY << ", " << actsValueZ << ")) (covfie: (" << covfieValueX << ", " << covfieValueY << ", " << covfieValueZ << "))" << std::endl;
            return false;
        }
        std::cout << "Success at value (" << x << ", " << y << ", " << z << ") (acts: (" << actsValueX << ", " << actsValueY << ", " << actsValueZ << ")) (covfie: (" << covfieValueX << ", " << covfieValueY << ", " << covfieValueZ << "))" << std::endl;
    }
    return true;
}

BOOST_AUTO_TEST_CASE(Covfie_Conversion_InterpolatedMagneticField1) {
    auto localToGlobalBin_xyz = [](std::array<std::size_t, 3> binsXYZ,
                                    std::array<std::size_t, 3> nBinsXYZ) {
        return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2)) +
                binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
    };

    std::vector<double> xPos = {0., 1., 2., 3.};
    std::vector<double> yPos = {0., 1., 2., 3.};
    std::vector<double> zPos = {0., 1., 2., 3.};

    std::vector<Acts::Vector3> bField_xyz;
    for (int i = 0; i < 64; i++) {
        bField_xyz.push_back(Acts::Vector3(i, i, i));
    }

    Acts::MagneticFieldContext fieldContext;
    auto actsField = Acts::fieldMapXYZ(localToGlobalBin_xyz, xPos, yPos, zPos,
                                    bField_xyz, 1, 1, false);
    auto cache = actsField.makeCache(fieldContext);

    auto field = Acts::CovfiePlugin::covfieField(actsField);
    auto view = typename decltype(field)::view_t(field);

    std::array<std::array<float, 3>, 13> points = {{
        {0.f, 0.f, 0.f},
        {1.f, 1.f, 1.f},
        {2.f, 2.f, 2.f},
        {1.2f, 2.5f, 0.8f},
        {0.7f, 1.9f, 2.3f},
        {2.1f, 0.3f, 1.5f},
        {0.4f, 2.8f, 2.9f},
        {1.6f, 1.2f, 0.5f},
        {2.3f, 0.6f, 2.2f},
        {1.1f, 2.7f, 1.3f},
        {0.9f, 1.4f, 2.7f},
        {2.4f, 1.8f, 0.9f},
        {0.6f, 2.2f, 2.1f},
    }};

    BOOST_TEST(MagneticFieldEqual(actsField, cache, view, points, 0.0001));
}

BOOST_AUTO_TEST_CASE(Covfie_Conversion_InterpolatedMagneticField2) {
    auto localToGlobalBin_xyz = [](std::array<std::size_t, 3> binsXYZ,
                                    std::array<std::size_t, 3> nBinsXYZ) {
        return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2)) +
                binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
    };

    std::vector<double> xPos = {8., 12., 16., 20.};
    std::vector<double> yPos = {8., 12., 16., 20.};
    std::vector<double> zPos = {8., 12., 16., 20.};

    std::vector<Acts::Vector3> bField_xyz;
    for (int i = 0; i < 64; i++) {
        bField_xyz.push_back(Acts::Vector3(i, i*i*0.01, i));
    }

    Acts::MagneticFieldContext fieldContext;
    auto actsField = Acts::fieldMapXYZ(localToGlobalBin_xyz, xPos, yPos, zPos,
                                    bField_xyz, 1, 1, false);
    auto cache = actsField.makeCache(fieldContext);

    auto field = Acts::CovfiePlugin::covfieField(actsField);
    auto view = typename decltype(field)::view_t(field);

    std::array<std::array<float, 3>, 13> points = {{
        {8.f, 8.f, 8.f},
        {12.f, 12.f, 12.f},
        {16.f, 16.f, 16.f},
        {8.1f, 10.2f, 12.3f},
        {9.4f, 11.5f, 13.6f},
        {10.7f, 12.8f, 14.9f},
        {11.0f, 13.1f, 15.2f},
        {12.3f, 14.4f, 16.5f},
        {13.6f, 15.7f, 17.8f},
        {14.9f, 16.0f, 18.1f},
        {16.2f, 17.3f, 19.4f},
        {17.5f, 18.6f, 19.7f},
        {18.8f, 19.9f, 14.0f},
    }};

    BOOST_TEST(MagneticFieldEqual(actsField, cache, view, points, 0.0001f));
}

BOOST_AUTO_TEST_CASE(Covfie_Conversion_ConstantMagneticField1) {

    Acts::ConstantBField actsField(Acts::Vector3{1.3f, 2.5f, 2.f});
    Acts::MagneticFieldContext ctx;
    auto cache = actsField.makeCache(ctx);

    auto field = Acts::CovfiePlugin::covfieField(actsField);
    auto view = typename decltype(field)::view_t(field);

    std::array<std::array<float, 3>, 13> points = {{
        {8.f, 8.f, 8.f},
        {12.f, 12.f, 12.f},
        {16.f, 16.f, 16.f},
        {8.1f, 10.2f, 12.3f},
        {9.4f, 11.5f, 13.6f},
        {10.7f, 12.8f, 14.9f},
        {11.0f, 13.1f, 15.2f},
        {12.3f, 14.4f, 16.5f},
        {13.6f, 15.7f, 17.8f},
        {14.9f, 16.0f, 18.1f},
        {16.2f, 17.3f, 19.4f},
        {17.5f, 18.6f, 19.7f},
        {18.8f, 19.9f, 14.0f},
    }};

    BOOST_TEST(MagneticFieldEqual(actsField, cache, view, points, 0.0901));
}

BOOST_AUTO_TEST_CASE(Covfie_Conversion_MagneticFieldProvider1) {

    Acts::SolenoidBField::Config cfg;
    cfg.length = 5.8_m;
    cfg.radius = (2.56 + 2.46) * 0.5 * 0.5_m;
    cfg.nCoils = 1154;
    cfg.bMagCenter = 2_T;
    Acts::SolenoidBField actsField(cfg);
    Acts::MagneticFieldContext ctx;
    auto cache = actsField.makeCache(ctx);

    auto field = Acts::CovfiePlugin::covfieField(actsField, cache, std::vector{20UL, 20UL, 20UL}, std::vector{0., 0., 0.}, std::vector{20., 20., 20.});
    auto view = typename decltype(field)::view_t(field);

    std::array<std::array<float, 3>, 13> points = {{
        {8.f, 8.f, 8.f},
        {12.f, 12.f, 12.f},
        {16.f, 16.f, 16.f},
        {8.1f, 10.2f, 12.3f},
        {9.4f, 11.5f, 13.6f},
        {10.7f, 12.8f, 14.9f},
        {11.0f, 13.1f, 15.2f},
        {12.3f, 14.4f, 16.5f},
        {13.6f, 15.7f, 17.8f},
        {14.9f, 16.0f, 18.1f},
        {16.2f, 17.3f, 19.4f},
        {17.5f, 18.6f, 19.7f},
        {18.8f, 19.9f, 14.0f},
    }};

    BOOST_TEST(MagneticFieldEqual(actsField, cache, view, points, 0.0001));
}
