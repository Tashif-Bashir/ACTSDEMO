// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ExaTrkX/TorchMetricLearning.hpp"

#include "Acts/Plugins/ExaTrkX/buildEdges.hpp"

#include <torch/script.h>
#include <torch/torch.h>

#include "printCudaMemInfo.hpp"

using namespace torch::indexing;

namespace Acts {

TorchMetricLearning::TorchMetricLearning(const Config &cfg,
                                         std::unique_ptr<const Logger> _logger)
    : m_logger(std::move(_logger)), m_cfg(cfg) {
  c10::InferenceMode guard(true);
  m_deviceType = torch::cuda::is_available() ? torch::kCUDA : torch::kCPU;
  ACTS_DEBUG("Using torch version " << TORCH_VERSION_MAJOR << "."
                                    << TORCH_VERSION_MINOR << "."
                                    << TORCH_VERSION_PATCH);
#ifndef ACTS_EXATRKX_CPUONLY
  if (not torch::cuda::is_available()) {
    ACTS_INFO("CUDA not available, falling back to CPU");
  }
#endif

  try {
    m_model = std::make_unique<torch::jit::Module>();
    *m_model = torch::jit::load(m_cfg.modelPath, m_deviceType);
    m_model->eval();
  } catch (const c10::Error &e) {
    throw std::invalid_argument("Failed to load models: " + e.msg());
  }
}

TorchMetricLearning::~TorchMetricLearning() {}

std::tuple<std::any, std::any> TorchMetricLearning::operator()(
    boost::multi_array<float, 2> &inputValues) {
  ACTS_DEBUG("Start graph construction");
  c10::InferenceMode guard(true);
  const torch::Device device(m_deviceType);

  const int64_t numSpacepoints = inputValues.shape()[0];
  const int64_t numAllFeatures = inputValues.shape()[1];

  // printout the r,phi,z of the first spacepoint
  ACTS_VERBOSE("First spacepoint information: " << [&]() {
    std::stringstream ss;
    for (int i = 0; i < numAllFeatures; ++i) {
      ss << inputValues[0][i] << "  ";
    }
    return ss.str();
  }());
  printCudaMemInfo(logger());

  auto opts = torch::TensorOptions().dtype(torch::kFloat32);
  torch::Tensor inputTensor =
      torch::from_blob(inputValues.data(), {numSpacepoints, numAllFeatures},
                       opts)
          .to(device);

  // **********
  // Embedding
  // **********

  if (m_cfg.numFeatures > numAllFeatures) {
    throw std::runtime_error("requested more features then available");
  }

  // Clone models (solve memory leak? members can be const...)
  auto model = m_model->clone();
  model.to(device);

  std::vector<torch::jit::IValue> inputTensors;
  inputTensors.push_back(
      m_cfg.numFeatures < numAllFeatures
          ? inputTensor.index({Slice{}, Slice{None, m_cfg.numFeatures}})
          : std::move(inputTensor));

  ACTS_DEBUG("embedding input tensor shape "
             << inputTensors[0].toTensor().size(0) << ", "
             << inputTensors[0].toTensor().size(1));

  auto output = model.forward(inputTensors).toTensor();

  ACTS_VERBOSE("Embedding space of the first SP:\n"
               << output.slice(/*dim=*/0, /*start=*/0, /*end=*/1));
  printCudaMemInfo(logger());

  // ****************
  // Building Edges
  // ****************
  auto edgeList = buildEdges(output, m_cfg.rVal, m_cfg.knnVal);

  ACTS_VERBOSE("Shape of built edges: (" << edgeList.size(0) << ", "
                                         << edgeList.size(1));
  ACTS_VERBOSE("Slice of edgelist:\n" << edgeList.slice(1, 0, 5));
  printCudaMemInfo(logger());

  return {std::move(inputTensors[0]).toTensor(), std::move(edgeList)};
}
}  // namespace Acts
