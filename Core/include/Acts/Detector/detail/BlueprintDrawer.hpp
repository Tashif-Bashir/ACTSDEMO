// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Utilities/BinningData.hpp"

namespace Acts {

namespace Experimental {

namespace detail {
namespace BlueprintDrawer {

/// @brief Nested options struct for the drawer
struct Options {
  struct Node {
    /// ROOT node definitions
    std::string shape = "circle";
    std::string color = "darkorange";

    /// Font properties
    std::string face = "sans-serif";
    int labelText = 12;
    int infoText = 10;

    /// Info properties
    int precision = 1;
  };

  /// @brief The name of the graph
  std::string graphName = "blueprint";

  // Main node types to come
  Node root = Node{};
  Node branch = Node{"diamond", "white"};
  Node leaf = Node{"box", "darkolivegreen1"};
  Node gap = Node{"box", "darkolivegreen3"};

  // Sub node types to come
  Node shape = Node{"cylinder", "lightgrey"};
  Node virtualShape = Node{"cylinder", "white"};
  Node internals = Node{"doubleoctagon", "cadetblue1"};
  Node geoID = Node{"box", "azure"};
  Node roots = Node{"box", "darkkhaki"};
};

/// @brief Geneate the shape string
/// @param s the shape of the object
/// @param c the color of the object
/// @return a string with the shape and color
std::string shapeStr(const Options::Node& node) {
  return "[shape=\"" + node.shape + "\";style=\"filled\";fillcolor=\"" +
         node.color + "\"];";
}

/// @brief Generate text output
///
/// @param node the node otions
/// @param label the label text
/// @param info the info text
std::string labelStr(const Options::Node& node, const std::string& label,
                     const std::vector<std::string>& info = {}) {
  std::string lText = "[label=<<font face=\"";
  lText += node.face;
  lText += "\" point-size=\"";
  lText += std::to_string(node.labelText);
  lText += "\">" + label;
  if (!info.empty()) {
    lText += "</font><br/>";
    lText += "<font face=\"";
    lText += node.face;
    lText += "\" point-size=\"";
    lText += std::to_string(node.infoText);
    lText += "\">";
    for (const auto& i : info) {
      lText += i;
      lText += "<br/>";
    }
  }
  lText += "</font>";
  lText += ">];";
  return lText;
}

/// @brief  Turn into a dot output
///
/// @tparam the stream type to be used
///
/// @param ss the stream into which the dot output should be written
/// @param node the node to be drawn
/// @param options the options for the drawer

template <typename stream_type>
void dotStream(stream_type& ss, const Blueprint::Node& node,
               const Options& options = Options{}) {
  // Root / leaf or branch
  if (node.isRoot()) {
    ss << "digraph " << options.graphName << " {" << '\n';
    ss << node.name << " " << labelStr(options.root, node.name, node.auxiliary)
       << '\n';
    ss << node.name << " " << shapeStr(options.root) << '\n';

  } else if (node.isLeaf()) {
    ss << node.name << " " << labelStr(options.leaf, node.name, node.auxiliary)
       << '\n';
    ss << node.name << " "
       << ((node.internalsBuilder != nullptr) ? shapeStr(options.leaf)
                                              : shapeStr(options.gap))
       << '\n';
  } else {
    ss << node.name << " "
       << labelStr(options.branch, node.name, node.auxiliary) << '\n';
    ss << node.name << " " << shapeStr(options.branch) << '\n';
  }
  // Recursive for children
  for (const auto& c : node.children) {
    ss << node.name << " -> " << c->name << ";" << '\n';
    dotStream(ss, *c, options);
  }

  // Shape
  Options::Node shape = node.isLeaf() ? options.shape : options.virtualShape;
  ss << node.name + "_shape " << shapeStr(shape) << '\n';
  ss << node.name + "_shape "
     << labelStr(shape, VolumeBounds::s_boundsTypeNames[node.boundsType],
                 {"t = " + toString(node.transform.translation(), 1),
                  "b = " + toString(node.boundaryValues, 1)})
     << '\n';
  ss << node.name << " -> " << node.name + "_shape [ arrowhead = \"none\" ];"
     << '\n';

  // Sub node dection
  if (node.internalsBuilder != nullptr) {
    ss << node.name + "_int " << shapeStr(options.internals) << '\n';
    ss << node.name << " -> " << node.name + "_int;" << '\n';
  }

  if (node.geoIdGenerator != nullptr) {
    ss << node.name + "_geoID " << shapeStr(options.geoID) << '\n';
    ss << node.name << " -> " << node.name + "_geoID;" << '\n';
  }

  if (node.rootVolumeFinderBuilder != nullptr) {
    ss << node.name + "_roots " << shapeStr(options.roots) << '\n';
    ss << node.name << " -> " << node.name + "_roots;" << '\n';
  }

  if (node.isRoot()) {
    ss << "}" << '\n';
  }
}

}  // namespace BlueprintDrawer
}  // namespace detail
}  // namespace Experimental
}  // namespace Acts
