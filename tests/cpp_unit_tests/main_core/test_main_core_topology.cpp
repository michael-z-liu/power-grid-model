// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/main_core/topology.hpp>

#include <power_grid_model/auxiliary/input.hpp>
#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/component_list.hpp>
#include <power_grid_model/component/base.hpp>
#include <power_grid_model/component/branch.hpp>
#include <power_grid_model/component/edge.hpp>
#include <power_grid_model/component/line.hpp>
#include <power_grid_model/component/link.hpp>
#include <power_grid_model/component/node.hpp>
#include <power_grid_model/container.hpp>
#include <power_grid_model/main_core/container_queries.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cstddef>
#include <iterator>

namespace power_grid_model::main_core {
namespace {
// Branch that only connects phase a, used to exercise a non-trivial phase incidence
class SinglePhaseTestBranch final : public Branch {
  public:
    explicit SinglePhaseTestBranch(BranchInput const& branch_input) : Branch{branch_input} {}

    double base_i_from() const override { return 1.0; }
    double base_i_to() const override { return 1.0; }
    double loading(double /*max_s*/, double /*max_i*/) const override { return 0.0; }
    double phase_shift() const override { return 0.0; }
    BranchPhaseIncidence phase_incidence() const override {
        return {.from = phase_mask_a, .to = phase_mask_a, .coupling = PhaseCoupling::per_phase};
    }

  private:
    BranchCalcParam<symmetric_t> sym_calc_param() const override { return calc_param_y_sym(1.0, 0.0, 1.0); }
    BranchCalcParam<asymmetric_t> asym_calc_param() const override {
        return calc_param_y_asym(1.0, 0.0, 1.0, 0.0, 1.0);
    }
};

using FullPhaseContainer = Container<ExtraRetrievableTypes<Base, Node, Edge, Branch>, Line, Link, Node>;
using MixedPhaseContainer =
    Container<ExtraRetrievableTypes<Base, Node, Edge, Branch>, Line, Link, Node, SinglePhaseTestBranch>;

constexpr double dummy_value = 1.0;

template <class Container> void add_nodes_and_line(Container& components) {
    emplace_component<Node>(components, 1, NodeInput{.id = 1, .u_rated = 10e3});
    emplace_component<Node>(components, 2, NodeInput{.id = 2, .u_rated = 10e3});
    emplace_component<Line>(components, 3,
                            LineInput{.id = 3, .from_node = 1, .to_node = 2, .from_status = 1, .to_status = 1},
                            dummy_value, 10e3, 10e3);
    emplace_component<Link>(
        components, 4, LinkInput{.id = 4, .from_node = 1, .to_node = 2, .from_status = 1, .to_status = 1}, 10e3, 10e3);
}
} // namespace

TEST_CASE("Test main core topology - branch phase incidence registration") {
    SUBCASE("All branches full abc: incidence stays empty") {
        FullPhaseContainer components;
        add_nodes_and_line(components);
        components.set_construction_complete();

        ComponentConnections conn;
        SUBCASE("Branch (v2) path") {
            detail::register_connections_components<Branch>(components, conn);
            CHECK(conn.branch_connected.size() == 1);
            CHECK(conn.branch_phase_shift.size() == 1);
        }
        SUBCASE("Edge (legacy) path") {
            detail::register_connections_components<Edge>(components, conn);
            CHECK(conn.branch_connected.size() == 2);
            CHECK(conn.branch_phase_shift.size() == 2);
        }
        CHECK(conn.branch_phase_incidence.empty());
    }

    SUBCASE("Previously filled incidence is cleared when all branches are abc") {
        FullPhaseContainer components;
        add_nodes_and_line(components);
        components.set_construction_complete();

        ComponentConnections conn;
        conn.branch_phase_incidence.assign(3, three_phase_winding);
        detail::register_connections_components<Branch>(components, conn);
        CHECK(conn.branch_phase_incidence.empty());
    }

    SUBCASE("A reduced-phase branch fills the incidence vector") {
        MixedPhaseContainer components;
        add_nodes_and_line(components);
        emplace_component<SinglePhaseTestBranch>(
            components, 5, BranchInput{.id = 5, .from_node = 1, .to_node = 2, .from_status = 1, .to_status = 1});
        components.set_construction_complete();

        ComponentConnections conn;
        detail::register_connections_components<Branch>(components, conn);
        REQUIRE(conn.branch_connected.size() == 2);
        REQUIRE(conn.branch_phase_incidence.size() == 2);
        CHECK(conn.branch_phase_incidence[0] == BranchPhaseIncidence{});
        CHECK(conn.branch_phase_incidence[1] ==
              BranchPhaseIncidence{.from = phase_mask_a, .to = phase_mask_a, .coupling = PhaseCoupling::per_phase});
    }

    SUBCASE("A reduced-phase branch fills the incidence vector on the legacy path with interleaved links") {
        MixedPhaseContainer components;
        add_nodes_and_line(components);
        // to side switched off, so that the branch is recognizable in branch_connected
        emplace_component<SinglePhaseTestBranch>(
            components, 5, BranchInput{.id = 5, .from_node = 1, .to_node = 2, .from_status = 1, .to_status = 0});
        components.set_construction_complete();

        ComponentConnections conn;
        detail::register_connections_components<Edge>(components, conn);
        // the index of the reduced-phase branch follows the iteration order of all Edges, links included
        auto const edges = components.citer<Edge>();
        auto const edge_it = std::ranges::find_if(edges, [](Edge const& edge) { return edge.id() == 5; });
        REQUIRE(edge_it != edges.end());
        auto const expected_idx = static_cast<size_t>(std::ranges::distance(edges.begin(), edge_it));
        CHECK(static_cast<Idx>(expected_idx) == get_component_sequence_idx<Edge>(components, ID{5}));

        REQUIRE(conn.branch_connected.size() == 3);
        REQUIRE(conn.branch_phase_incidence.size() == conn.branch_connected.size());
        REQUIRE(expected_idx < conn.branch_connected.size());
        CHECK(conn.branch_connected[expected_idx] == BranchConnected{1, 0});
        for (size_t idx = 0; idx < conn.branch_phase_incidence.size(); ++idx) {
            if (idx == expected_idx) {
                CHECK(conn.branch_phase_incidence[idx] == BranchPhaseIncidence{.from = phase_mask_a,
                                                                               .to = phase_mask_a,
                                                                               .coupling = PhaseCoupling::per_phase});
            } else {
                CHECK(conn.branch_connected[idx] == BranchConnected{1, 1});
                CHECK(conn.branch_phase_incidence[idx] == BranchPhaseIncidence{});
            }
        }
    }
}

} // namespace power_grid_model::main_core
