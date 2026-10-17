// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/main_core/topology.hpp>

#include <power_grid_model/auxiliary/input.hpp>
#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/component_list.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/component/appliance.hpp>
#include <power_grid_model/component/base.hpp>
#include <power_grid_model/component/branch.hpp>
#include <power_grid_model/component/edge.hpp>
#include <power_grid_model/component/line.hpp>
#include <power_grid_model/component/link.hpp>
#include <power_grid_model/component/load_gen.hpp>
#include <power_grid_model/component/node.hpp>
#include <power_grid_model/container.hpp>
#include <power_grid_model/main_core/container_queries.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <vector>

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

using LoadGenContainer =
    Container<ExtraRetrievableTypes<Base, Node, Appliance, GenericLoadGen>, Node, SymLoad, AsymLoad, AsymGenerator>;

constexpr double dummy_value = 1.0;

// a sym_load, an asym_load and an asym_gen at one node; the asym_gen gets the given connection
LoadGenContainer make_load_gens(LoadGenConnection gen_connection) {
    LoadGenContainer components;
    emplace_component<Node>(components, 1, NodeInput{.id = 1, .u_rated = 10e3});
    emplace_component<SymLoad>(components, 2,
                               SymLoadGenInput{.id = 2, .node = 1, .status = 1, .type = LoadGenType::const_pq}, 10e3);
    emplace_component<AsymLoad>(
        components, 3,
        AsymLoadGenInput{
            .id = 3, .node = 1, .status = 1, .type = LoadGenType::const_i, .connection = LoadGenConnection::wye},
        10e3);
    emplace_component<AsymGenerator>(
        components, 4,
        AsymLoadGenInput{.id = 4, .node = 1, .status = 0, .type = LoadGenType::const_y, .connection = gen_connection},
        10e3);
    components.set_construction_complete();
    return components;
}

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
        CHECK(conn.branch_id.empty());
    }

    SUBCASE("Previously filled incidence is cleared when all branches are abc") {
        FullPhaseContainer components;
        add_nodes_and_line(components);
        components.set_construction_complete();

        ComponentConnections conn;
        conn.branch_phase_incidence.assign(3, three_phase_winding);
        conn.branch_id.assign(3, ID{1});
        detail::register_connections_components<Branch>(components, conn);
        CHECK(conn.branch_phase_incidence.empty());
        CHECK(conn.branch_id.empty());
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
        // the IDs are registered together with the incidence, for the error messages of the topology
        CHECK(conn.branch_id == std::vector<ID>{3, 5});
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
        auto const link_idx = static_cast<size_t>(get_component_sequence_idx<Edge>(components, ID{4}));
        CHECK(link_idx != expected_idx);

        REQUIRE(conn.branch_connected.size() == 3);
        REQUIRE(conn.branch_phase_incidence.size() == conn.branch_connected.size());
        REQUIRE(conn.branch_id.size() == conn.branch_connected.size());
        CHECK(conn.branch_id[expected_idx] == 5);
        CHECK(conn.branch_id[link_idx] == 4);
        REQUIRE(expected_idx < conn.branch_connected.size());
        CHECK(conn.branch_connected[expected_idx] == BranchConnected{1, 0});
        for (size_t idx = 0; idx < conn.branch_phase_incidence.size(); ++idx) {
            if (idx == expected_idx) {
                CHECK(conn.branch_phase_incidence[idx] == BranchPhaseIncidence{.from = phase_mask_a,
                                                                               .to = phase_mask_a,
                                                                               .coupling = PhaseCoupling::per_phase});
            } else if (idx == link_idx) {
                // on the legacy path the link is an Edge, which reports its own coupling
                CHECK(conn.branch_connected[idx] == BranchConnected{1, 1});
                CHECK(conn.branch_phase_incidence[idx] == link_incidence);
            } else {
                CHECK(conn.branch_connected[idx] == BranchConnected{1, 1});
                CHECK(conn.branch_phase_incidence[idx] == BranchPhaseIncidence{});
            }
        }
    }
}

TEST_CASE("Test main core topology - load_gen connection registration") {
    ComponentTopology comp_topo;

    SUBCASE("All load_gens wye: the connection vector stays empty") {
        auto const components = make_load_gens(LoadGenConnection::wye);
        detail::register_topology_components<GenericLoadGen>(components, comp_topo, false);
        CHECK(comp_topo.load_gen_type.size() == 3);
        CHECK(comp_topo.load_gen_connection.empty());
    }

    SUBCASE("Previously filled connection vector is cleared when all load_gens are wye") {
        auto const components = make_load_gens(LoadGenConnection::wye);
        comp_topo.load_gen_connection.assign(5, LoadGenConnection::delta);
        detail::register_topology_components<GenericLoadGen>(components, comp_topo, false);
        CHECK(comp_topo.load_gen_connection.empty());
    }

    SUBCASE("One delta load_gen fills the connection vector, status ignored") {
        auto const components = make_load_gens(LoadGenConnection::delta);
        detail::register_topology_components<GenericLoadGen>(components, comp_topo, false);
        REQUIRE(comp_topo.load_gen_connection.size() == 3);
        for (ID const id : {2, 3, 4}) {
            CAPTURE(id);
            auto const idx = static_cast<size_t>(get_component_sequence_idx<GenericLoadGen>(components, id));
            CHECK(comp_topo.load_gen_connection[idx] == (id == 4 ? LoadGenConnection::delta : LoadGenConnection::wye));
        }
    }
}

} // namespace power_grid_model::main_core
