// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/topology.hpp>

#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/exception.hpp>
#include <power_grid_model/common/grouped_index_vector.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/*
 *  [0]   = Node / Bus
 * --0--> = Branch (from --id--> to)
 * -(b0)- = Virtual node, representing Branch3,
 *          the sequence within the branch 3 is indicated by $0 $1 $2
 *  -X-   = Open switch / not connected
 *  s0    = source
 *  s0X   = disconnected source
 *  lg0   = load_gen
 *  h0    = shunt
 *  +v0   = voltage sensor
 *  +p0   = power sensor
 *  +p?0  = components power sensor (e.g. +ps1 = source power sensor 1)
 *          (f = branch_from, t = branch_to, s = source, h = shunt, l = load, g = generator, b = bus)
 *
 * Topology:
 *
 *                                               7 -> [5+v4+p17:          ]
 *                                              /     [s1+p1+p12,lg2+p4+p8]     [6:h1+p5+p9] -X-4-> [7] -3-> [8+v5]
 *      0 ----->+p13 [1+v1:lg3+p7]             /     /                     \        /                  \        /
 *     /          +p14\            5 ---X--- [4] <- 6                      $2     $1                   $1     $2
 *    /                $0         /          ^  \                            \    /                      \    /
 *   /                  \        v          /    $0                           (b2)                        (b1)
 *  [0:s0,lg0]         (b0)-$2- [2+v0+v2]  /      \                            |                           |
 *   +p0+p11            /    +p15         X       (b3)                         $0                          $0
 *    \                $1                /         \  \                        X                           |
 *     \          +p16/                 /           $1 $2                      [10]                        [11:lg1+p6]
 *      1 -->+p2+p10 [3+v3:s3X,h0] -- 2              \  \
 *                    |         |                     [12]
 *                    |         |                                                                 /--$0--\
 *                    \----8----/                                   [9:s2X+p3,h2]          [13:s3]---$1---(b4)
 *                                                                                                \--$2--/
 *
 *
 * Math model #0:                       Math model #1:
 *
 *      0 ----->+pt0 [2+v3:lg0+pg0]             1 -> [5+v0+pb0:s0+ps0+ps1,lg0+pl0+pl1]         [0:h1+ps0+ps1]
 *     /          +pf2\             3 --X      /     /                                \       /
 *    /                5           /         [4] <- 0                                  4     3
 *   /                  v         v            \                                        v   v
 *  [4:s0,lg1]          [3] <-7- [0+v0+v1]      5                                        [1]
 *   +pf0+pf1           ^     +pf4               v                                        ^
 *    \                6                         [3]                                      2
 *     \          +pf3/                           ^ ^                                     X
 *      1 ->+pt1+pt2 [1+v2:h0] -- 2 --X            6 7
 *                    |     |                       \ \
 *                    |     |                        [2]
 *                    \--4--/
 *
 *                                      Math model #2:
 *                                                       /--0-->
 *                                                 [1:s0]---1-->[0]
 *                                                       \--2-->
 *
 * Extra fill-in:
 * (3, 4)  by removing node 1
 *
 *
 * Topology for cycle reodering
 *
 *
 *   [5]  <---4--[4] <--3- [3]
 *    ^ \         ^       /  ^
 *    |   9----   |     /    |
 *    5        \  6   10     2
 *    |         v |  v       |
 * [0:s0] --0--> [1] --1--> [2]
 *    ^        ^    <- 12-   ^
 *    |   -11-/     parallel |
 *    7  /                   |
 *    | /                    |
 *   [6] -----------------8--
 *
 * Math model after reodering
 *
 *   [1]  <---4--[6] <--3- [2]
 *    ^ \         ^       /  ^
 *    |   9----   |     /    |
 *    5        \  6   10     2
 *    |         v |  v       |
 * [3:s0] --0--> [5] --1--> [4]
 *    ^        ^    <- 12-   ^
 *    |   -11-/     parallel |
 *    7  /                   |
 *    | /                    |
 *   [0] -----------------8--
 *
 * Extra fill-in:
 * (3, 4)  by removing node 0
 * (3, 6)  by removing node 1
 * (4, 6)  by removing node 2
 */

namespace power_grid_model::topology {

namespace {

template <grouped_idx_vector_type T> void check_equal(T const& first, T const& second) {
    REQUIRE(first.size() == second.size());
    for ([[maybe_unused]] auto const& [index, first_element, second_element] : enumerated_zip_sequence(first, second)) {
        CHECK(std::ranges::equal(first_element, second_element));
    }
}

// phase mask and branch phase incidence
static_assert(PhaseMask{}.none());
static_assert(!PhaseMask{}.all());
static_assert(phase_mask_abc.all());
static_assert(!phase_mask_abc.none());
static_assert(PhaseMask{0b1111} == phase_mask_abc); // bits beyond c are dropped
static_assert(phase_mask_a.bits() == 0b001 && phase_mask_b.bits() == 0b010 && phase_mask_c.bits() == 0b100);
static_assert(phase_mask_ab.has(0) && phase_mask_ab.has(1) && !phase_mask_ab.has(2));
static_assert(phase_mask_ac.has(0) && !phase_mask_ac.has(1) && phase_mask_ac.has(2));
static_assert(!phase_mask_bc.has(0) && phase_mask_bc.has(1) && phase_mask_bc.has(2));
static_assert(!phase_mask_abc.has(-1) && !phase_mask_abc.has(3));
static_assert(!phase_mask_none.has(0) && !phase_mask_none.has(1) && !phase_mask_none.has(2));
static_assert(phase_mask_abc.contains(phase_mask_ab) && phase_mask_ab.contains(phase_mask_a));
static_assert(!phase_mask_ab.contains(phase_mask_c) && !phase_mask_a.contains(phase_mask_ab));
static_assert(phase_mask_a.contains(phase_mask_none) && phase_mask_none.contains(phase_mask_none));
static_assert((phase_mask_a | phase_mask_b) == phase_mask_ab);
static_assert((phase_mask_ab | phase_mask_bc) == phase_mask_abc);
static_assert((phase_mask_ab & phase_mask_bc) == phase_mask_b);
static_assert((phase_mask_a & phase_mask_b) == phase_mask_none);
static_assert(phase_letters(phase_mask_none) == "-");
static_assert(phase_letters(phase_mask_a) == "a");
static_assert(phase_letters(phase_mask_b) == "b");
static_assert(phase_letters(phase_mask_ab) == "ab");
static_assert(phase_letters(phase_mask_c) == "c");
static_assert(phase_letters(phase_mask_ac) == "ac");
static_assert(phase_letters(phase_mask_bc) == "bc");
static_assert(phase_letters(phase_mask_abc) == "abc");
static_assert(BranchPhaseIncidence{}.all_phases());
static_assert(BranchPhaseIncidence{}.coupling == PhaseCoupling::per_phase);
static_assert(three_phase_winding.all_phases());
static_assert(three_phase_winding.coupling == PhaseCoupling::winding);
static_assert(three_phase_winding != BranchPhaseIncidence{});
static_assert(!BranchPhaseIncidence{.from = phase_mask_a, .to = phase_mask_a}.all_phases());
static_assert(link_incidence.all_phases());
static_assert(link_incidence.coupling == PhaseCoupling::link);
static_assert(link_incidence != BranchPhaseIncidence{});
static_assert(phase_mask_abc.without(phase_mask_a) == phase_mask_bc);
static_assert(phase_mask_ab.without(phase_mask_bc) == phase_mask_a);
static_assert(phase_mask_a.without(phase_mask_abc) == phase_mask_none);
static_assert(phase_mask_ac.without(phase_mask_none) == phase_mask_ac);
static_assert(phase_mask_none.without(phase_mask_b) == phase_mask_none);

using TopologyResult = std::pair<std::vector<std::shared_ptr<MathModelTopology const>>,
                                 std::shared_ptr<TopologicalComponentToMathCoupling const>>;

// per phase branch with the same phases at both terminals
constexpr BranchPhaseIncidence pp(PhaseMask mask) {
    return {.from = mask, .to = mask, .coupling = PhaseCoupling::per_phase};
}
// winding branch
constexpr BranchPhaseIncidence winding(PhaseMask from, PhaseMask to) {
    return {.from = from, .to = to, .coupling = PhaseCoupling::winding};
}
// winding branch with a centre-tapped (split-phase) to side
constexpr BranchPhaseIncidence winding_ct(PhaseMask from) {
    return {
        .from = from, .to = phase_mask_ab, .coupling = PhaseCoupling::winding, .to_pattern = SlotPattern::split_phase};
}

struct PhaseTestBranch {
    BranchIdx nodes;
    BranchPhaseIncidence incidence{};
    BranchConnected connected{1, 1};
    double phase_shift{}; // theta_from - theta_to
};

// small topology fixture with the phase incidence of each branch set explicitly
struct PhaseTestNetwork {
    ReducedComponentTopology comp_topo{};
    ComponentConnections comp_conn{};

    PhaseTestNetwork(Idx n_node, std::vector<Idx> sources, std::vector<PhaseTestBranch> const& branches) {
        comp_topo.n_node = n_node;
        comp_topo.source_node_idx = std::move(sources);
        comp_conn.source_connected = std::vector<IntS>(comp_topo.source_node_idx.size(), IntS{1});
        std::vector<BranchIdx> branch_node_idx;
        for (auto const& branch : branches) {
            branch_node_idx.push_back(branch.nodes);
            comp_conn.branch_connected.push_back(branch.connected);
            comp_conn.branch_phase_shift.push_back(branch.phase_shift);
            comp_conn.branch_phase_incidence.push_back(branch.incidence);
        }
        comp_topo.branch_node_idx = std::move(branch_node_idx);
    }

    void add_branch3(Branch3Idx nodes, Branch3Connected connected = {1, 1, 1}) {
        std::vector<Branch3Idx> branch3_node_idx(comp_topo.branch3_node_idx.begin(), comp_topo.branch3_node_idx.end());
        branch3_node_idx.push_back(nodes);
        comp_topo.branch3_node_idx = std::move(branch3_node_idx);
        comp_conn.branch3_connected.push_back(connected);
        comp_conn.branch3_phase_shift.push_back({0.0, 0.0, 0.0});
    }

    TopologyResult build() const { return Topology{comp_topo, comp_conn}.build_topology(); }
};

PhaseMask node_phase(TopologyResult const& result, Idx node) {
    auto const [group, pos] = result.second->node[static_cast<size_t>(node)];
    return result.first[static_cast<size_t>(group)]->bus_phase(pos);
}

std::vector<std::string> node_phase_letters(TopologyResult const& result, std::vector<Idx> const& nodes) {
    std::vector<std::string> letters;
    letters.reserve(nodes.size());
    for (Idx const node : nodes) {
        letters.emplace_back(phase_letters(node_phase(result, node)));
    }
    return letters;
}

using Letters = std::vector<std::string>;

SlotPattern node_slot_pattern(TopologyResult const& result, Idx node) {
    auto const [group, pos] = result.second->node[static_cast<size_t>(node)];
    return result.first[static_cast<size_t>(group)]->slot_pattern(pos);
}

double node_phase_shift(TopologyResult const& result, Idx node) {
    auto const [group, pos] = result.second->node[static_cast<size_t>(node)];
    return result.first[static_cast<size_t>(group)]->phase_shift[static_cast<size_t>(pos)];
}

} // namespace

TEST_CASE("Test topology") {
    // component topology
    ReducedComponentTopology comp_topo{};
    comp_topo.n_node = 14;

    comp_topo.branch_node_idx = {
        {0, 1}, // 0
        {0, 3}, // 1
        {3, 4}, // 2
        {7, 8}, // 3
        {6, 7}, // 4
        {4, 2}, // 5
        {5, 4}, // 6
        {4, 5}, // 7
        {1, 1}  // 8 (branch into itself)
    };
    comp_topo.branch3_node_idx = {
        {1, 3, 2},   // b0
        {11, 7, 8},  // b1
        {10, 6, 5},  // b2
        {4, 12, 12}, // b3 (branch3 with 2 branches into itself)
        {13, 13, 13} // b4 (branch3 with 3 branches into itself)
    };
    comp_topo.source_node_idx = {0, 5, 9, 3, 13};
    comp_topo.load_gen_node_idx = {0, 11, 5, 1};
    comp_topo.load_gen_type = {LoadGenType::const_pq, LoadGenType::const_pq, LoadGenType::const_i,
                               LoadGenType::const_y};
    comp_topo.shunt_node_idx = {3, 6, 9};
    comp_topo.voltage_sensor_node_idx = {2, 1, 2, 3, 5, 8};
    comp_topo.power_sensor_object_idx = {1, 1, 1, 2, 2, 1, 1, 3, 2, 1, 1, 1, 1, 0, 0, 0, 0, 5};
    comp_topo.power_sensor_terminal_type = {
        MeasuredTerminalType::branch_from, // 0 (branch   1)
        MeasuredTerminalType::source,      // 1 (source   1)
        MeasuredTerminalType::branch_to,   // 2 (branch   1)
        MeasuredTerminalType::source,      // 3 (source   2)
        MeasuredTerminalType::load,        // 4 (load_gen 2)
        MeasuredTerminalType::shunt,       // 5 (shunt    1)
        MeasuredTerminalType::load,        // 6 (load_gen 1)
        MeasuredTerminalType::generator,   // 7 (load_gen 3)
        MeasuredTerminalType::load,        // 8 (load_gen 2)
        MeasuredTerminalType::shunt,       // 9 (shunt    1)
        MeasuredTerminalType::branch_to,   // 10 (branch   1)
        MeasuredTerminalType::branch_from, // 11 (branch   1)
        MeasuredTerminalType::source,      // 12 (source   1)
        MeasuredTerminalType::branch_to,   // 13 (branch   0)
        MeasuredTerminalType::branch3_1,   // 14 (branch3  0)
        MeasuredTerminalType::branch3_3,   // 15 (branch3  0)
        MeasuredTerminalType::branch3_2,   // 16 (branch3  0)
        MeasuredTerminalType::node,        // 17 (node     5)
    };

    // component connection
    ComponentConnections comp_conn{};
    comp_conn.branch_connected = {
        // {from, to}
        {1, 1}, // 0
        {1, 1}, // 1
        {1, 0}, // 2
        {1, 1}, // 3
        {0, 1}, // 4
        {0, 1}, // 5
        {1, 1}, // 6
        {1, 1}, // 7
        {1, 1}, // 8
    };
    comp_conn.branch3_connected = {
        {1, 1, 1}, // b0
        {1, 1, 1}, // b1
        {0, 1, 1}, // b2
        {1, 1, 1}, // b3
        {1, 1, 1}, // b4
    };
    comp_conn.branch_phase_shift = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    comp_conn.branch3_phase_shift = {
        {0.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
    };
    comp_conn.source_connected = {1, 1, 0, 0, 1};

    // result
    TopologicalComponentToMathCoupling comp_coup_ref{};
    comp_coup_ref.node = {
        // 0 1 2 3
        {.group = 0, .pos = 1}, // Topological node 0 has become node 4 in mathematical model (group) 0
        {.group = 0, .pos = 2},
        {.group = 0, .pos = 0},
        {.group = 0, .pos = 4},
        // 4 5 6
        {.group = 1, .pos = 4}, // Topological node 4 has become node 4 in mathematical model (group) 1
        {.group = 1, .pos = 5},
        {.group = 1, .pos = 0},
        // 7, 8, 9, 10, 11, 13
        {.group = -1, .pos = -1}, // Topological node 7 is not included in the mathematical model, because it was not
                                  // connected to any power source
        {.group = -1, .pos = -1},
        {.group = -1, .pos = -1},
        {.group = -1, .pos = -1},
        {.group = -1, .pos = -1},
        {.group = 1, .pos = 2},
        // 13
        {.group = 2, .pos = 1},
        // b0, b1, b2, b3, b4
        {.group = 0, .pos = 3}, // Branch3 b0 is replaced by a virtual node 3, in mathematical model 0
        {.group = -1, .pos = -1},
        {.group = 1, .pos = 1},
        {.group = 1, .pos = 3},
        {.group = 2, .pos = 0},
    };
    comp_coup_ref.source = {
        {.group = 0, .pos = 0},   // 0
        {.group = 1, .pos = 0},   // 1
        {.group = -1, .pos = -1}, // 2
        {.group = -1, .pos = -1}, // 3
        {.group = 2, .pos = 0},   // 4
    };
    comp_coup_ref.branch = {
        {.group = 0, .pos = 0},   // 0
        {.group = 0, .pos = 1},   // 1
        {.group = 0, .pos = 2},   // 2
        {.group = -1, .pos = -1}, // 3
        {.group = -1, .pos = -1}, // 4
        {.group = 0, .pos = 3},   // 5
        {.group = 1, .pos = 0},   // 6
        {.group = 1, .pos = 1},   // 7
        {.group = 0, .pos = 4},   // 8
    };
    comp_coup_ref.branch3 = {
        {.group = 0, .pos = {5, 6, 7}},     // b0
        {.group = -1, .pos = {-1, -1, -1}}, // b1
        {.group = 1, .pos = {2, 3, 4}},     // b2
        {.group = 1, .pos = {5, 6, 7}},     // b3
        {.group = 2, .pos = {0, 1, 2}},     // b4
    };
    comp_coup_ref.load_gen = {
        {.group = 0, .pos = 0}, {.group = -1, .pos = -1}, {.group = 1, .pos = 0}, {.group = 0, .pos = 1}};
    comp_coup_ref.shunt = {{.group = 0, .pos = 0}, {.group = 1, .pos = 0}, {.group = -1, .pos = -1}};
    comp_coup_ref.voltage_sensor = {{.group = 0, .pos = 0}, {.group = 0, .pos = 2}, {.group = 0, .pos = 1},
                                    {.group = 0, .pos = 3}, {.group = 1, .pos = 0}, {.group = -1, .pos = -1}};
    comp_coup_ref.power_sensor = {
        {.group = 0, .pos = 0},   // 0 branch_from
        {.group = 1, .pos = 0},   // 1 source
        {.group = 0, .pos = 1},   // 2 branch_to
        {.group = -1, .pos = -1}, // 3 source
        {.group = 1, .pos = 0},   // 4 load       = load power sensor 0 in math model 1
        {.group = 1, .pos = 0},   // 5 shunt      = shunt power sensor 0 in math model 1
        {.group = -1, .pos = -1}, // 6 load
        {.group = 0, .pos = 0},   // 7 generator
        {.group = 1, .pos = 1},   // 8 load
        {.group = 1, .pos = 1},   // 9 shunt
        {.group = 0, .pos = 2},   // 10 branch_to
        {.group = 0, .pos = 1},   // 11 branch_from
        {.group = 1, .pos = 1},   // 12 source
        {.group = 0, .pos = 0},   // 13 branch_to
        {.group = 0, .pos = 2},   // 14 branch_from
        {.group = 0, .pos = 4},   // 15 branch_from
        {.group = 0, .pos = 3},   // 16 branch_from
        {.group = 1, .pos = 0}    // 17 node
    };

    // Sub graph / math model 0
    MathModelTopology math0;
    math0.slack_bus = 1;
    math0.sources_per_bus = {from_dense, {1}, 5};
    math0.branch_bus_idx = {{1, 2}, {1, 4}, {4, -1}, {-1, 0}, {2, 2}, {2, 3}, {4, 3}, {0, 3}};
    math0.phase_shift = {0.0, 0.0, 0.0, 0.0, -1.0};
    math0.load_gens_per_bus = {from_dense, {1, 2}, 5};
    math0.load_gen_type = {LoadGenType::const_pq, LoadGenType::const_y};
    math0.shunts_per_bus = {from_dense, {4}, 5};
    math0.voltage_sensors_per_bus = {from_dense, {0, 0, 2, 4}, 5};
    math0.power_sensors_per_bus = {from_dense, {}, 5};
    math0.power_sensors_per_source = {from_dense, {}, 1};
    math0.power_sensors_per_shunt = {from_dense, {}, 1};
    math0.power_sensors_per_load_gen = {from_dense, {1}, 2};
    math0.power_sensors_per_branch_from = {from_dense, {1, 1, 5, 6, 7}, 8};
    // 7 branches, 3 branch-to power sensors
    // sensor 0 is connected to branch 0
    // sensor 1 and 2 are connected to branch 1
    math0.power_sensors_per_branch_to = {from_dense, {0, 1, 1}, 8};
    math0.fill_in = {{2, 4}};

    // Sub graph / math model 1
    MathModelTopology math1;
    math1.slack_bus = 5;
    math1.sources_per_bus = {from_dense, {5}, 6};
    math1.branch_bus_idx = {{5, 4}, {4, 5}, {-1, 1}, {0, 1}, {5, 1}, {4, 3}, {2, 3}, {2, 3}};
    math1.phase_shift = {0, 0, 0, 0, 0, 0};
    math1.load_gens_per_bus = {from_dense, {5}, 6};
    math1.load_gen_type = {LoadGenType::const_i};
    math1.shunts_per_bus = {from_dense, {0}, 6};
    math1.voltage_sensors_per_bus = {from_dense, {5}, 6};
    math1.power_sensors_per_bus = {from_dense, {5}, 6};
    math1.power_sensors_per_source = {from_dense, {0, 0}, 1};
    math1.power_sensors_per_shunt = {from_dense, {0, 0}, 1};
    math1.power_sensors_per_load_gen = {from_dense, {0, 0}, 1};
    math1.power_sensors_per_branch_from = {from_dense, {}, 8};
    math1.power_sensors_per_branch_to = {from_dense, {}, 8};

    // Sub graph / math model 2
    MathModelTopology math2;
    math2.slack_bus = 1;
    math2.sources_per_bus = {from_dense, {1}, 2};
    math2.branch_bus_idx = {{1, 0}, {1, 0}, {1, 0}};
    math2.phase_shift = {0, 0};
    math2.load_gens_per_bus = {from_dense, {}, 2};
    math2.load_gen_type = {};
    math2.shunts_per_bus = {from_dense, {}, 2};
    math2.voltage_sensors_per_bus = {from_dense, {}, 2};
    math2.power_sensors_per_bus = {from_dense, {}, 2};
    math2.power_sensors_per_source = {from_dense, {}, 1};
    math2.power_sensors_per_shunt = {from_dense, {}, 0};
    math2.power_sensors_per_load_gen = {from_dense, {}, 0};
    math2.power_sensors_per_branch_from = {from_dense, {}, 3};
    math2.power_sensors_per_branch_to = {from_dense, {}, 3};

    std::vector<MathModelTopology> math_topology_ref = {math0, math1, math2};

    SUBCASE("Test topology result") {
        Topology topo{comp_topo, comp_conn};
        auto const [math_topology, topo_comp_coup_ptr] = topo.build_topology();

        REQUIRE(topo_comp_coup_ptr != nullptr);
        auto const& topo_comp_coup = *topo_comp_coup_ptr;

        CHECK(math_topology.size() == math_topology_ref.size());
        // test component coupling
        CHECK(topo_comp_coup.node == comp_coup_ref.node);
        CHECK(topo_comp_coup.source == comp_coup_ref.source);
        CHECK(topo_comp_coup.branch == comp_coup_ref.branch);
        CHECK(topo_comp_coup.branch3 == comp_coup_ref.branch3);
        CHECK(topo_comp_coup.load_gen == comp_coup_ref.load_gen);
        CHECK(topo_comp_coup.voltage_sensor == comp_coup_ref.voltage_sensor);
        CHECK(topo_comp_coup.power_sensor == comp_coup_ref.power_sensor);

        for (size_t i = 0; i < math_topology.size(); i++) {
            auto const& math = *math_topology[i];
            auto const& math_ref = math_topology_ref[i];
            CHECK(math.slack_bus == math_ref.slack_bus);
            CHECK(math.n_bus() == math_ref.n_bus());
            CHECK_FALSE(math.has_reduced_phase);
            CHECK(math.bus_phase_mask.empty());
            CHECK(math.bus_slot_pattern.empty());
            check_equal(math.sources_per_bus, math_ref.sources_per_bus);
            CHECK(math.branch_bus_idx == math_ref.branch_bus_idx);
            CHECK(math.phase_shift == math_ref.phase_shift);
            check_equal(math.load_gens_per_bus, math_ref.load_gens_per_bus);
            CHECK(math.load_gen_type == math_ref.load_gen_type);
            CHECK(math.load_gen_connection == math_ref.load_gen_connection);
            CHECK_FALSE(math.has_delta_load_gen());
            check_equal(math.shunts_per_bus, math_ref.shunts_per_bus);
            check_equal(math.voltage_sensors_per_bus, math_ref.voltage_sensors_per_bus);
            check_equal(math.power_sensors_per_bus, math_ref.power_sensors_per_bus);
            check_equal(math.power_sensors_per_source, math_ref.power_sensors_per_source);
            check_equal(math.power_sensors_per_shunt, math_ref.power_sensors_per_shunt);
            check_equal(math.power_sensors_per_load_gen, math_ref.power_sensors_per_load_gen);
            check_equal(math.power_sensors_per_branch_from, math_ref.power_sensors_per_branch_from);
            check_equal(math.power_sensors_per_branch_to, math_ref.power_sensors_per_branch_to);
            CHECK(math.fill_in == math_ref.fill_in);
        }
    }

    SUBCASE("Test delta load_gen connection") {
        using enum LoadGenConnection;
        // load_gen 0 -> math model 0 position 0, 1 -> disconnected, 2 -> math model 1 position 0, 3 -> math model 0
        // position 1
        SUBCASE("A connected delta load_gen fills the vector of its math model only") {
            comp_topo.load_gen_connection = std::vector<LoadGenConnection>{wye, delta, wye, delta};
            math_topology_ref[0].load_gen_connection = {wye, delta};
        }
        SUBCASE("A disconnected delta load_gen is skipped") {
            comp_topo.load_gen_connection = std::vector<LoadGenConnection>{wye, delta, wye, wye};
        }

        Topology topo{comp_topo, comp_conn};
        auto const [math_topology, topo_comp_coup_ptr] = topo.build_topology();
        REQUIRE(math_topology.size() == math_topology_ref.size());
        for (size_t i = 0; i < math_topology.size(); i++) {
            CAPTURE(i);
            auto const& math = *math_topology[i];
            CHECK(math.load_gen_type == math_topology_ref[i].load_gen_type);
            CHECK(math.load_gen_connection == math_topology_ref[i].load_gen_connection);
            for (Idx load_gen = 0; load_gen != math.n_load_gen(); ++load_gen) {
                CAPTURE(load_gen);
                CHECK(math.is_delta_load_gen(load_gen) == (i == 0 && load_gen == 1 && math.has_delta_load_gen()));
            }
        }
        CHECK_FALSE(math_topology[1]->has_delta_load_gen());
        CHECK_FALSE(math_topology[2]->has_delta_load_gen());
    }
}

TEST_CASE("Test topology phase mask") {
    constexpr Idx n_node = 4;

    SUBCASE("Explicit incidence of full phase branches gives the same topology") {
        // same network as in "Test topology", with a transformer-like branch 1
        // the result must be identical to that of the reference, which has no incidence
        ReducedComponentTopology comp_topo{};
        comp_topo.n_node = 14;
        comp_topo.branch_node_idx = {{0, 1}, {0, 3}, {3, 4}, {7, 8}, {6, 7}, {4, 2}, {5, 4}, {4, 5}, {1, 1}};
        comp_topo.branch3_node_idx = {{1, 3, 2}, {11, 7, 8}, {10, 6, 5}, {4, 12, 12}, {13, 13, 13}};
        comp_topo.source_node_idx = {0, 5, 9, 3, 13};

        ComponentConnections comp_conn_ref{};
        comp_conn_ref.branch_connected = {{1, 1}, {1, 1}, {1, 0}, {1, 1}, {0, 1}, {0, 1}, {1, 1}, {1, 1}, {1, 1}};
        comp_conn_ref.branch3_connected = {{1, 1, 1}, {1, 1, 1}, {0, 1, 1}, {1, 1, 1}, {1, 1, 1}};
        comp_conn_ref.branch_phase_shift = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        comp_conn_ref.branch3_phase_shift = {
            {0.0, -1.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
        };
        comp_conn_ref.source_connected = {1, 1, 0, 0, 1};

        ComponentConnections comp_conn = comp_conn_ref;
        comp_conn.branch_phase_incidence = std::vector<BranchPhaseIncidence>(9, BranchPhaseIncidence{});
        comp_conn.branch_phase_incidence[1] = three_phase_winding;

        auto const [math_ref, coup_ref] = Topology{comp_topo, comp_conn_ref}.build_topology();
        auto const [math, coup] = Topology{comp_topo, comp_conn}.build_topology();

        CHECK(coup->node == coup_ref->node);
        CHECK(coup->source == coup_ref->source);
        CHECK(coup->branch == coup_ref->branch);
        CHECK(coup->branch3 == coup_ref->branch3);
        REQUIRE(math.size() == math_ref.size());
        for (size_t i = 0; i < math.size(); i++) {
            CHECK_FALSE(math[i]->has_reduced_phase);
            CHECK(math[i]->bus_phase_mask.empty());
            CHECK(math[i]->bus_slot_pattern.empty());
            CHECK(math[i]->slack_bus == math_ref[i]->slack_bus);
            CHECK(math[i]->n_bus() == math_ref[i]->n_bus());
            CHECK(math[i]->phase_shift == math_ref[i]->phase_shift);
            CHECK(math[i]->branch_bus_idx == math_ref[i]->branch_bus_idx);
            CHECK(math[i]->fill_in == math_ref[i]->fill_in);
        }
    }

    SUBCASE("Empty incidence leaves all math models at full phases") {
        PhaseTestNetwork network{3, {0}, {{.nodes = {0, 1}}, {.nodes = {1, 2}}}};
        network.comp_conn.branch_phase_incidence.clear();
        auto const result = network.build();
        REQUIRE(result.first.size() == 1);
        CHECK_FALSE(result.first[0]->has_reduced_phase);
        CHECK(result.first[0]->bus_phase_mask.empty());
        CHECK(result.first[0]->bus_slot_pattern.empty());
        CHECK(node_phase(result, 2) == phase_mask_abc);
        CHECK(node_slot_pattern(result, 2) == SlotPattern::three_phase);
    }

    SUBCASE("Per phase branches restrict the phases of the far bus") {
        // 0 -a- 1 -a- 2
        //  \-abc- 3
        PhaseTestNetwork const network{n_node,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {0, 3}}}};
        auto const result = network.build();
        REQUIRE(result.first.size() == 1);
        CHECK(result.first[0]->has_reduced_phase);
        CHECK(result.first[0]->bus_phase_mask.size() == 4);
        CHECK(result.first[0]->bus_slot_pattern == std::vector<SlotPattern>(4, SlotPattern::three_phase));
        CHECK(node_phase_letters(result, {0, 1, 2, 3}) == Letters{"abc", "a", "a", "abc"});
    }

    SUBCASE("Parallel paths add up their phases") {
        PhaseTestNetwork const network{3,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {0, 2}},
                                        {.nodes = {2, 1}, .incidence = pp(phase_mask_b)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "ab", "abc"});
    }

    SUBCASE("Windings transform the phases") {
        PhaseTestNetwork const network{5,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = winding(phase_mask_ab, phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {2, 3}, .incidence = winding(phase_mask_a, phase_mask_ab)},
                                        {.nodes = {3, 4}, .incidence = pp(phase_mask_b)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2, 3, 4}) == Letters{"abc", "a", "a", "ab", "b"});
    }

    SUBCASE("Windings are directional in their branch direction") {
        PhaseTestNetwork const network{2, {0}, {{.nodes = {1, 0}, .incidence = winding(phase_mask_a, phase_mask_ab)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "a"});
    }

    SUBCASE("The result does not depend on the order of the branches") {
        // in the reversed order node 1 is popped with mask a before it has reached ab, which makes this test
        // meaningful; this relies on the stable placement of the edges in the CSR graph of Boost
        std::vector<PhaseTestBranch> branches{{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                              {.nodes = {1, 3}, .incidence = winding(phase_mask_ab, phase_mask_a)},
                                              {.nodes = {0, 2}},
                                              {.nodes = {2, 1}, .incidence = pp(phase_mask_b)}};
        PhaseTestNetwork const forward{4, {0}, branches};
        CHECK(node_phase_letters(forward.build(), {0, 1, 2, 3}) == Letters{"abc", "ab", "abc", "a"});

        std::ranges::reverse(branches);
        PhaseTestNetwork const reversed{4, {0}, branches};
        CHECK(node_phase_letters(reversed.build(), {0, 1, 2, 3}) == Letters{"abc", "ab", "abc", "a"});
    }

    SUBCASE("A winding without any of its phases present throws") {
        PhaseTestNetwork const network{3,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = winding(phase_mask_bc, phase_mask_a)}}};
        CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), MissingPhaseForWinding);
    }

    SUBCASE("A conductor that is not present at its connected nodes throws") {
        SUBCASE("Branch between two separately fed nodes") {
            // 0 -a- 1 -abc- 2 -a- 0: the b and c conductors of branch 2 are present at neither node
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {0, 2}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #2"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("(bc)"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("from node phases: a, to node phases: a"),
                                 DanglingConductor);
        }
        SUBCASE("Single side stub") {
            PhaseTestNetwork const network{
                3, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 2}, .connected = {1, 0}}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("(bc)"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("to node phases: open"), DanglingConductor);
        }
        SUBCASE("Single side stub, on the to side") {
            PhaseTestNetwork const network{
                3, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {2, 1}, .connected = {0, 1}}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("from node phases: open, to node phases: a"),
                                 DanglingConductor);
        }
        SUBCASE("Branch disjoint from its feeder") {
            // replaces the M1 case of a bus that no phase reaches: such a bus can no longer exist in a math model
            PhaseTestNetwork const network{
                3,
                {0},
                {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 2}, .incidence = pp(phase_mask_b)}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("on phases (b)"), DanglingConductor);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("from node phases: a, to node phases: -"),
                                 DanglingConductor);
        }
        SUBCASE("Self-loop, which is not in the graph but is still validated") {
            PhaseTestNetwork const network{
                2,
                {0},
                {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 1}, .incidence = pp(phase_mask_b)}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), DanglingConductor);
        }
    }

    SUBCASE("A conductor is not validated where it is not needed") {
        SUBCASE("Per phase branch at a node outside every math model") {
            // node 2 has no connected branch and no source, so it is in no math model
            PhaseTestNetwork const network{
                3, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 2}, .connected = {0, 1}}}};
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(result.second->node[2].group == disconnected);
            CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "a"});
        }
        SUBCASE("Link between two separately fed single phase nodes") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {0, 2}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = link_incidence}}};
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "a", "a"});
        }
        SUBCASE("Link feeds the phases of its node") {
            PhaseTestNetwork const network{
                3,
                {0},
                {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 2}, .incidence = link_incidence}}};
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "a", "a"});
        }
    }

    SUBCASE("A winding with only some of its phases present throws") {
        SUBCASE("Three phase winding on a single phase bus") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = three_phase_winding}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), MissingPhaseForWinding);
        }
        SUBCASE("Winding at a side that is the only connected side") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = three_phase_winding, .connected = {1, 0}}}};
            CHECK_THROWS_AS(network.build(), MissingPhaseForWinding);
        }
        SUBCASE("Winding at a side that is the only connected side, on the to side") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {2, 1}, .incidence = three_phase_winding, .connected = {0, 1}}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), MissingPhaseForWinding);
        }
        SUBCASE("Self-loop winding, which is not in the graph but is still validated") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 1}, .incidence = three_phase_winding}}};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch #1"), MissingPhaseForWinding);
        }
        SUBCASE("Partial winding") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = winding(phase_mask_ab, phase_mask_a)}}};
            CHECK_THROWS_AS(network.build(), MissingPhaseForWinding);
        }
        SUBCASE("Branch3 leg") {
            PhaseTestNetwork network{n_node, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}}};
            network.add_branch3({1, 2, 3});
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch3 #0"), MissingPhaseForWinding);
        }
    }

    SUBCASE("Errors name the component ID when the connections carry IDs") {
        SUBCASE("Dangling conductor") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                      {.nodes = {0, 2}, .incidence = pp(phase_mask_a)},
                                      {.nodes = {1, 2}}}};
            network.comp_conn.branch_id = {10, 11, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), DanglingConductor);
        }
        SUBCASE("Branch winding") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                      {.nodes = {1, 2}, .incidence = three_phase_winding}}};
            network.comp_conn.branch_id = {10, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), MissingPhaseForWinding);
        }
        SUBCASE("Branch3 leg") {
            PhaseTestNetwork network{n_node, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}}};
            network.add_branch3({1, 2, 3});
            network.comp_conn.branch_id = {10};
            network.comp_conn.branch3_id = {7};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch3 with ID 7"), MissingPhaseForWinding);
        }
    }

    SUBCASE("A winding at a node outside every math model is not validated") {
        SUBCASE("Branch with only its to side connected") {
            // node 2 has no connected branch and no source, so it is in no math model
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = three_phase_winding, .connected = {0, 1}}}};
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(result.second->node[2].group == disconnected);
            CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "a"});
        }
        SUBCASE("Branch3 with the leg at the reduced bus switched off") {
            // the legs at nodes 2 and 3 are connected to each other but not to a source
            PhaseTestNetwork network{n_node, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}}};
            network.add_branch3({1, 2, 3}, {0, 1, 1});
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(result.second->node[2].group == disconnected);
            CHECK(result.second->node[3].group == disconnected);
            CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "a"});
        }
    }

    SUBCASE("Every source starts with all phases in its own math model") {
        PhaseTestNetwork const network{
            n_node, {0, 2}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {2, 3}}}};
        auto const result = network.build();
        REQUIRE(result.first.size() == 2);
        auto const group0 = static_cast<size_t>(result.second->node[0].group);
        auto const group1 = static_cast<size_t>(result.second->node[2].group);
        CHECK(group0 != group1);
        CHECK(result.first[group0]->has_reduced_phase);
        CHECK(result.first[group0]->bus_phase_mask.size() == 2);
        CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "a"});
        CHECK_FALSE(result.first[group1]->has_reduced_phase);
        CHECK(result.first[group1]->bus_phase_mask.empty());
        CHECK(node_phase_letters(result, {2, 3}) == Letters{"abc", "abc"});
    }

    SUBCASE("A disconnected branch does not carry phases") {
        PhaseTestNetwork const network{2, {0}, {{.nodes = {0, 1}, .incidence = pp(phase_mask_a), .connected = {1, 0}}}};
        auto const result = network.build();
        CHECK(result.second->node[1].group == disconnected);
        REQUIRE(result.first.size() == 1);
        CHECK_FALSE(result.first[0]->has_reduced_phase);
        CHECK(result.first[0]->bus_phase_mask.empty());
    }

    SUBCASE("Multiple sources in one math model all feed all phases") {
        PhaseTestNetwork const network{
            3,
            {0, 2},
            {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)}, {.nodes = {1, 2}, .incidence = pp(phase_mask_a)}}};
        auto const result = network.build();
        REQUIRE(result.first.size() == 1);
        CHECK(result.first[0]->has_reduced_phase);
        CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "a", "abc"});
    }
}

TEST_CASE("Test topology single-phase transformer windings") {
    SUBCASE("Propagation") {
        SUBCASE("A phase-to-phase primary at an abc node feeds a phase-to-earth secondary") {
            PhaseTestNetwork const network{
                2, {0}, {{.nodes = {0, 1}, .incidence = winding(phase_mask_ab, phase_mask_a)}}};
            CHECK(node_phase_letters(network.build(), {0, 1}) == Letters{"abc", "a"});
        }
        SUBCASE("A phase-to-earth transformer at a SWER node feeds the LV node") {
            // 22 kV abc -ab/a- 12.7 kV SWER -a- SWER -a/a- 240 V
            PhaseTestNetwork const network{4,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = winding(phase_mask_ab, phase_mask_a)},
                                            {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {2, 3}, .incidence = winding(phase_mask_a, phase_mask_a)}}};
            CHECK(node_phase_letters(network.build(), {0, 1, 2, 3}) == Letters{"abc", "a", "a", "a"});
        }
        SUBCASE("A phase-to-phase secondary") {
            PhaseTestNetwork const network{
                2, {0}, {{.nodes = {0, 1}, .incidence = winding(phase_mask_bc, phase_mask_ab)}}};
            CHECK(node_phase_letters(network.build(), {0, 1}) == Letters{"abc", "ab"});
        }
        SUBCASE("Reverse feeding gives the primary its winding phases only") {
            // only the secondary side is fed: the primary node gets ab, not abc
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                            {.nodes = {2, 1}, .incidence = winding(phase_mask_ab, phase_mask_a)}}};
            CHECK(node_phase_letters(network.build(), {0, 1, 2}) == Letters{"abc", "a", "ab"});
        }
    }

    SUBCASE("Missing winding phases (D1)") {
        SUBCASE("A phase-to-phase primary at a single-phase node") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                      {.nodes = {1, 2}, .incidence = winding(phase_mask_ab, phase_mask_a)}}};
            network.comp_conn.branch_id = {10, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), MissingPhaseForWinding);
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("(ab)"), MissingPhaseForWinding);
        }
        SUBCASE("A phase-to-earth winding at a node of another phase") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_b)},
                                      {.nodes = {1, 2}, .incidence = winding(phase_mask_a, phase_mask_a)}}};
            network.comp_conn.branch_id = {10, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), MissingPhaseForWinding);
        }
        SUBCASE("A winding with fewer terminals than its node has phases is fine") {
            PhaseTestNetwork const network{3,
                                           {0},
                                           {{.nodes = {0, 1}, .incidence = winding(phase_mask_a, phase_mask_bc)},
                                            {.nodes = {0, 2}, .incidence = winding(phase_mask_ac, phase_mask_c)}}};
            TopologyResult result;
            REQUIRE_NOTHROW(result = network.build());
            CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "bc", "c"});
        }
    }

    SUBCASE("The phase shift of the windings accumulates along the feeder") {
        // ab -> a: theta_from - theta_to = -pi/6; a -> a with reversed polarity: pi
        constexpr double shift_ab_a = -pi / 6.0;
        PhaseTestNetwork const network{
            4,
            {0},
            {{.nodes = {0, 1}, .incidence = winding(phase_mask_ab, phase_mask_a), .phase_shift = shift_ab_a},
             {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
             {.nodes = {2, 3}, .incidence = winding(phase_mask_a, phase_mask_a), .phase_shift = pi}}};
        auto const result = network.build();
        CHECK(node_phase_shift(result, 0) == 0.0);
        CHECK(node_phase_shift(result, 1) == doctest::Approx(pi / 6.0));
        CHECK(node_phase_shift(result, 2) == doctest::Approx(pi / 6.0));
        CHECK(node_phase_shift(result, 3) == doctest::Approx(pi / 6.0 - pi));
    }
}

TEST_CASE("Test topology centre-tapped split-phase slot patterns") {
    using enum SlotPattern;
    // the split-phase pattern of a centre-tapped winding spreads along per-phase branches and links, never through a
    // winding, and only holds at nodes without phase c
    SUBCASE("A phase-to-earth primary feeds a split-phase node") {
        PhaseTestNetwork const network{2, {0}, {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a)}}};
        auto const result = network.build();
        REQUIRE(result.first.size() == 1);
        CHECK(result.first[0]->bus_slot_pattern.size() == 2);
        CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "ab"});
        CHECK(node_slot_pattern(result, 0) == three_phase);
        CHECK(node_slot_pattern(result, 1) == split_phase);
    }
    SUBCASE("An ab line off a split-phase node is split") {
        PhaseTestNetwork const network{3,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_bc)},
                                        {.nodes = {1, 2}, .incidence = pp(phase_mask_ab)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "ab", "ab"});
        CHECK(node_slot_pattern(result, 1) == split_phase);
        CHECK(node_slot_pattern(result, 2) == split_phase);
    }
    SUBCASE("A single-leg line off a split-phase node is split") {
        // the line is fed from the to side of a second line, so the pattern also spreads against the branch direction
        PhaseTestNetwork const network{4,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {3, 1}, .incidence = pp(phase_mask_b)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2, 3}) == Letters{"abc", "ab", "a", "b"});
        CHECK(node_slot_pattern(result, 2) == split_phase);
        CHECK(node_slot_pattern(result, 3) == split_phase);
    }
    SUBCASE("A winding off a split-phase node blocks the pattern") {
        PhaseTestNetwork const network{3,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = winding(phase_mask_a, phase_mask_a)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "ab", "a"});
        CHECK(node_slot_pattern(result, 1) == split_phase);
        CHECK(node_slot_pattern(result, 2) == three_phase);
    }
    SUBCASE("A centre-tapped to side at a node that is also fed three-phase stays three-phase") {
        // node 1 is abc through the line from the source; the a lines to nodes 2 and 3 make the math model reduced
        PhaseTestNetwork const network{4,
                                       {0},
                                       {{.nodes = {0, 1}},
                                        {.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a)},
                                        {.nodes = {1, 2}, .incidence = pp(phase_mask_a)},
                                        {.nodes = {0, 3}, .incidence = pp(phase_mask_a)}}};
        auto const result = network.build();
        REQUIRE(result.first[0]->has_reduced_phase);
        CHECK(node_phase_letters(result, {0, 1, 2, 3}) == Letters{"abc", "abc", "a", "a"});
        CHECK(node_slot_pattern(result, 1) == three_phase);
        // a three-phase node does not pass the pattern on
        CHECK(node_slot_pattern(result, 2) == three_phase);
        CHECK(node_slot_pattern(result, 3) == three_phase);
    }
    SUBCASE("The pattern does not leak through a three-phase node into another lateral") {
        // 0 (abc) -abc- 1 (abc) -a/ct- 2 (ab); node 0 also feeds an a lateral to node 3, which stays three-phase
        // although the per-phase arcs link it to the split-phase node 2 through nodes 1 and 0
        PhaseTestNetwork const network{4,
                                       {0},
                                       {{.nodes = {0, 1}},
                                        {.nodes = {1, 2}, .incidence = winding_ct(phase_mask_a)},
                                        {.nodes = {2, 1}, .incidence = pp(phase_mask_ab)},
                                        {.nodes = {0, 3}, .incidence = pp(phase_mask_a)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2, 3}) == Letters{"abc", "abc", "ab", "a"});
        CHECK(node_slot_pattern(result, 2) == split_phase);
        CHECK(node_slot_pattern(result, 1) == three_phase);
        CHECK(node_slot_pattern(result, 3) == three_phase);
    }
    SUBCASE("A link carries the pattern") {
        PhaseTestNetwork const network{
            3,
            {0},
            {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a)}, {.nodes = {1, 2}, .incidence = link_incidence}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1, 2}) == Letters{"abc", "ab", "ab"});
        CHECK(node_slot_pattern(result, 2) == split_phase);
    }
    SUBCASE("A disconnected to side does not seed the pattern") {
        // node 1 is fed by the ab line only
        PhaseTestNetwork const network{2,
                                       {0},
                                       {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_a), .connected = {1, 0}},
                                        {.nodes = {0, 1}, .incidence = pp(phase_mask_ab)}}};
        auto const result = network.build();
        CHECK(node_phase_letters(result, {0, 1}) == Letters{"abc", "ab"});
        CHECK(node_slot_pattern(result, 1) == three_phase);
    }
    SUBCASE("Missing winding phases (D1)") {
        SUBCASE("A phase-to-phase primary at a single-phase node") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_a)},
                                      {.nodes = {1, 2}, .incidence = winding_ct(phase_mask_ab)}}};
            network.comp_conn.branch_id = {10, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), MissingPhaseForWinding);
        }
        SUBCASE("A centre-tapped to side at a node without a and b") {
            PhaseTestNetwork network{3,
                                     {0},
                                     {{.nodes = {0, 1}, .incidence = pp(phase_mask_c)},
                                      {.nodes = {2, 1}, .incidence = winding_ct(phase_mask_a)}}};
            network.comp_conn.branch_id = {10, 42};
            CHECK_THROWS_WITH_AS(network.build(), doctest::Contains("branch with ID 42"), MissingPhaseForWinding);
        }
    }
    SUBCASE("The phase shift of a centre-tapped winding accumulates") {
        // ab -> ct: theta_from - theta_to = 0 - pi/6
        PhaseTestNetwork const network{
            3,
            {0},
            {{.nodes = {0, 1}, .incidence = winding_ct(phase_mask_ab), .phase_shift = -pi / 6.0},
             {.nodes = {1, 2}, .incidence = pp(phase_mask_ab)}}};
        auto const result = network.build();
        CHECK(node_phase_shift(result, 1) == doctest::Approx(pi / 6.0));
        CHECK(node_phase_shift(result, 2) == doctest::Approx(pi / 6.0));
        CHECK(node_slot_pattern(result, 2) == split_phase);
    }
}

TEST_CASE("Test cycle reorder") {
    SUBCASE("9 nodes") {
        // {
        //     0: [3, 5],
        //     1: [4, 5, 8],
        //     2: [4, 5, 6],
        //     3: [6, 7],
        //     4: [6, 8],
        //     6: [7, 8, 9],
        //     7: [8, 9],
        //     8: [9]
        // }

        // component topology
        ReducedComponentTopology comp_topo{};
        comp_topo.n_node = 10;
        comp_topo.branch_node_idx = {
            {0, 3}, {0, 5}, {1, 4}, {1, 5}, {1, 8}, {2, 4}, {2, 5}, {2, 6}, {3, 6},
            {3, 7}, {4, 6}, {4, 8}, {6, 7}, {6, 8}, {6, 9}, {7, 8}, {7, 9}, {8, 9},
        };
        comp_topo.source_node_idx = {0};
        // component connection
        ComponentConnections comp_conn{};
        comp_conn.branch_connected = std::vector<BranchConnected>(18, {1, 1});
        comp_conn.branch_phase_shift = std::vector<double>(18, 0.0);
        comp_conn.source_connected = {1};
        // result
        TopologicalComponentToMathCoupling comp_coup_ref{};
        comp_coup_ref.node = {{.group = 0, .pos = 0}, {.group = 0, .pos = 1}, {.group = 0, .pos = 2},
                              {.group = 0, .pos = 3}, {.group = 0, .pos = 4}, {.group = 0, .pos = 5},
                              {.group = 0, .pos = 6}, {.group = 0, .pos = 7}, {.group = 0, .pos = 8},
                              {.group = 0, .pos = 9}};
        std::vector<BranchIdx> const fill_in_ref{{3, 5}, {4, 5}, {5, 8}, {5, 6}, {5, 7}};

        Topology topo{comp_topo, comp_conn};
        auto pair = topo.build_topology();
        auto const& topo_comp_coup = *pair.second;
        auto const& math_topo = *pair.first[0];
        CHECK(topo_comp_coup.node == comp_coup_ref.node);
        CHECK(math_topo.fill_in == fill_in_ref);
    }

    SUBCASE("7 nodes") {
        // {
        //     0: [1, 5, 6],
        //     1: [2, 4, 5, 3, 6],
        //     2: [3, 6]
        //     3: [4],
        //     4: [5],
        // }

        // component topology
        ReducedComponentTopology comp_topo{};
        comp_topo.n_node = 7;
        comp_topo.branch_node_idx = {
            {0, 1}, // 0
            {1, 2}, // 1
            {2, 3}, // 2
            {3, 4}, // 3
            {4, 5}, // 4
            {0, 5}, // 5
            {1, 4}, // 6
            {6, 0}, // 7
            {6, 2}, // 8
            {5, 1}, // 9
            {3, 1}, // 10
            {6, 1}, // 11
            {2, 1}, // 12
        };
        comp_topo.source_node_idx = {0};
        // component connection
        ComponentConnections comp_conn{};
        comp_conn.branch_connected = std::vector<BranchConnected>(13, {1, 1});
        comp_conn.branch_phase_shift = std::vector<double>(13, 0.0);
        comp_conn.source_connected = {1};
        // result
        TopologicalComponentToMathCoupling comp_coup_ref{};
        comp_coup_ref.node = {{.group = 0, .pos = 0}, {.group = 0, .pos = 3}, {.group = 0, .pos = 1},
                              {.group = 0, .pos = 2}, {.group = 0, .pos = 4}, {.group = 0, .pos = 5},
                              {.group = 0, .pos = 6}};
        std::vector<BranchIdx> const fill_in_ref{{5, 6}, {2, 6}, {4, 6}};

        Topology topo{comp_topo, comp_conn};
        auto pair = topo.build_topology();
        auto const& topo_comp_coup = *pair.second;
        auto const& math_topo = *pair.first[0];
        CHECK(topo_comp_coup.node == comp_coup_ref.node);
        CHECK(math_topo.fill_in == fill_in_ref);
    }
}

} // namespace power_grid_model::topology
