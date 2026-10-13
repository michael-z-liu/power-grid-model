// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include "power_grid_model/common/counting_iterator.hpp"
#include <power_grid_model/math_solver/y_bus.hpp>

#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/grouped_index_vector.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>

#include <doctest/doctest.h>

#include <complex>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

namespace power_grid_model {

namespace {
using math_solver::YBusStructure;

// asym input
// Symmetrical parameters and admittances are converted to asymmetrical tensors,
// i.e. each parameter/admittance x is converted to:
//   x 0 0
//   0 x 0
//   0 0 x
auto convert_to_asymmetrical(MathModelParam<symmetric_t> const& param_sym,
                             ComplexTensorVector<symmetric_t> const& admittance_sym) {
    std::pair<MathModelParam<asymmetric_t>, ComplexTensorVector<asymmetric_t>> result{
        {}, ComplexTensorVector<asymmetric_t>(admittance_sym.size())};

    auto& [param_asym, admittance_asym] = result;

    // value
    param_asym.branch_param.resize(param_sym.branch_param.size());
    for (size_t i = 0; i < param_sym.branch_param.size(); i++) {
        for (size_t j = 0; j < 4; j++) {
            param_asym.branch_param[i].value[j] = ComplexTensor<asymmetric_t>{param_sym.branch_param[i].value[j]};
        }
    }
    param_asym.shunt_param.resize(param_sym.shunt_param.size());
    for (size_t i = 0; i < param_sym.shunt_param.size(); i++) {
        param_asym.shunt_param[i] = ComplexTensor<asymmetric_t>{param_sym.shunt_param[i]};
    }
    // admittance_sym
    for (size_t i = 0; i < admittance_sym.size(); i++) {
        admittance_asym[i] = ComplexTensor<asymmetric_t>{admittance_sym[i]};
    }
    return result;
}
} // namespace

TEST_CASE("Test y bus") {
    //     test Y bus struct
    //     [
    //             x, x, 0, 0
    //             x, x, x, 0
    //             0, x, x, x
    //             0, 0, x, x
    //     ]

    //      [0]   = Node
    // --0--> = Branch (from --id--> to)
    //  -X-   = Open switch / not connected

    //     Topology:

    //   --- 4 ---            /-3-\  }
    //  |         |           |   |  } loop from 2 to 2
    //  |         v           v   |  }
    // [0]       [1] --- 1 --> [2] --- 2 --> [3]
    //  ^         |           /  ^            |
    //  |         |           5  |            |
    //   --- 0 ---            |   \--- 6 ----/
    //                        X
    MathModelTopology topo{};
    MathModelParam<symmetric_t> param_sym;
    topo.phase_shift.resize(4, 0.0);
    topo.branch_bus_idx = {
        {1, 0}, // branch 0 from node 1 to 0
        {1, 2}, // branch 1 from node 1 to 2
        {2, 3}, // branch 2 from node 2 to 3
        {2, 2}, // branch 3 from node 2 to 2 (loop into itself)
        {3, 2}, // branch 6 from node 3 to 2
        {0, 1}, // branch 4 from node 0 to 1
        {2, -1} // branch 5 from node 2 to "not connected"
    };
    param_sym.branch_param = {
        // ff, ft, tf, tt
        {1.0i, 2.0i, 3.0i, 4.0i}, {5.0, 6.0, 7.0, 8.0},     {9.0i, 10.0i, 11.0i, 12.0i}, {21.0i, 22.0i, 22.0i, 21.0i},
        {13.0, 14.0, 15.0, 16.0}, {17.0, 18.0, 19.0, 20.0}, {1000i, 0.0, 0.0, 0.0},
    };
    topo.shunts_per_bus = {from_sparse, {0, 1, 1, 1, 2}}; // 4 buses, 2 shunts -> shunt connected to bus 0 and bus 3
    param_sym.shunt_param = {100.0i, 200.0i};

    // output
    IdxVector const row_indptr = {0, 2, 5, 8, 10};

    // Use col_indices to find the location in Y bus
    //  e.g. col_indices = {0, 1, 0} results in Y bus:
    // [
    //	x, x
    //   x, 0
    // ]
    IdxVector const col_indices = {// Column col_indices for each non-zero element in Y bus.
                                   0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    Idx nnz = 10; // Number of non-zero elements in Y bus
    IdxVector const bus_entry = {0, 3, 6, 9};
    IdxVector const lu_transpose_entry = {// Flip the id's of non-diagonal elements
                                          0, 2, 1, 3, 5, 4, 6, 8, 7, 9};
    IdxVector const y_bus_entry_indptr = {0,  3,      // 0, 1, 2 belong to element [0,0] in Ybus /  3,4 to element [0,1]
                                          5,  7,  10, // 5,6 to [1,0] / 7, 8, 9 to [1,1] / 10 to [1,2]
                                          11, 12,     // 11 to [2,1] / 12, 13, 14, 15, 16, 17, 18, 19 to [2,2]
                                          20, 22,     // 20, 21 to [2,3] / 22, 23 to [3, 2]
                                          24, 27};    // 24, 25, 26  to [3,3]
    IdxVector map_lu_y_bus = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    ComplexTensorVector<symmetric_t> admittance_sym = {
        4.0i + 17.0 + 100.0i,                 // 0, 0 -> {1, 0}tt + {0, 1}ff + shunt(0)
        18.0 + 3.0i,                          // 0, 1 -> {0, 1}ft + {1, 0}tf
        19.0 + 2.0i,                          // 1, 0 -> {0, 1}tf + {1, 0}ft
        20.0 + 1.0i + 5.0,                    // 1, 1 -> {0, 1}tt + {1, 0}ff + {1,2}ff
        6.0,                                  // 1, 2 -> {1,2}ft
        7.0,                                  // 2, 1 -> {1,2}tf
        8.0                                   // 2, 2 -> {1,2}tt
            + (21.0i + 22.0i + 22.0i + 21.0i) // + ({3,3}ff + {3,3}ft + {3,3}{tf} + {3,3}tt)
            + 9.0i + 16.0 + 1000.0i,          //    + {2,3}ff + {3, 2}tt + {2,-1}ff
        10.0i + 15.0,                         // 2, 3 -> {2,3}ft + {3,2}tf
        11.0i + 14.0,                         // 3, 2 -> {2,3}tf + {3,2}ft
        12.0i + 13.0 + 200.0i                 // 3, 3 -> {2,3}tt + {3,2}ff + shunt(1)
    };

    auto const [param_asym, admittance_asym] = convert_to_asymmetrical(param_sym, admittance_sym);

    SUBCASE("Test y bus construction (symmetrical)") {
        YBus<symmetric_t> const ybus{topo, param_sym};
        CHECK(ybus.size() == 4);
        CHECK(ybus.nnz() == nnz);
        CHECK(row_indptr == ybus.row_indptr());
        CHECK(col_indices == ybus.col_indices());
        CHECK(bus_entry == ybus.bus_entry());
        CHECK(lu_transpose_entry == ybus.lu_transpose_entry());
        CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr());
        CHECK(ybus.admittance().size() == admittance_sym.size());
        for (size_t i = 0; i < admittance_sym.size(); i++) {
            CAPTURE(i);
            CHECK(cabs(ybus.admittance()[i] - admittance_sym[i]) < numerical_tolerance);
        }

        // check lu
        CHECK(ybus.map_lu_y_bus() == map_lu_y_bus);

        SUBCASE("Test y bus structure getter") {
            YBusStructure ybus_struct_ref{topo};
            auto const& ybus_struct = ybus.y_bus_structure();
            CHECK(ybus_struct.bus_entry == ybus_struct_ref.bus_entry);
            CHECK(ybus_struct.col_indices == ybus_struct_ref.col_indices);
            CHECK(ybus_struct.col_indices_lu == ybus_struct_ref.col_indices_lu);
            CHECK(ybus_struct.diag_lu == ybus_struct_ref.diag_lu);
            CHECK(ybus_struct.lu_transpose_entry == ybus_struct_ref.lu_transpose_entry);
            CHECK(ybus_struct.map_lu_y_bus == ybus_struct_ref.map_lu_y_bus);
            CHECK(ybus_struct.row_indptr == ybus_struct_ref.row_indptr);
            CHECK(ybus_struct.row_indptr_lu == ybus_struct_ref.row_indptr_lu);
            CHECK(ybus_struct.y_bus_element.size() == ybus_struct_ref.y_bus_element.size());
            CHECK(ybus_struct.y_bus_entry_indptr == ybus_struct_ref.y_bus_entry_indptr);
        }

        SUBCASE("Test shared y bus structure getter") {
            YBusStructure ybus_struct_ref{topo};
            auto const& ybus_struct = ybus.shared_y_bus_structure();
            CHECK(ybus_struct->bus_entry == ybus_struct_ref.bus_entry);
            CHECK(ybus_struct->col_indices == ybus_struct_ref.col_indices);
            CHECK(ybus_struct->col_indices_lu == ybus_struct_ref.col_indices_lu);
            CHECK(ybus_struct->diag_lu == ybus_struct_ref.diag_lu);
            CHECK(ybus_struct->lu_transpose_entry == ybus_struct_ref.lu_transpose_entry);
            CHECK(ybus_struct->map_lu_y_bus == ybus_struct_ref.map_lu_y_bus);
            CHECK(ybus_struct->row_indptr == ybus_struct_ref.row_indptr);
            CHECK(ybus_struct->row_indptr_lu == ybus_struct_ref.row_indptr_lu);
            CHECK(ybus_struct->y_bus_element.size() == ybus_struct_ref.y_bus_element.size());
            CHECK(ybus_struct->y_bus_entry_indptr == ybus_struct_ref.y_bus_entry_indptr);
        }
    }

    SUBCASE("Test y bus construction (asymmetrical)") {
        YBus<symmetric_t> const ybus_sym{topo, param_sym};
        // construct from existing structure
        YBus<asymmetric_t> const ybus{topo, param_asym, ybus_sym.shared_y_bus_structure()};
        CHECK(ybus.size() == 4);
        CHECK(ybus.nnz() == nnz);
        CHECK(row_indptr == ybus.row_indptr());
        CHECK(col_indices == ybus.col_indices());
        CHECK(bus_entry == ybus.bus_entry());
        CHECK(lu_transpose_entry == ybus.lu_transpose_entry());
        CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr());
        CHECK(ybus.admittance().size() == admittance_asym.size());
        for (size_t i = 0; i < admittance_asym.size(); i++) {
            CHECK((cabs(ybus.admittance()[i] - admittance_asym[i]) < numerical_tolerance).all());
        }
    }

    SUBCASE("Test branch flow calculation") {
        YBus<symmetric_t> const ybus{topo, param_sym};
        ComplexVector const u{1.0, 2.0, 3.0, 4.0};
        auto branch_flow = ybus.calculate_branch_flow<BranchSolverOutput<symmetric_t>>(u);

        // branch 2, bus 2->3
        // if = 3 * 9i + 4 * 10i = 67i
        // it = 3 * 11i + 4 * 12i = 81i
        // sf = 3 * conj(67i) = -201i
        // st = 4 * conj(81i) = -324i
        CHECK(cabs(branch_flow[2].i_f - 67.0i) < numerical_tolerance);
        CHECK(cabs(branch_flow[2].i_t - 81.0i) < numerical_tolerance);
        CHECK(cabs(branch_flow[2].s_f - (-201.0i)) < numerical_tolerance);
        CHECK(cabs(branch_flow[2].s_t - (-324.0i)) < numerical_tolerance);
    }

    SUBCASE("Test shunt flow calculation") {
        YBus<symmetric_t> const ybus{topo, param_sym};
        ComplexVector const u{1.0, 2.0, 3.0, 4.0};
        auto shunt_flow = ybus.template calculate_shunt_flow<ApplianceSolverOutput<symmetric_t>>(u);

        // shunt 1
        // i = -4 * 200i
        // s = 4 * conj(-800i) = 3200i
        CHECK(cabs(shunt_flow[1].i - (-800.0i)) < numerical_tolerance);
        CHECK(cabs(shunt_flow[1].s - 3200.0i) < numerical_tolerance);
    }
}

TEST_CASE("Test one bus system") {
    MathModelTopology topo{};
    MathModelParam<symmetric_t> param_sym;
    topo.phase_shift = {0.0};

    // output
    IdxVector const indptr = {0, 1};
    IdxVector const col_indices = {0};
    Idx const nnz = 1;
    IdxVector const bus_entry = {0};
    IdxVector const lu_transpose_entry = {0};

    SUBCASE("One shunt") {
        topo.shunts_per_bus = {from_sparse, {0, 0}};
        IdxVector const y_bus_entry_indptr = {0, 0};

        YBus<symmetric_t> const ybus{topo, param_sym};

        CHECK(ybus.size() == 1);
        CHECK(ybus.nnz() == nnz);
        CHECK(indptr == ybus.row_indptr());
        CHECK(col_indices == ybus.col_indices());
        CHECK(bus_entry == ybus.bus_entry());
        CHECK(lu_transpose_entry == ybus.lu_transpose_entry());
        CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr());
    }
    SUBCASE("Branch into itself") {
        topo.branch_bus_idx = {{0, 0}};
        param_sym.branch_param = {// ff, ft,   tf,   tt
                                  {1.0i, 2.0i, 2.0i, 1.0i}};

        IdxVector const y_bus_entry_indptr = {0, 4};

        YBus<symmetric_t> const ybus{topo, param_sym};

        CHECK(ybus.size() == 1);
        CHECK(ybus.nnz() == nnz);
        CHECK(indptr == ybus.row_indptr());
        CHECK(col_indices == ybus.col_indices());
        CHECK(bus_entry == ybus.bus_entry());
        CHECK(lu_transpose_entry == ybus.lu_transpose_entry());
        CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr());
    }
}

TEST_CASE("Test cyclic three bus system with branch into itself") {

    //     test Y bus struct
    //     [
    //             x, x, x
    //             x, x, x
    //             x, x, x
    //     ]

    //      [0]   = Node
    // --0--> = Branch (from --id--> to)
    //  -X-   = Open switch / not connected

    //     Topology:

    //                  /- 3 -\  }
    //                  |     |  } loop from 1 to 1
    //      /--- 0 --> [1] <--/  }
    //      |           |
    //     [0]          1
    //      ^           v
    //      \--- 2 --- [2] <--\  }
    //                  |     |  } loop from 2 to 2
    //                  \- 4 -/  }
    MathModelTopology topo{};
    MathModelParam<symmetric_t> param_sym;
    topo.phase_shift.resize(3, 0.0);
    topo.branch_bus_idx = {
        {0, 1}, // branch 0 from node 0 to 1
        {1, 2}, // branch 1 from node 1 to 2
        {2, 0}, // branch 2 from node 2 to 0
        {1, 1}, // branch 3 from node 1 to 1 (loop into itself)
        {2, 2}, // branch 4 from node 2 to 2 (loop into itself)
    };
    param_sym.branch_param = {
        // ff, ft, tf, tt
        {1.0i, 2.0i, 3.0i, 4.0i},     // 0 -> 1
        {5.0, 6.0, 7.0, 8.0},         // 1 -> 2
        {9.0i, 10.0i, 11.0i, 12.0i},  // 2 -> 0
        {21.0i, 22.0i, 22.0i, 21.0i}, // 1 -> 1
        {13.0, 14.0, 15.0, 16.0},     // 2 -> 2
    };
    topo.shunts_per_bus = {from_dense, {}, 3};
    param_sym.shunt_param = {};

    // output
    IdxVector const row_indptr = {0, 3, 6, 9};

    // Use col_indices to find the location in Y bus
    //  e.g. col_indices = {0, 1, 0} results in Y bus:
    // [
    //	x, x
    //   x, 0
    // ]
    IdxVector const col_indices = {// Column col_indices for each non-zero element in Y bus.
                                   0, 1, 2, 0, 1, 2, 0, 1, 2};
    Idx const nnz = 9; // Number of non-zero elements in Y bus
    IdxVector const bus_entry = {0, 4, 8};
    IdxVector const lu_transpose_entry = {// Flip the id's of non-diagonal elements
                                          0, 3, 6, 1, 4, 7, 2, 5, 8};
    IdxVector const y_bus_entry_indptr = {
        0,  2, // 0, 1 belong to element [0,0] in Ybus
        3,     // 2 to element [0,1]
        4,     // 3 to [0, 2]
        5,     // 4 to [1, 0]
        11,    // 5, 6, (7, 8, 9, 10) to [1,1]
        12,    // 11 to [1,2]
        13,    // 12 to [2,0]
        14,    // 13 to [2,1]
        20,    // 14, 15, (16, 17, 18, 19) to [2,2]
    };
    IdxVector const map_lu_y_bus = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    ComplexTensorVector<symmetric_t> const admittance_sym = {
        1.0i + 12.0i, // 0, 0 -> {0, 1}ff + {2, 0}tt
        2.0i,         // 0, 1 -> {0, 1}ft
        11.0i,        // 0, 2 -> {2, 0}tf
        3.0i,         // 1, 0 -> {0, 1}tf
        4.0i + 5.0 +
            (21.0i + 22.0i + 22.0i + 21.0i), // 1, 1 -> {0, 1}tt + {1, 2}ff + ({1,1}ff + {1,1}ft + {1,1}tf + {1,1}tt)
        6.0,                                 // 1, 2 -> {1,2}ft
        10.0i,                               // 2, 0 -> {2, 0}ft
        7.0,                                 // 2, 1 -> {1,2}tf
        8.0 + 9.0i +
            (13.0 + 14.0 + 15.0 + 16.0), // 2, 2 -> {1, 2}tt + {2, 0}ff + ({2,2}ff + {2,2}ft + {2,2}tf + {2,2}tt)
    };

    YBus<symmetric_t> const ybus{topo, param_sym};

    CHECK(ybus.size() == 3);
    CHECK(ybus.nnz() == nnz);
    CHECK(row_indptr == ybus.row_indptr());
    CHECK(col_indices == ybus.col_indices());
    CHECK(bus_entry == ybus.bus_entry());
    CHECK(lu_transpose_entry == ybus.lu_transpose_entry());
    CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr());
    CHECK(ybus.admittance().size() == admittance_sym.size());
    for (size_t i = 0; i < admittance_sym.size(); i++) {
        CAPTURE(i);
        CHECK(cabs(ybus.admittance()[i] - admittance_sym[i]) < numerical_tolerance);
    }
}

TEST_CASE("Test fill-in y bus") {
    // [1] --0--> [0] --1--> [2]
    // extra fill-in: (1, 2) by removing node 0
    //
    // [
    //   0, 1, 2
    //   3, 4, f
    //   5, f, 6
    // ]
    MathModelTopology topo{};
    topo.phase_shift.resize(3, 0.0);
    topo.branch_bus_idx = {
        {1, 0}, // branch 0 from node 1 to 0
        {0, 2}, // branch 1 from node 0 to 2
    };
    topo.shunts_per_bus = {from_sparse, {0, 0, 0, 0}};
    topo.fill_in = {{1, 2}};

    IdxVector const row_indptr = {0, 3, 5, 7};
    IdxVector const col_indices = {0, 1, 2, 0, 1, 0, 2};
    IdxVector const bus_entry = {0, 4, 6};
    IdxVector const lu_transpose_entry = {0, 3, 6, 1, 4, 7, 2, 5, 8};
    IdxVector const y_bus_entry_indptr = {0, 2,              // 0, 1 belong to element [0,0] in Ybus
                                          3, 4, 5, 6, 7, 8}; // everything else has only one entry
    // lu matrix
    IdxVector const row_indptr_lu = {0, 3, 6, 9};
    IdxVector const col_indices_lu = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    IdxVector const map_lu_y_bus = {0, 1, 2, 3, 4, -1, 5, -1, 6};
    IdxVector const diag_lu = {0, 4, 8};

    YBusStructure const ybus{topo};

    CHECK(row_indptr == ybus.row_indptr);
    CHECK(col_indices == ybus.col_indices);
    CHECK(bus_entry == ybus.bus_entry);
    CHECK(lu_transpose_entry == ybus.lu_transpose_entry);
    CHECK(y_bus_entry_indptr == ybus.y_bus_entry_indptr);
    // check lu
    CHECK(ybus.row_indptr_lu == row_indptr_lu);
    CHECK(ybus.col_indices_lu == col_indices_lu);
    CHECK(ybus.diag_lu == diag_lu);
    CHECK(ybus.map_lu_y_bus == map_lu_y_bus);
}

TEST_CASE("Incremental update y-bus") {
    // test Y bus struct
    // [
    //         x, x, 0, 0
    //         x, x, x, 0
    //         0, x, x, x
    //         0, 0, x, x
    // ]
    //
    //  [0]   = Node
    // --0--> = Branch (from --id--> to)
    //  -X-   = Open switch / not connected
    //
    //     Topology:
    //
    //   --- 4 ---               ----- 3 -----
    //  |         |             |             |
    //  |         v             v             |
    // [0]       [1] --- 1 --> [2] --- 2 --> [3]
    //  ^         |             |
    //  |         |             5
    //   --- 0 ---              |
    //                          X
    MathModelTopology topo{};
    MathModelParam<symmetric_t> param_sym;
    topo.phase_shift.resize(4, 0.0);
    topo.branch_bus_idx = {
        {1, 0}, // branch 0 from node 1 to 0
        {1, 2}, // branch 1 from node 1 to 2
        {2, 3}, // branch 2 from node 2 to 3
        {3, 2}, // branch 3 from node 3 to 2
        {0, 1}, // branch 4 from node 0 to 1
        {2, -1} // branch 5 from node 2 to "not connected"
    };
    param_sym.branch_param = {
        //  ff,    ft,    tf,   tt
        {1.0i, 2.0i, 3.0i, 4.0i},    // (1, 0)
        {5.0, 6.0, 7.0, 8.0},        // (1, 2)
        {9.0i, 10.0i, 11.0i, 12.0i}, // (2, 3)
        {13.0, 14.0, 15.0, 16.0},    // (3, 2)
        {17.0, 18.0, 19.0, 20.0},    // (0, 1)
        {1000i, 0.0, 0.0, 0.0}       // (2, -1)
    };
    topo.shunts_per_bus = {from_sparse, {0, 1, 1, 1, 2}}; // 4 buses, 2 shunts -> shunt connected to bus 0 and bus 3
    param_sym.shunt_param = {100.0i, 200.0i};

    const ComplexTensorVector<symmetric_t> admittance_sym = {
        17.0 + 104.0i,  // 0, 0 -> {1, 0}tt + {0, 1}ff + shunt(0) = 4.0i + 17.0 + 100.0i
        18.0 + 3.0i,    // 0, 1 -> {0, 1}ft + {1, 0}tf = 18.0 + 3.0i
        19.0 + 2.0i,    // 1, 0 -> {0, 1}tf + {1, 0}ft = 19.0 + 2.0i
        25.0 + 1.0i,    // 1, 1 -> {0, 1}tt + {1, 0}ff + {1,2}ff = 20.0 + 1.0i + 5.0
        6.0,            // 1, 2 -> {1,2}ft = 6.0
        7.0,            // 2, 1 -> {1,2}tf = 7.0
        24.0 + 1009.0i, // 2, 2 -> {1,2}tt + {2,3}ff + {3, 2}tt + {2,-1}ff = 8.0 + 9.0i + 16.0 + 1000.0i = 24.0 + 1009i
        15.0 + 10.0i,   // 2, 3 -> {2,3}ft + {3,2}tf = 10.0i + 15.0
        14.0 + 11.0i,   // 3, 2 -> {2,3}tf + {3,2}ft = 11.0i + 14.0
        13.0 + 212.0i   // 3, 3 -> {2,3}tt + {3,2}ff + shunt(1) = 12.0i + 13.0 + 200.0i
    };

    const ComplexTensorVector<symmetric_t> admittance_sym_state_1 = {
        34.0 + 208.0i, 36.0 + 6.0i,    38.0 + 4.0i,  50.0 + 2.0i,  12.0,
        14.0,          48.0 + 2018.0i, 30.0 + 20.0i, 14.0 + 22.0i, 26.0 + 424.0i};
    topo.branch_bus_idx = {
        {1, 0}, // branch 0 from node 1 to 0
        {1, 2}, // branch 1 from node 1 to 2
        {2, 3}, // branch 2 from node 2 to 3
        {3, 2}, // branch 3 from node 3 to 2
        {0, 1}, // branch 4 from node 0 to 1
        {2, -1} // branch 5 from node 2 to "not connected"
    };

    MathModelParam<symmetric_t> param_sym_update;
    // Swap based params
    param_sym_update.branch_param = {
        //   ff,    ft,   tf,   tt
        {2.0i, 2.0i, 3.0i, 4.0i},    // (1, 0)
        {5.0, 7.0, 7.0, 8.0},        // (1, 2)
        {9.0i, 10.0i, 11.0i, 14.0i}, // (2, 3)
        {13.0, 14.0, 17.0, 16.0},    // (3, 2)
        {17.0, 18.0, 19.0, 20.0},    // (0, 1)
        {1001.0i, 0.0, 0.0, 0.0}     // (2,-1)
    };
    param_sym_update.shunt_param = {101.0i, 200.0i};

    const ComplexTensorVector<symmetric_t> admittance_sym_2 = {
        // 17.0 + 104.0i,                                                                        [v]
        17.0 + 105.0i, // 0, 0 -> += {1, 0}tt + {0, 1}ff + shunt(0) = 0.0 + 0.0 + 1.0i
        // 18.0 + 3.0i,                                                                          [v]
        18.0 + 3.0i, // 0, 1 -> += {0, 1}ft + {1, 0}tf = 0.0 + 0.0
        // 19.0 + 2.0i,                                                                          [v]
        19.0 + 2.0i, // 1, 0 -> += {0, 1}tf + {1, 0}ft = 0.0 + 0.0
        // 25.0 + 1.0i,                                                                          [v]
        25.0 + 2.0i, // 1, 1 -> += {0, 1}tt + {1, 0}ff + {1,2}ff = 0.0 + 1.0i + 0.0
        // 6.0,                                                                                  [v]
        7.0, // 1, 2 -> += {1,2}ft = 1.0
        // 7.0,                                                                                  [v]
        7.0, // 2, 1 -> += {1,2}tf = 0.0
        // 24.0 + 1009.0i,                                                                       [v]
        24.0 + 1010.0i, // 2, 2 -> += {1,2}tt + {2,3}ff + {3, 2}tt + {2,-1}ff = 0.0 + 0.0 + 0.0 + 1.0i
        // 15.0 + 10.0i,                                                                         [v]
        17.0 + 10.0i, // 2, 3 -> += {2,3}ft + {3,2}tf = 0.0 + 2.0
        // 14.0 + 11.0i,                                                                         [v]
        14.0 + 11.0i, // 3, 2 -> += {2,3}tf + {3,2}ft = 0.0 + 0.0
        // 13.0 + 212.0i                                                                         [v]
        13.0 + 214.0i // 3, 3 -> += {2,3}tt + {3,2}ff + shunt(1) = 2.0i + 0.0 + 0.0
    };

    auto verify_admittance = [](ComplexTensorVector<symmetric_t> const& admittance,
                                ComplexTensorVector<symmetric_t> const& admittance_ref) {
        CHECK(admittance.size() == admittance_ref.size());
        for (size_t i = 0; i < admittance.size(); i++) {
            CHECK(cabs(admittance[i] - admittance_ref[i]) < numerical_tolerance);
        }
    };
    SUBCASE("Test whole scale update") {
        YBus<symmetric_t> ybus{topo, param_sym};
        verify_admittance(ybus.admittance(), admittance_sym);

        ybus.update_admittance(param_sym);
        verify_admittance(ybus.admittance(), admittance_sym);
    }

    SUBCASE("Test progressive update") {
        YBus<symmetric_t> ybus{topo, param_sym};
        verify_admittance(ybus.admittance(), admittance_sym);

        auto branch_param_to_change_views =
            IdxRange{std::ssize(param_sym_update.branch_param)} | std::views::filter([&param_sym_update](Idx i) {
                return param_sym_update.branch_param[i].yff() != ComplexTensor<symmetric_t>{0.0} ||
                       param_sym_update.branch_param[i].yft() != ComplexTensor<symmetric_t>{0.0} ||
                       param_sym_update.branch_param[i].ytf() != ComplexTensor<symmetric_t>{0.0} ||
                       param_sym_update.branch_param[i].ytt() != ComplexTensor<symmetric_t>{0.0};
            }) |
            std::ranges::to<IdxVector>();
        auto shunt_param_to_change_views =
            IdxRange{std::ssize(param_sym_update.shunt_param)} | std::views::filter([&param_sym_update](Idx i) {
                return param_sym_update.shunt_param[i] != ComplexTensor<symmetric_t>{0.0};
            }) |
            std::ranges::to<IdxVector>();

        MathModelParamIncrement<symmetric_t> const math_model_param_incrmt{
            .branch_param = branch_param_to_change_views | std::views::transform([&param_sym_update](Idx i) {
                                return param_sym_update.branch_param[i];
                            }) |
                            std::ranges::to<std::vector>(),
            .shunt_param = shunt_param_to_change_views | std::views::transform([&param_sym_update](Idx i) {
                               return param_sym_update.shunt_param[i];
                           }) |
                           std::ranges::to<std::vector>(),
            .source_param = {},
            .branch_param_to_change = std::move(branch_param_to_change_views),
            .shunt_param_to_change = std::move(shunt_param_to_change_views),
            .source_param_to_change = {},
        };

        ybus.update_admittance_increment(math_model_param_incrmt);
        verify_admittance(ybus.admittance(), admittance_sym_2);
    }

    SUBCASE("Test source param incremental update") {
        YBus<symmetric_t> ybus{topo, param_sym};
        verify_admittance(ybus.admittance(), admittance_sym);

        auto source_param_to_change_views =
            IdxRange{std::ssize(param_sym_update.source_param)} | std::views::filter([&param_sym_update](Idx i) {
                return cabs(param_sym_update.source_param[i].y0) > numerical_tolerance ||
                       cabs(param_sym_update.source_param[i].y1) > numerical_tolerance;
            }) |
            std::ranges::to<IdxVector>();

        MathModelParamIncrement<symmetric_t> const math_model_param_incrmt{
            .branch_param = {},
            .shunt_param = {},
            .source_param = source_param_to_change_views | std::views::transform([&param_sym_update](Idx i) {
                                return param_sym_update.source_param[i];
                            }) |
                            std::ranges::to<std::vector>(),
            .branch_param_to_change = {},
            .shunt_param_to_change = {},
            .source_param_to_change = source_param_to_change_views,
        };

        ybus.update_admittance_increment(math_model_param_incrmt);

        verify_admittance(ybus.admittance(), admittance_sym);
        CHECK(std::ssize(ybus.math_model_param().source_param) == std::ssize(param_sym.source_param));
        for (Idx const i : source_param_to_change_views) {
            CHECK(cabs(ybus.math_model_param().source_param[i].y0 - param_sym_update.source_param[i].y0) <
                  numerical_tolerance);
            CHECK(cabs(ybus.math_model_param().source_param[i].y1 - param_sym_update.source_param[i].y1) <
                  numerical_tolerance);
        }
    }
}

TEST_CASE("Test y bus - absent phases are masked") {
    // bus0 (abc) --branch0-- bus1 (a) -- shunt0
    // expected values are derived by hand from the masking rule Y_ij <- M_i Y_ij M_j and Y_ii(p, p) = 1 for each
    // absent phase p, with exact equality where no arithmetic other than the summation of the contributions is involved
    using Tensor = ComplexTensor<asymmetric_t>;

    MathModelTopology topo{};
    topo.slack_bus = 0;
    topo.phase_shift = {0.0, 0.0};
    topo.branch_bus_idx = {{0, 1}};
    topo.sources_per_bus = {from_sparse, {0, 1, 1}};
    topo.shunts_per_bus = {from_sparse, {0, 0, 1}};
    topo.load_gens_per_bus = {from_sparse, {0, 0, 0}};
    topo.has_reduced_phase = true;
    topo.bus_phase_mask = {phase_mask_abc, phase_mask_a};

    constexpr DoubleComplex y{2.0, -4.0};
    constexpr DoubleComplex ys{0.1, 0.3};
    constexpr DoubleComplex ys0{0.4, -0.2};
    // sequence-based shunt: diagonal (2 ys + ys0) / 3, off-diagonal (ys0 - ys) / 3
    Tensor const y_shunt{(2.0 * ys + ys0) / 3.0, (ys0 - ys) / 3.0};
    Tensor y_a{0.0};
    y_a(0, 0) = y;

    MathModelParam<asymmetric_t> param;
    param.branch_param = {{y_a, -y_a, -y_a, y_a}};
    param.shunt_param = {y_shunt};
    param.source_param = {SourceCalcParam{.y1 = 1.0, .y0 = 1.0}};

    auto const entry_of = [](YBus<asymmetric_t> const& ybus, Idx row, Idx col) {
        for (Idx k = ybus.row_indptr()[row]; k != ybus.row_indptr()[row + 1]; ++k) {
            if (ybus.col_indices()[k] == col) {
                return k;
            }
        }
        FAIL("entry not found");
        return Idx{-1};
    };
    auto const check_exact = [](Tensor const& actual, Tensor const& expected) {
        for (Idx i = 0; i != 3; ++i) {
            for (Idx j = 0; j != 3; ++j) {
                CHECK(actual(i, j) == expected(i, j));
            }
        }
    };
    // bus 1 diagonal: [[y + (2 ys + ys0) / 3, 0, 0], [0, 1, 0], [0, 0, 1]]
    auto const bus1_diagonal = [&y_a](Tensor const& shunt) {
        Tensor result{0.0};
        result(0, 0) = y_a(0, 0) + shunt(0, 0);
        result(1, 1) = 1.0;
        result(2, 2) = 1.0;
        return result;
    };

    SUBCASE("Single phase branch and sequence shunt") {
        YBus<asymmetric_t> const ybus{topo, param};
        auto const& admittance = ybus.admittance();

        check_exact(admittance[entry_of(ybus, 0, 0)], y_a);
        check_exact(admittance[entry_of(ybus, 0, 1)], -y_a);
        check_exact(admittance[entry_of(ybus, 1, 0)], -y_a);
        check_exact(admittance[entry_of(ybus, 1, 1)], bus1_diagonal(y_shunt));
    }

    SUBCASE("Three phase parameters are masked per side") {
        // not a valid network (the conductors on b and c would dangle), but it exercises the row and column masks
        Tensor const y_full{(2.0 * y + 0.5) / 3.0, (0.5 - y) / 3.0};
        param.branch_param = {{y_full, -y_full, -y_full, y_full}};
        YBus<asymmetric_t> const ybus{topo, param};
        auto const& admittance = ybus.admittance();

        Tensor expected_01{0.0}; // row mask abc, column mask a
        expected_01.col(0) = -y_full.col(0);
        Tensor expected_10{0.0}; // row mask a, column mask abc
        expected_10.row(0) = -y_full.row(0);
        Tensor expected_11{0.0};
        expected_11(0, 0) = y_full(0, 0) + y_shunt(0, 0);
        expected_11(1, 1) = 1.0;
        expected_11(2, 2) = 1.0;

        check_exact(admittance[entry_of(ybus, 0, 0)], y_full);
        check_exact(admittance[entry_of(ybus, 0, 1)], expected_01);
        check_exact(admittance[entry_of(ybus, 1, 0)], expected_10);
        check_exact(admittance[entry_of(ybus, 1, 1)], expected_11);

        // flows into bus 1 are exactly zero on b and c
        ComplexValueVector<asymmetric_t> const u{ComplexValue<asymmetric_t>{1.02 * std::polar(1.0, -0.01)},
                                                 {0.95 * std::polar(1.0, -0.07), 0.0, 0.0}};
        auto const branch_flow = ybus.calculate_branch_flow<BranchSolverOutput<asymmetric_t>>(u);
        REQUIRE(branch_flow.size() == 1);
        CHECK(branch_flow[0].i_t(1) == 0.0);
        CHECK(branch_flow[0].i_t(2) == 0.0);
        CHECK(branch_flow[0].s_t(1) == 0.0);
        CHECK(branch_flow[0].s_t(2) == 0.0);
        CHECK(cabs(branch_flow[0].i_t(0) - dot(y_full, u[1])(0) + dot(y_full, u[0])(0)) < numerical_tolerance);
        // the from side is at a three phase bus and is not masked
        CHECK(branch_flow[0].i_f(1) != 0.0);
        auto const sc_branch_flow = ybus.calculate_branch_flow<BranchShortCircuitSolverOutput<asymmetric_t>>(u);
        CHECK(sc_branch_flow[0].i_t(1) == 0.0);
        CHECK(sc_branch_flow[0].i_t(2) == 0.0);
    }

    SUBCASE("Shunt flow on absent phases is exactly zero") {
        YBus<asymmetric_t> const ybus{topo, param};
        DoubleComplex const u1a = 0.95 * std::polar(1.0, -0.07);
        ComplexValueVector<asymmetric_t> const u{ComplexValue<asymmetric_t>{1.0}, {u1a, 0.0, 0.0}};

        auto const shunt_flow = ybus.calculate_shunt_flow<ApplianceSolverOutput<asymmetric_t>>(u);
        REQUIRE(shunt_flow.size() == 1);
        // without masking, the mutual term (ys0 - ys) / 3 would give a current on b and c
        CHECK(y_shunt(1, 0) != 0.0);
        CHECK(shunt_flow[0].i(1) == 0.0);
        CHECK(shunt_flow[0].i(2) == 0.0);
        CHECK(shunt_flow[0].s(1) == 0.0);
        CHECK(shunt_flow[0].s(2) == 0.0);
        CHECK(cabs(shunt_flow[0].i(0) + y_shunt(0, 0) * u1a) < numerical_tolerance);

        auto const sc_shunt_flow = ybus.calculate_shunt_flow<ApplianceShortCircuitSolverOutput<asymmetric_t>>(u);
        CHECK(sc_shunt_flow[0].i(1) == 0.0);
        CHECK(sc_shunt_flow[0].i(2) == 0.0);

        // the injection is exactly zero on the absent phases without any masking
        auto const injection = ybus.calculate_injection(u);
        CHECK(injection[1](1) == 0.0);
        CHECK(injection[1](2) == 0.0);
    }

    SUBCASE("Non-finite values on absent phases are discarded") {
        Tensor y_shunt_nan = y_shunt;
        y_shunt_nan(1, 2) = nan;
        y_shunt_nan(2, 0) = nan;
        param.shunt_param = {y_shunt_nan};
        YBus<asymmetric_t> const ybus{topo, param};
        check_exact(ybus.admittance()[entry_of(ybus, 1, 1)], bus1_diagonal(y_shunt));
    }

    SUBCASE("Incremental update keeps the identity") {
        YBus<asymmetric_t> ybus{topo, param};
        constexpr DoubleComplex ys_new{0.7, 0.9};
        Tensor const y_shunt_new{ys_new, 0.5 * ys_new};
        MathModelParamIncrement<asymmetric_t> const increment{
            .branch_param = {},
            .shunt_param = {y_shunt_new},
            .source_param = {},
            .branch_param_to_change = {},
            .shunt_param_to_change = {0},
            .source_param_to_change = {},
        };
        ybus.update_admittance_increment(increment);
        check_exact(ybus.admittance()[entry_of(ybus, 1, 1)], bus1_diagonal(y_shunt_new));
        check_exact(ybus.admittance()[entry_of(ybus, 0, 0)], y_a);
    }

    SUBCASE("All buses with all phases give the unmasked admittance") {
        param.branch_param[0] = {y_shunt, -y_shunt, -y_shunt, y_shunt};
        topo.bus_phase_mask = {phase_mask_abc, phase_mask_abc};
        MathModelTopology topo_ref = topo;
        topo_ref.has_reduced_phase = false;
        topo_ref.bus_phase_mask.clear();
        YBus<asymmetric_t> const ybus_flagged{topo, param};
        YBus<asymmetric_t> const ybus_ref{topo_ref, param};
        REQUIRE(ybus_flagged.nnz() == ybus_ref.nnz());
        for (Idx k = 0; k != ybus_ref.nnz(); ++k) {
            check_exact(ybus_flagged.admittance()[k], ybus_ref.admittance()[k]);
        }
    }
}

TEST_CASE("Test counting_sort_element") {
    using math_solver::counting_sort_element;
    using math_solver::YBusElementMap;
    using enum power_grid_model::YBusElementType;

    SUBCASE("Test basic sorting") {
        // Create test data: elements at various matrix positions
        std::vector<YBusElementMap> vec = {
            {.pos = {2, 1}, .element = {.element_type = bft, .idx = 5}},   // pos (2,1)
            {.pos = {0, 0}, .element = {.element_type = bff, .idx = 0}},   // pos (0,0)
            {.pos = {1, 2}, .element = {.element_type = btf, .idx = 3}},   // pos (1,2)
            {.pos = {0, 1}, .element = {.element_type = bft, .idx = 1}},   // pos (0,1)
            {.pos = {2, 1}, .element = {.element_type = shunt, .idx = 6}}, // pos (2,1) - same position as first
            {.pos = {1, 0}, .element = {.element_type = btf, .idx = 2}},   // pos (1,0)
            {.pos = {2, 2}, .element = {.element_type = btt, .idx = 7}},   // pos (2,2)
        };

        // Expected sorted order: by row first, then by column
        // (0,0), (0,1), (1,0), (1,2), (2,1), (2,1), (2,2)
        std::vector<std::pair<Idx, Idx>> expected_positions = {{0, 0}, {0, 1}, {1, 0}, {1, 2}, {2, 1}, {2, 1}, {2, 2}};

        Idx const n_bus = 3;
        counting_sort_element(vec, n_bus);

        // Verify sorting
        CHECK(vec.size() == 7);
        for (size_t i = 0; i < vec.size(); ++i) {
            CHECK(vec[i].pos.first == expected_positions[i].first);
            CHECK(vec[i].pos.second == expected_positions[i].second);
        }

        // Verify specific elements are preserved correctly
        CHECK(vec[0].element.element_type == bff);
        CHECK(vec[0].element.idx == 0);
        CHECK(vec[1].element.element_type == bft);
        CHECK(vec[1].element.idx == 1);
    }

    SUBCASE("Test with single bus") {
        std::vector<YBusElementMap> vec = {{.pos = {0, 0}, .element = {.element_type = shunt, .idx = 10}}};

        counting_sort_element(vec, 1);

        CHECK(vec.size() == 1);
        CHECK(vec[0].pos == std::make_pair(0, 0));
        CHECK(vec[0].element.element_type == shunt);
        CHECK(vec[0].element.idx == 10);
    }

    SUBCASE("Test with empty vector") {
        std::vector<YBusElementMap> vec;
        counting_sort_element(vec, 5);
        CHECK(vec.empty());
    }

    SUBCASE("Test stability - elements with same position maintain relative order") {
        std::vector<YBusElementMap> vec = {
            {.pos = {1, 1}, .element = {.element_type = bff, .idx = 100}},
            {.pos = {1, 1}, .element = {.element_type = bft, .idx = 200}},
            {.pos = {1, 1}, .element = {.element_type = shunt, .idx = 300}},
        };

        counting_sort_element(vec, 2);

        CHECK(vec.size() == 3);
        // All should be at position (1,1)
        for (const auto& element : vec) {
            CHECK(element.pos == std::make_pair(1, 1));
        }
        // Original relative order should be preserved (stable sort)
        CHECK(vec[0].element.idx == 100);
        CHECK(vec[1].element.idx == 200);
        CHECK(vec[2].element.idx == 300);
    }

    SUBCASE("Test large sparse matrix scenario") {
        std::vector<YBusElementMap> vec;
        Idx const n_bus = 10;

        // Add elements in reverse order to test sorting thoroughly
        for (Idx const i : IdxRange{n_bus * n_bus} | std::views::reverse) {
            Idx const row = i / n_bus;
            Idx const col = i % n_bus;
            if ((row + col) % 3 == 0) { // Sparse pattern
                vec.push_back({.pos = {row, col}, .element = {.element_type = bff, .idx = row * n_bus + col}});
            }
        }

        size_t original_size = vec.size();
        counting_sort_element(vec, n_bus);

        CHECK(vec.size() == original_size);

        // Verify sorted order
        for (size_t i = 1; i < vec.size(); ++i) {
            auto [prev_row, prev_col] = vec[i - 1].pos;
            auto [curr_row, curr_col] = vec[i].pos;

            // Check row-major order
            if (prev_row == curr_row) {
                CHECK(prev_col <= curr_col); // Same row, column should be non-decreasing
            } else {
                CHECK(prev_row < curr_row); // Different row, previous row should be smaller
            }
        }
    }

    SUBCASE("Test all YBusElementType values") {
        std::vector<YBusElementMap> vec = {
            {.pos = {1, 1}, .element = {.element_type = fill_in_tf, .idx = 6}},
            {.pos = {0, 1}, .element = {.element_type = bft, .idx = 1}},
            {.pos = {1, 0}, .element = {.element_type = btf, .idx = 2}},
            {.pos = {0, 0}, .element = {.element_type = bff, .idx = 0}},
            {.pos = {1, 1}, .element = {.element_type = btt, .idx = 3}},
            {.pos = {2, 2}, .element = {.element_type = shunt, .idx = 4}},
            {.pos = {1, 2}, .element = {.element_type = fill_in_ft, .idx = 5}},
        };

        counting_sort_element(vec, 3);

        // Verify positions are sorted correctly
        std::vector<std::pair<Idx, Idx>> expected_positions = {{0, 0}, {0, 1}, {1, 0}, {1, 1}, {1, 1}, {1, 2}, {2, 2}};

        CHECK(vec.size() == 7);
        for (size_t i = 0; i < vec.size(); ++i) {
            CHECK(vec[i].pos == expected_positions[i]);
        }
    }
}

} // namespace power_grid_model
