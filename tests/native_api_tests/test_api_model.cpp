// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include "load_dataset.hpp"

#include <power_grid_model_cpp/basics.hpp>
#include <power_grid_model_cpp/buffer.hpp>
#include <power_grid_model_cpp/dataset.hpp>
#include <power_grid_model_cpp/handle.hpp>
#include <power_grid_model_cpp/model.hpp>
#include <power_grid_model_cpp/options.hpp>

#include <power_grid_model_c/basics.h>
#include <power_grid_model_c/dataset_definitions.h>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <exception> // NOLINT(misc-include-cleaner)
#include <initializer_list>
#include <map>
#include <numbers>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

/*
Testing network-

source_1 -- node_0 |---- line_5 ----| node_4
              |    |---- line_6 ----|
              |
              load_2

source_1: 1.0 p.u., 100.0 V, internal_impedance(j10.0 ohm, sk=1000.0 VA, rx_ratio=0.0)
node_0: 100.0 V
oad_2: const_i, -j5.0A, 0.0 W, 500.0 var

update_0 voltage calculation:
    u_ref = 0.5 p.u. (50.0 V)
    q_specified = 100 var (-j1.0A)
u0 = 50.0 V - (j10.0 ohm * -j1.0 A) = 40.0 V

update_1 voltage calculation:
    q_specified = 300 var (-j3.0A)
u0 = 100.0 V - (j10.0 ohm * -j3.0 A) = 70.0 V

Dataset created with the following buffers:

|                | Row Based | Columnar | Dense | Sparse |
|----------------|-----------|----------|-------|--------|
| input data     |           |          |       |        |
| - node         |           |    Y     |       |        |
| - line         |           |    Y     |       |        |
| - load         |     Y     |          |       |        |
| - source       |     Y     |          |       |        |
| single update  |           |          |       |        |
| - line         |           |    Y     |       |        |
| - load         |     Y     |          |       |        |
| - source       |     Y     |          |       |        |
| batch update   |           |          |       |        |
| - line         |           |    Y     |   Y   |        |
| - load         |     Y     |          |   Y   |        |
| - source       |     Y     |          |       |    Y   |

*/

namespace power_grid_model_cpp {
namespace {
using namespace std::string_literals;
using power_grid_model_cpp_test::load_dataset;

enum class MeasuredTerminalType : IntS {
    branch_from = 0,
    branch_to = 1,
    source = 2,
    shunt = 3,
    load = 4,
    generator = 5,
    branch3_1 = 6,
    branch3_2 = 7,
    branch3_3 = 8,
    node = 9
};

void check_exception(PowerGridError const& e, PGM_ErrorCode const& reference_error,
                     std::string_view reference_err_msg) {
    CHECK(e.error_code() == reference_error);
    std::string const err_msg{e.what()};
    REQUIRE(std::in_range<doctest::String::size_type>(reference_err_msg.size()));
    doctest::String const ref_err_msg{reference_err_msg.data(),
                                      static_cast<doctest::String::size_type>(reference_err_msg.size())};
    REQUIRE(err_msg.c_str() == doctest::Contains(ref_err_msg));
}

// we pass functor by value
// it should be cheap to copy
template <typename Func, class... Args>
void check_throws_with(Func func, PGM_ErrorCode const& reference_error, std::string_view reference_err_msg,
                       Args&&... args) {
    try {
        func(std::forward<Args>(args)...);
        FAIL("Expected error not thrown.");
    } catch (PowerGridError const& e) {
        check_exception(e, reference_error, reference_err_msg);
    }
}

auto complete_state_json() {
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "sym_load": [
      {"id": 2, "node": 0, "status": 1, "type": 2, "p_specified": 0, "q_specified": 500}
    ],
    "source": [
      {"id": 1, "node": 0, "status": 1, "u_ref": 1, "sk": 1000, "rx_ratio": 0}
    ],
    "node": [
      {"id": 0, "u_rated": 100},
      {"id": 4, "u_rated": 100}
    ],
    "line": [
      {"id": 5, "from_node": 0, "to_node": 4, "from_status": 0, "to_status": 1},
      {"id": 6, "from_node": 4, "to_node": 0, "from_status": 0, "to_status": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122
}

auto single_update_json() {
    return R"json({
  "version": "1.0",
  "type": "update",
  "is_batch": false,
  "attributes": {},
  "data": {
    "source": [
      {"id": 1, "u_ref": 0.5}
    ],
    "sym_load": [
      {"id": 2, "q_specified": 100}
    ],
    "line": [
      {"id": 5, "from_status": 0, "to_status": 1},
      {"id": 6, "from_status": 0, "to_status": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122
}

auto batch_update_json() {
    return R"json({
  "version": "1.0",
  "type": "update",
  "is_batch": true,
  "attributes": {},
  "data": [
    {
      "source": [
        {"id": 1, "u_ref": 0.5}
      ],
      "sym_load": [
        {"id": 2, "q_specified": 100}
      ],
      "line": [
        {"id": 5, "from_status": 0, "to_status": 1},
        {"id": 6, "from_status": 0, "to_status": 0}
      ]
    },
    {
      "sym_load": [
        {"id": 2, "q_specified": 300}
      ],
      "line": [
        {"id": 5, "from_status": 0, "to_status": 0},
        {"id": 6, "from_status": 0, "to_status": 0}
      ]
    }
  ]
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122
}
} // namespace

TEST_CASE("API Model") {
    using namespace std::string_literals;

    Options options{};

    auto const owning_input_dataset = load_dataset(complete_state_json());
    auto const& input_dataset = owning_input_dataset.dataset;

    auto const single_owning_update_dataset = load_dataset(single_update_json());
    auto const& single_update_dataset = single_owning_update_dataset.dataset;

    auto const batch_owning_update_dataset = load_dataset(batch_update_json());
    auto const& batch_update_dataset = batch_owning_update_dataset.dataset;

    // output data
    Buffer node_output{PGM_def_sym_output_node, 2};
    node_output.set_nan();
    DatasetMutable single_output_dataset{"sym_output", false, 1};
    single_output_dataset.add_buffer("node", 2, 2, nullptr, node_output);
    Buffer node_batch_output{PGM_def_sym_output_node, 4};
    node_batch_output.set_nan();
    DatasetMutable batch_output_dataset{"sym_output", true, 2};
    batch_output_dataset.add_buffer("node", 2, 4, nullptr, node_batch_output);

    std::vector<ID> node_result_id(2);
    std::vector<int8_t> node_result_energized(2);
    std::vector<double> node_result_u(2);
    std::vector<double> node_result_u_pu(2);
    std::vector<double> node_result_u_angle(2);
    std::vector<ID> batch_node_result_id(4);
    std::vector<int8_t> batch_node_result_energized(4);
    std::vector<double> batch_node_result_u(4);
    std::vector<double> batch_node_result_u_pu(4);
    std::vector<double> batch_node_result_u_angle(4);

    // create model
    Model model{50.0, input_dataset};

    SUBCASE("Test Movability") {
        Model model_dummy{std::move(model)};
        model = std::move(model_dummy);
    }
    SUBCASE("Test Copyability") {
        Model const model_dummy{std::move(model)};
        model = Model{model_dummy};
    }
    SUBCASE("Single power flow") {
        // Common checker used for all single power flow test subcases
        auto check_common_node_results = [&]() {
            node_output.get_value(PGM_def_sym_output_node_id, node_result_id.data(), -1);
            node_output.get_value(PGM_def_sym_output_node_energized, node_result_energized.data(), 0, -1);
            node_output.get_value(PGM_def_sym_output_node_u_angle, node_result_u_angle.data(), 0, 1, -1);

            CHECK(node_result_id[0] == 0);
            CHECK(node_result_energized[0] == 1);
            CHECK(node_result_u_angle[0] == doctest::Approx(0.0));
            CHECK(node_result_id[1] == 4);
            CHECK(node_result_energized[1] == 0);
            CHECK(node_result_u[1] == doctest::Approx(0.0));
            CHECK(node_result_u_pu[1] == doctest::Approx(0.0));
            CHECK(node_result_u_angle[1] == doctest::Approx(0.0));
        };

        SUBCASE("Single power flow") {
            model.calculate(options, single_output_dataset);
            node_output.get_value(PGM_def_sym_output_node_u, node_result_u.data(), -1);
            node_output.get_value(PGM_def_sym_output_node_u_pu, node_result_u_pu.data(), -1);
            CHECK(node_result_u[0] == doctest::Approx(50.0));
            CHECK(node_result_u_pu[0] == doctest::Approx(0.5));
            check_common_node_results();
        }

        SUBCASE("Permanent update") {
            model.update(single_update_dataset);
            model.calculate(options, single_output_dataset);
            node_output.get_value(PGM_def_sym_output_node_u, node_result_u.data(), -1);
            node_output.get_value(PGM_def_sym_output_node_u_pu, node_result_u_pu.data(), -1);
            CHECK(node_result_u[0] == doctest::Approx(40.0));
            CHECK(node_result_u_pu[0] == doctest::Approx(0.4));
            check_common_node_results();
        }

        SUBCASE("Copy model") {
            auto model_copy = Model{model};
            model_copy.calculate(options, single_output_dataset);
            node_output.get_value(PGM_def_sym_output_node_u, node_result_u.data(), -1);
            node_output.get_value(PGM_def_sym_output_node_u_pu, node_result_u_pu.data(), -1);
            CHECK(node_result_u[0] == doctest::Approx(50.0));
            CHECK(node_result_u_pu[0] == doctest::Approx(0.5));
            check_common_node_results();
        }
    }

    SUBCASE("Test get indexer") {
        std::vector<ID> const node_id_2{1, 2, 3};
        std::vector<double> const node_u_rated_2{10.0e3, 10.0e3, 10.0e3};

        DatasetConst input_dataset_2{"input", false, 1};
        input_dataset_2.add_buffer("node", std::ssize(node_id_2), std::ssize(node_id_2), nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("node", "id", node_id_2.data());
        input_dataset_2.add_attribute_buffer("node", "u_rated", node_u_rated_2.data());

        auto model_2 = Model{50.0, input_dataset_2};

        SUBCASE("Good weather") {
            std::vector<ID> const ids_to_index{2, 1, 3, 2};
            std::vector<Idx> const expected_indexer{1, 0, 2, 1};
            std::vector<Idx> indexer(ids_to_index.size());
            model_2.get_indexer("node", std::ssize(ids_to_index), ids_to_index.data(), indexer.data());
            CHECK(indexer == expected_indexer);
        }
        SUBCASE("Bad weather: wrong id") {
            std::vector<ID> const ids_to_index{2, 1, 3, 4};
            std::vector<Idx> indexer(ids_to_index.size());
            CHECK_THROWS_WITH_AS(
                model_2.get_indexer("node", std::ssize(ids_to_index), ids_to_index.data(), indexer.data()),
                doctest::Contains("The id cannot be found: 4"), PowerGridRegularError);
        }
        SUBCASE("Bad weather: wrong type") {
            std::vector<ID> const ids_to_index{2, 1, 3, 2};
            std::vector<Idx> indexer(ids_to_index.size());
            CHECK_THROWS_WITH_AS(
                model_2.get_indexer("sym_load", std::ssize(ids_to_index), ids_to_index.data(), indexer.data()),
                doctest::Contains("Wrong type for object with id 2"), PowerGridRegularError);
        }
    }

    SUBCASE("Batch power flow") {
        model.calculate(options, batch_output_dataset, batch_update_dataset);
        node_batch_output.get_value(PGM_def_sym_output_node_id, batch_node_result_id.data(), -1);
        node_batch_output.get_value(PGM_def_sym_output_node_energized, batch_node_result_energized.data(), -1);
        node_batch_output.get_value(PGM_def_sym_output_node_u, batch_node_result_u.data(), -1);
        node_batch_output.get_value(PGM_def_sym_output_node_u_pu, batch_node_result_u_pu.data(), -1);
        node_batch_output.get_value(PGM_def_sym_output_node_u_angle, batch_node_result_u_angle.data(), -1);
        CHECK(batch_node_result_id[0] == 0);
        CHECK(batch_node_result_energized[0] == 1);
        CHECK(batch_node_result_u[0] == doctest::Approx(40.0));
        CHECK(batch_node_result_u_pu[0] == doctest::Approx(0.4));
        CHECK(batch_node_result_u_angle[0] == doctest::Approx(0.0));
        CHECK(batch_node_result_id[1] == 4);
        CHECK(batch_node_result_energized[1] == 0);
        CHECK(batch_node_result_u[1] == doctest::Approx(0.0));
        CHECK(batch_node_result_u_pu[1] == doctest::Approx(0.0));
        CHECK(batch_node_result_u_angle[1] == doctest::Approx(0.0));
        CHECK(batch_node_result_id[2] == 0);
        CHECK(batch_node_result_energized[2] == 1);
        CHECK(batch_node_result_u[2] == doctest::Approx(70.0));
        CHECK(batch_node_result_u_pu[2] == doctest::Approx(0.7));
        CHECK(batch_node_result_u_angle[2] == doctest::Approx(0.0));
        CHECK(batch_node_result_id[3] == 4);
        CHECK(batch_node_result_energized[3] == 0);
        CHECK(batch_node_result_u[3] == doctest::Approx(0.0));
        CHECK(batch_node_result_u_pu[3] == doctest::Approx(0.0));
        CHECK(batch_node_result_u_angle[3] == doctest::Approx(0.0));
    }

    SUBCASE("Input error handling") {
        SUBCASE("Construction error") {
            auto const bad_load_id_state_json = R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "sym_load": [
      {"id": 0, "node": 0, "status": 1, "type": 2, "p_specified": 0, "q_specified": 500}
    ],
    "source": [
      {"id": 1, "node": 0, "status": 1, "u_ref": 1, "sk": 1000, "rx_ratio": 0}
    ],
    "node": [
      {"id": 0, "u_rated": 100},
      {"id": 4, "u_rated": 100}
    ],
    "line": [
      {"id": 5, "from_node": 0, "to_node": 4, "from_status": 0, "to_status": 1},
      {"id": 6, "from_node": 4, "to_node": 0, "from_status": 0, "to_status": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

            auto const bad_owning_input_dataset = load_dataset(bad_load_id_state_json);
            auto const& bad_input_dataset = bad_owning_input_dataset.dataset;

            auto const bad_model_lambda = [&bad_input_dataset]() { Model const wrong_model{50.0, bad_input_dataset}; };

            check_throws_with(bad_model_lambda, PGM_regular_error, "Conflicting id detected:"s);
        }

        SUBCASE("Update error") {
            auto const bad_source_id_single_update_json = R"json({
  "version": "1.0",
  "type": "update",
  "is_batch": false,
  "attributes": {},
  "data": {
    "source": [
      {"id": 99, "u_ref": 0.5}
    ],
    "sym_load": [
      {"id": 2, "q_specified": 100}
    ],
    "line": [
      {"id": 5, "from_status": 0, "to_status": 1},
      {"id": 6, "from_status": 0, "to_status": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

            auto const bad_single_owning_update_dataset = load_dataset(bad_source_id_single_update_json);
            auto const bad_update_lambda = [&model, &bad_single_owning_update_dataset]() {
                model.update(bad_single_owning_update_dataset.dataset);
            };

            check_throws_with(bad_update_lambda, PGM_regular_error, "The id cannot be found:"s);
        }

        SUBCASE("Update error in calculation") {
            auto const bad_load_id_batch_update_json = R"json({
  "version": "1.0",
  "type": "update",
  "is_batch": true,
  "attributes": {},
  "data": [
    {
      "source": [
        {"id": 1, "u_ref": 0.5}
      ],
      "sym_load": [
        {"id": 2, "q_specified": 100}
      ],
      "line": [
        {"id": 99, "from_status": 0, "to_status": 1},
        {"id": 999, "from_status": 0, "to_status": 0}
      ]
    },
    {
      "sym_load": [
        {"id": 2, "q_specified": 300}
      ],
      "line": [
        {"id": 9999, "from_status": 0, "to_status": 0},
        {"id": 99999, "from_status": 0, "to_status": 0}
      ]
    }
  ]
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

            auto const bad_batch_owning_update_dataset = load_dataset(bad_load_id_batch_update_json);
            auto const bad_calc_with_update_lambda = [&model, &options, &batch_output_dataset,
                                                      &bad_batch_owning_update_dataset]() {
                model.calculate(options, batch_output_dataset, bad_batch_owning_update_dataset.dataset);
            };
            check_throws_with(bad_calc_with_update_lambda, PGM_batch_error, "The id cannot be found:"s);
        }

        SUBCASE("Invalid calculation type error") {
            auto const bad_calc_type_lambda = [&options, &model, &single_output_dataset]() {
                options.set_calculation_type(-128);
                model.calculate(options, single_output_dataset);
            };
            check_throws_with(bad_calc_type_lambda, PGM_regular_error, "CalculationType is not implemented for"s);
        }

        SUBCASE("Invalid tap changing strategy error") {
            auto const bad_tap_strat_lambda = [&options, &model, &single_output_dataset]() {
                options.set_tap_changing_strategy(-128);
                model.calculate(options, single_output_dataset);
            };
            check_throws_with(bad_tap_strat_lambda, PGM_regular_error, "get_optimizer_type is not implemented for"s);
        }

        SUBCASE("Tap changing strategy") {
            options.set_tap_changing_strategy(PGM_tap_changing_strategy_min_voltage_tap);
            CHECK_NOTHROW(model.calculate(options, single_output_dataset));
        }
    }

    SUBCASE("Calculation error") {
        SUBCASE("Single calculation error") {
            // not converging
            options.set_max_iter(1);
            options.set_err_tol(1e-100);
            options.set_symmetric(0);
            options.set_threading(1);
            auto const calc_error_lambda = [&model, &single_output_dataset](auto const& opt) {
                model.calculate(opt, single_output_dataset);
            };
            check_throws_with(calc_error_lambda, PGM_regular_error, "Iteration failed to converge after"s, options);

            // wrong method
            options.set_calculation_type(PGM_state_estimation);
            options.set_calculation_method(PGM_iterative_current);
            check_throws_with(calc_error_lambda, PGM_regular_error,
                              "The calculation method is invalid for this calculation!"s, options);
        }

        SUBCASE("Batch calculation error") {
            SUBCASE("Line bad line id") {
                auto const bad_line_id_batch_update_json = R"json({
  "version": "1.0",
  "type": "update",
  "is_batch": true,
  "attributes": {},
  "data": [
    {
      "source": [
        {"id": 1, "u_ref": 0.5}
      ],
      "sym_load": [
        {"id": 2, "q_specified": 100}
      ],
      "line": [
        {"id": 5, "from_status": 0, "to_status": 1},
        {"id": 6, "from_status": 0, "to_status": 0}
      ]
    },
    {
      "sym_load": [
        {"id": 999, "q_specified": 300}
      ],
      "line": [
        {"id": 5, "from_status": 0, "to_status": 0},
        {"id": 6, "from_status": 0, "to_status": 0}
      ]
    }
  ]
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

                auto const bad_batch_owning_update_dataset = load_dataset(bad_line_id_batch_update_json);

                // failed in batch scenario 1
                try {
                    model.calculate(options, batch_output_dataset, bad_batch_owning_update_dataset.dataset);
                    FAIL("Expected batch calculation error not thrown.");
                } catch (PowerGridBatchError const& e) {
                    CHECK(e.error_code() == PGM_batch_error);
                    auto const& failed_scenarios = e.failed_scenarios();
                    CHECK(failed_scenarios.size() == 1);
                    CHECK(failed_scenarios[0].scenario == 1);
                    std::string const err_msg{failed_scenarios[0].error_message};
                    CHECK(err_msg.find("The id cannot be found:"s) != std::string::npos);
                }
                // valid results for batch 0
                node_batch_output.get_value(PGM_def_sym_output_node_id, batch_node_result_id.data(), -1);
                node_batch_output.get_value(PGM_def_sym_output_node_energized, batch_node_result_energized.data(), -1);
                node_batch_output.get_value(PGM_def_sym_output_node_u, batch_node_result_u.data(), -1);
                node_batch_output.get_value(PGM_def_sym_output_node_u_pu, batch_node_result_u_pu.data(), -1);
                node_batch_output.get_value(PGM_def_sym_output_node_u_angle, batch_node_result_u_angle.data(), -1);
                CHECK(batch_node_result_id[0] == 0);
                CHECK(batch_node_result_energized[0] == 1);
                CHECK(batch_node_result_u[0] == doctest::Approx(40.0));
                CHECK(batch_node_result_u_pu[0] == doctest::Approx(0.4));
                CHECK(batch_node_result_u_angle[0] == doctest::Approx(0.0));
                CHECK(batch_node_result_id[1] == 4);
                CHECK(batch_node_result_energized[1] == 0);
                CHECK(batch_node_result_u[1] == doctest::Approx(0.0));
                CHECK(batch_node_result_u_pu[1] == doctest::Approx(0.0));
                CHECK(batch_node_result_u_angle[1] == doctest::Approx(0.0));
            }
        }
    }

    SUBCASE("Model update optional id") {
        std::vector<ID> const input_node_id{0};
        std::vector<double> const input_node_u_rated{100.0};
        Buffer input_node_buffer{PGM_def_input_node, 1};
        input_node_buffer.set_nan();
        input_node_buffer.set_value(PGM_def_input_node_id, input_node_id.data(), -1);
        input_node_buffer.set_value(PGM_def_input_node_u_rated, input_node_u_rated.data(), -1);

        std::vector<ID> const input_source_id{1};
        std::vector<ID> const input_source_node{0};
        std::vector<int8_t> const input_source_status{1};
        std::vector<double> const input_source_u_ref{1.0};
        std::vector<double> const input_source_sk{1000.0};
        std::vector<double> const input_source_rx_ratio{0.0};
        Buffer input_source_buffer{PGM_def_input_source, 1};
        input_source_buffer.set_nan();
        input_source_buffer.set_value(PGM_def_input_source_id, input_source_id.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_node, input_source_node.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_status, input_source_status.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_u_ref, input_source_u_ref.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_sk, input_source_sk.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_rx_ratio, input_source_rx_ratio.data(), -1);

        std::vector<ID> const input_sym_load_id{2};
        std::vector<ID> const input_sym_load_node{0};
        std::vector<int8_t> const input_sym_load_status{1};
        std::vector<int8_t> const input_sym_load_type{2};
        std::vector<double> const input_sym_load_p_specified{0.0};
        std::vector<double> const input_sym_load_q_specified{500.0};
        Buffer input_sym_load_buffer{PGM_def_input_sym_load, 1};
        input_sym_load_buffer.set_nan();
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_id, input_sym_load_id.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_node, input_sym_load_node.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_status, input_sym_load_status.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_type, input_sym_load_type.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_p_specified, input_sym_load_p_specified.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_q_specified, input_sym_load_q_specified.data(), -1);

        // input dataset - row
        DatasetConst input_dataset_row{"input", false, 1};
        input_dataset_row.add_buffer("node", 1, 1, nullptr, input_node_buffer);
        input_dataset_row.add_buffer("source", 1, 1, nullptr, input_source_buffer);
        input_dataset_row.add_buffer("sym_load", 1, 1, nullptr, input_sym_load_buffer);

        // input dataset - col
        DatasetConst input_dataset_col{"input", false, 1};
        input_dataset_col.add_buffer("node", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("node", "id", input_node_id.data());
        input_dataset_col.add_attribute_buffer("node", "u_rated", input_node_u_rated.data());

        input_dataset_col.add_buffer("source", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("source", "id", input_source_id.data());
        input_dataset_col.add_attribute_buffer("source", "node", input_source_node.data());
        input_dataset_col.add_attribute_buffer("source", "status", input_source_status.data());
        input_dataset_col.add_attribute_buffer("source", "u_ref", input_source_u_ref.data());
        input_dataset_col.add_attribute_buffer("source", "sk", input_source_sk.data());
        input_dataset_col.add_attribute_buffer("source", "rx_ratio", input_source_rx_ratio.data());

        input_dataset_col.add_buffer("sym_load", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("sym_load", "id", input_sym_load_id.data());
        input_dataset_col.add_attribute_buffer("sym_load", "node", input_sym_load_node.data());
        input_dataset_col.add_attribute_buffer("sym_load", "status", input_sym_load_status.data());
        input_dataset_col.add_attribute_buffer("sym_load", "type", input_sym_load_type.data());
        input_dataset_col.add_attribute_buffer("sym_load", "p_specified", input_sym_load_p_specified.data());
        input_dataset_col.add_attribute_buffer("sym_load", "q_specified", input_sym_load_q_specified.data());

        // update dataset
        std::vector<Idx> update_source_indptr{0, 1, 2};
        std::vector<ID> const update_source_id{1, 1};
        std::vector<double> const update_source_u_ref{0.5, 1.0};
        Buffer update_source_buffer{PGM_def_update_source, 2};
        update_source_buffer.set_nan();
        update_source_buffer.set_value(PGM_def_update_source_id, update_source_id.data(), -1);
        update_source_buffer.set_value(PGM_def_update_source_u_ref, update_source_u_ref.data(), -1);

        std::vector<Idx> update_sym_load_indptr{0, 1, 2};
        std::vector<ID> const update_sym_load_id{2, 5};
        std::vector<double> const update_sym_load_q_specified{100.0, 300.0};
        Buffer update_sym_load_buffer{PGM_def_update_sym_load, 2};
        update_sym_load_buffer.set_nan();
        update_sym_load_buffer.set_value(PGM_def_update_sym_load_id, update_sym_load_id.data(), -1);
        update_sym_load_buffer.set_value(PGM_def_update_sym_load_q_specified, update_sym_load_q_specified.data(), -1);

        // update dataset buffers - no ids
        Buffer update_source_buffer_no_id{PGM_def_update_source, 2};
        update_source_buffer_no_id.set_nan();
        update_source_buffer_no_id.set_value(PGM_def_update_source_u_ref, update_source_u_ref.data(), -1);

        Buffer update_sym_load_buffer_no_id{PGM_def_update_sym_load, 2};
        update_sym_load_buffer_no_id.set_nan();
        update_sym_load_buffer_no_id.set_value(PGM_def_update_sym_load_q_specified, update_sym_load_q_specified.data(),
                                               -1);

        // update dataset - row
        DatasetConst update_dataset_row{"update", true, 2};
        update_dataset_row.add_buffer("source", -1, 2, update_source_indptr.data(), update_source_buffer);
        update_dataset_row.add_buffer("sym_load", -1, 2, update_sym_load_indptr.data(), update_sym_load_buffer);

        // update dataset - col
        DatasetConst update_dataset_col{"update", true, 2};

        update_dataset_col.add_buffer("source", -1, 2, update_source_indptr.data(), nullptr);
        update_dataset_col.add_attribute_buffer("source", "id", update_source_id.data());
        update_dataset_col.add_attribute_buffer("source", "u_ref", update_source_u_ref.data());

        update_dataset_col.add_buffer("sym_load", -1, 2, update_sym_load_indptr.data(), nullptr);
        update_dataset_col.add_attribute_buffer("sym_load", "id", update_sym_load_id.data());
        update_dataset_col.add_attribute_buffer("sym_load", "q_specified", update_sym_load_q_specified.data());

        // update dataset - row no ids
        DatasetConst update_dataset_row_no_id{"update", true, 2};
        update_dataset_row_no_id.add_buffer("source", -1, 2, update_source_indptr.data(), update_source_buffer_no_id);
        update_dataset_row_no_id.add_buffer("sym_load", -1, 2, update_sym_load_indptr.data(),
                                            update_sym_load_buffer_no_id);

        // update dataset - col no ids
        DatasetConst update_dataset_col_no_id{"update", true, 2};
        update_dataset_col_no_id.add_buffer("source", -1, 2, update_source_indptr.data(), nullptr);

        update_dataset_col_no_id.add_attribute_buffer("source", "u_ref", update_source_u_ref.data());

        update_dataset_col_no_id.add_buffer("sym_load", -1, 2, update_sym_load_indptr.data(), nullptr);
        update_dataset_col_no_id.add_attribute_buffer("sym_load", "q_specified", update_sym_load_q_specified.data());

        // output data
        Buffer batch_node_output{PGM_def_sym_output_node, 2};
        batch_node_output.set_nan();
        DatasetMutable batch_output{"sym_output", true, 2};
        batch_output.add_buffer("node", 1, 2, nullptr, batch_node_output);

        // options
        Options const batch_options{};

        SUBCASE("Row-based input dataset") {
            Model row_model{50.0, input_dataset_row};

            SUBCASE("Row-based update dataset error") {
                CHECK_THROWS_AS(row_model.calculate(batch_options, batch_output, update_dataset_row),
                                PowerGridBatchError);
            }
            SUBCASE("Row-based update dataset wo id") {
                row_model.calculate(batch_options, batch_output_dataset, update_dataset_row_no_id);
            }
            SUBCASE("Columnar update dataset error") {
                CHECK_THROWS_AS(row_model.calculate(batch_options, batch_output, update_dataset_col),
                                PowerGridBatchError);
            }
            SUBCASE("Columnar update dataset wo id") {
                row_model.calculate(batch_options, batch_output, update_dataset_col_no_id);
            }
            SUBCASE("Columnar update dataset wo id - non-uniform") {
                update_source_indptr = {0, 1, 1};
                CHECK_THROWS_AS(row_model.calculate(batch_options, batch_output, update_dataset_col_no_id),
                                PowerGridBatchError);
            }
        }

        SUBCASE("Columnar input dataset") {
            Model col_model{50.0, input_dataset_col};

            SUBCASE("Row-based update dataset error") {
                CHECK_THROWS_AS(col_model.calculate(batch_options, batch_output, update_dataset_row),
                                PowerGridBatchError);
            }
            SUBCASE("Row-based update dataset wo id") {
                col_model.calculate(batch_options, batch_output_dataset, update_dataset_row_no_id);
            }
            SUBCASE("Columnar update dataset error") {
                CHECK_THROWS_AS(col_model.calculate(batch_options, batch_output, update_dataset_col),
                                PowerGridBatchError);
            }
            SUBCASE("Columnar update dataset wo id") {
                col_model.calculate(batch_options, batch_output, update_dataset_col_no_id);
            }
            SUBCASE("Columnar update dataset wo id - non-uniform") {
                update_source_indptr = {0, 1, 1};
                CHECK_THROWS_AS(col_model.calculate(batch_options, batch_output, update_dataset_col_no_id),
                                PowerGridBatchError);
            }
        }
    }

    SUBCASE("Self contained model update error") {
        std::vector<ID> const input_node_id{0};
        std::vector<double> const input_node_u_rated{100.0};
        Buffer input_node_buffer{PGM_def_input_node, 1};
        input_node_buffer.set_nan();
        input_node_buffer.set_value(PGM_def_input_node_id, input_node_id.data(), -1);
        input_node_buffer.set_value(PGM_def_input_node_u_rated, input_node_u_rated.data(), -1);

        std::vector<ID> const input_source_id{1};
        std::vector<ID> const input_source_node{0};
        std::vector<int8_t> const input_source_status{1};
        std::vector<double> const input_source_u_ref{1.0};
        std::vector<double> const input_source_sk{1000.0};
        std::vector<double> const input_source_rx_ratio{0.0};
        Buffer input_source_buffer{PGM_def_input_source, 1};
        input_source_buffer.set_nan();
        input_source_buffer.set_value(PGM_def_input_source_id, input_source_id.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_node, input_source_node.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_status, input_source_status.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_u_ref, input_source_u_ref.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_sk, input_source_sk.data(), -1);
        input_source_buffer.set_value(PGM_def_input_source_rx_ratio, input_source_rx_ratio.data(), -1);

        std::vector<ID> const input_sym_load_id{2};
        std::vector<ID> const input_sym_load_node{0};
        std::vector<int8_t> const input_sym_load_status{1};
        std::vector<int8_t> const input_sym_load_type{2};
        std::vector<double> const input_sym_load_p_specified{0.0};
        std::vector<double> const input_sym_load_q_specified{500.0};
        Buffer input_sym_load_buffer{PGM_def_input_sym_load, 1};
        input_sym_load_buffer.set_nan();
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_id, input_sym_load_id.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_node, input_sym_load_node.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_status, input_sym_load_status.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_type, input_sym_load_type.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_p_specified, input_sym_load_p_specified.data(), -1);
        input_sym_load_buffer.set_value(PGM_def_input_sym_load_q_specified, input_sym_load_q_specified.data(), -1);

        // input dataset - row
        DatasetConst input_dataset_row{"input", false, 1};
        input_dataset_row.add_buffer("node", 1, 1, nullptr, input_node_buffer);
        input_dataset_row.add_buffer("source", 1, 1, nullptr, input_source_buffer);
        input_dataset_row.add_buffer("sym_load", 1, 1, nullptr, input_sym_load_buffer);

        // input dataset - col
        DatasetConst input_dataset_col{"input", false, 1};
        input_dataset_col.add_buffer("node", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("node", "id", input_node_id.data());
        input_dataset_col.add_attribute_buffer("node", "u_rated", input_node_u_rated.data());

        input_dataset_col.add_buffer("source", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("source", "id", input_source_id.data());
        input_dataset_col.add_attribute_buffer("source", "node", input_source_node.data());
        input_dataset_col.add_attribute_buffer("source", "status", input_source_status.data());
        input_dataset_col.add_attribute_buffer("source", "u_ref", input_source_u_ref.data());
        input_dataset_col.add_attribute_buffer("source", "sk", input_source_sk.data());
        input_dataset_col.add_attribute_buffer("source", "rx_ratio", input_source_rx_ratio.data());

        input_dataset_col.add_buffer("sym_load", 1, 1, nullptr, nullptr);
        input_dataset_col.add_attribute_buffer("sym_load", "id", input_sym_load_id.data());
        input_dataset_col.add_attribute_buffer("sym_load", "node", input_sym_load_node.data());
        input_dataset_col.add_attribute_buffer("sym_load", "status", input_sym_load_status.data());
        input_dataset_col.add_attribute_buffer("sym_load", "type", input_sym_load_type.data());
        input_dataset_col.add_attribute_buffer("sym_load", "p_specified", input_sym_load_p_specified.data());
        input_dataset_col.add_attribute_buffer("sym_load", "q_specified", input_sym_load_q_specified.data());

        // update dataset
        std::vector<Idx> source_indptr{0, 1, 1};
        std::vector<ID> const update_source_id{1};
        std::vector<double> const update_source_u_ref{0.5};
        Buffer update_source_buffer{PGM_def_update_source, 1};
        update_source_buffer.set_nan();
        update_source_buffer.set_value(PGM_def_update_source_id, update_source_id.data(), -1);
        update_source_buffer.set_value(PGM_def_update_source_u_ref, update_source_u_ref.data(), -1);

        std::vector<Idx> sym_load_indptr{0, 1, 2};
        std::vector<ID> const update_sym_load_id{2, 5};
        std::vector<double> const update_sym_load_q_specified{100.0, 300.0};
        Buffer update_sym_load_buffer{PGM_def_update_sym_load, 2};
        update_sym_load_buffer.set_nan();
        update_sym_load_buffer.set_value(PGM_def_update_sym_load_id, update_sym_load_id.data(), -1);
        update_sym_load_buffer.set_value(PGM_def_update_sym_load_q_specified, update_sym_load_q_specified.data(), -1);

        // update dataset - row
        DatasetConst update_dataset_row{"update", true, 2};
        update_dataset_row.add_buffer("source", -1, 1, source_indptr.data(), update_source_buffer);
        update_dataset_row.add_buffer("sym_load", -1, 2, sym_load_indptr.data(), update_sym_load_buffer);

        // update dataset - col
        DatasetConst update_dataset_col{"update", true, 2};

        update_dataset_col.add_buffer("source", -1, 1, source_indptr.data(), nullptr);
        update_dataset_col.add_attribute_buffer("source", "id", update_source_id.data());
        update_dataset_col.add_attribute_buffer("source", "u_ref", update_source_u_ref.data());

        update_dataset_col.add_buffer("sym_load", -1, 2, sym_load_indptr.data(), nullptr);
        update_dataset_col.add_attribute_buffer("sym_load", "id", update_sym_load_id.data());
        update_dataset_col.add_attribute_buffer("sym_load", "q_specified", update_sym_load_q_specified.data());

        // output data
        Buffer output_node_batch{PGM_def_sym_output_node, 2};
        output_node_batch.set_nan();
        DatasetMutable output_batch_dataset{"sym_output", true, 2};
        output_batch_dataset.add_buffer("node", 1, 2, nullptr, output_node_batch);

        // options
        Options const opt{};

        SUBCASE("Row-based input dataset") {
            Model row_model{50.0, input_dataset_row};

            SUBCASE("Row-based update dataset") {
                CHECK_THROWS_AS(row_model.calculate(opt, output_batch_dataset, update_dataset_row),
                                PowerGridBatchError);
            }
        }

        SUBCASE("Columnar input dataset") {
            Model col_model{50.0, input_dataset_col};

            SUBCASE("Row-based update dataset") {
                CHECK_THROWS_AS(col_model.calculate(opt, output_batch_dataset, update_dataset_row),
                                PowerGridBatchError);
            }
        }
    }

    SUBCASE("Wrong calculation methods") {
        using namespace std::string_view_literals;

        constexpr auto invalid_calculation_method_pattern = "The calculation method is invalid for this calculation!";
        constexpr auto all_types = std::array{PGM_power_flow, PGM_state_estimation, PGM_short_circuit};
        constexpr auto all_methods =
            std::array{PGM_default_method,    PGM_linear,           PGM_newton_raphson, PGM_linear_current,
                       PGM_iterative_current, PGM_iterative_linear, PGM_iec60909};

        auto supported_methods = std::map<PGM_CalculationType, std::vector<PGM_CalculationMethod>>{
            {PGM_power_flow, std::vector{PGM_default_method, PGM_newton_raphson, PGM_linear, PGM_linear_current,
                                         PGM_iterative_current}},
            {PGM_state_estimation, std::vector{PGM_default_method, PGM_iterative_linear, PGM_newton_raphson}},
            {PGM_short_circuit, std::vector{PGM_default_method, PGM_iec60909}}};

        auto output_dataset_types = std::map<PGM_CalculationType, std::string>{
            {PGM_power_flow, "sym_output"s}, {PGM_state_estimation, "sym_output"s}, {PGM_short_circuit, "sc_output"s}};

        for (auto calculation_type : all_types) {
            CAPTURE(calculation_type);
            auto const& supported_type_methods = supported_methods.at(calculation_type);

            DatasetMutable const output_dataset{output_dataset_types.at(calculation_type), false, 1};

            for (auto calculation_method : all_methods) {
                CAPTURE(calculation_method);
                options.set_calculation_type(calculation_type);
                options.set_calculation_method(calculation_method);

                if (std::ranges::find(supported_type_methods, calculation_method) == std::end(supported_type_methods)) {
                    CHECK_THROWS_WITH_AS(model.calculate(options, output_dataset), invalid_calculation_method_pattern,
                                         PowerGridRegularError);
                } else {
                    try {
                        model.calculate(options, output_dataset);
                    } catch (std::exception const& e) { // NOSONAR(S1181)
                        CHECK(e.what() != doctest::Contains(invalid_calculation_method_pattern));
                    }
                }
            }
        }
    }

    SUBCASE("Forbid link power measurements") {
        SUBCASE("SymPowerSensor") {
            auto const input_data_se_json = R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 10000},
      {"id": 2, "u_rated": 10000}
    ],
    "link": [
      {"id": 3, "from_node": 1, "to_node": 2}
    ],
    "sym_power_sensor": [
      {"id": 4, "measured_object": 3, "measured_terminal_type": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

            auto const owning_input_dataset_se = load_dataset(input_data_se_json);
            auto const& input_dataset_se = owning_input_dataset_se.dataset;

            auto const construct_model = [&input_dataset_se] { return Model{50.0, input_dataset_se}; };

            CHECK_THROWS_WITH_AS(construct_model(), "PowerSensor measurement is not supported for object of type Link",
                                 PowerGridRegularError);
        }

        SUBCASE("AsymPowerSensor") {
            auto const input_data_se_json = R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 10000},
      {"id": 2, "u_rated": 10000}
    ],
    "link": [
      {"id": 3, "from_node": 1, "to_node": 2}
    ],
    "asym_power_sensor": [
      {"id": 4, "measured_object": 3, "measured_terminal_type": 0}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

            auto const owning_input_dataset_se = load_dataset(input_data_se_json);
            auto const& input_dataset_se = owning_input_dataset_se.dataset;

            auto const construct_model = [&input_dataset_se] { return Model{50.0, input_dataset_se}; };

            CHECK_THROWS_WITH_AS(construct_model(), "PowerSensor measurement is not supported for object of type Link",
                                 PowerGridRegularError);
        }
    }

    SUBCASE("Test duplicated id") {
        auto const input_data_2_json = R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1},
      {"id": 1},
      {"id": 3}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122

        auto const owning_input_dataset_2 = load_dataset(input_data_2_json);
        auto const& input_dataset_2 = owning_input_dataset_2.dataset;

        auto construct_model = [&] { Model{50.0, input_dataset_2}; };
        CHECK_THROWS_WITH_AS(construct_model(), "Conflicting id detected: 1\n", PowerGridRegularError);
    }

    SUBCASE("Test non-existing id") {
        auto const input_data_2_json = R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 10000},
      {"id": 2, "u_rated": 10000}
    ],
    "link": [
      {"id": 5, "from_node": 99, "to_node": 2}
    ]
  }
})json"s; // NOLINT(misc-include-cleaner) https://github.com/llvm/llvm-project/issues/98122
        auto const owning_input_dataset_2 = load_dataset(input_data_2_json);
        auto const& input_dataset_2 = owning_input_dataset_2.dataset;

        auto construct_model = [&] { Model{50.0, input_dataset_2}; };
        CHECK_THROWS_WITH_AS(construct_model(), "The id cannot be found: 99\n", PowerGridRegularError);
    }

    SUBCASE("Test id for wrong type") {
        std::vector<ID> const node_id_2{1, 2, 3};
        std::vector<double> const node_u_rated_2{10.0e3, 10.0e3, 10.0e3};

        std::vector<ID> line_id_2{9};
        std::vector<ID> line_from_node_2{1};
        std::vector<ID> line_to_node_2{2};

        std::vector<ID> link_id{5};
        std::vector<ID> link_from_node{2};
        std::vector<ID> link_to_node{3};

        std::vector<ID> sym_voltage_sensor_id{25};
        std::vector<ID> sym_voltage_sensor_measured_object{1};

        std::vector<ID> sym_power_sensor_id{28};
        std::vector<ID> sym_power_sensor_measured_object{3};
        std::vector<MeasuredTerminalType> sym_power_sensor_measured_terminal_type{MeasuredTerminalType::node};

        DatasetConst input_dataset_2{"input", false, 1};

        input_dataset_2.add_buffer("node", std::ssize(node_id_2), std::ssize(node_id_2), nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("node", "id", node_id_2.data());
        input_dataset_2.add_attribute_buffer("node", "u_rated", node_u_rated_2.data());

        input_dataset_2.add_buffer("line", std::ssize(line_id_2), std::ssize(line_id_2), nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("line", "id", line_id_2.data());
        input_dataset_2.add_attribute_buffer("line", "from_node", line_from_node_2.data());
        input_dataset_2.add_attribute_buffer("line", "to_node", line_to_node_2.data());

        input_dataset_2.add_buffer("link", std::ssize(link_id), std::ssize(link_id), nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("link", "id", link_id.data());
        input_dataset_2.add_attribute_buffer("link", "from_node", link_from_node.data());
        input_dataset_2.add_attribute_buffer("link", "to_node", link_to_node.data());

        input_dataset_2.add_buffer("sym_voltage_sensor", std::ssize(sym_voltage_sensor_id),
                                   std::ssize(sym_voltage_sensor_id), nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("sym_voltage_sensor", "id", sym_voltage_sensor_id.data());
        input_dataset_2.add_attribute_buffer("sym_voltage_sensor", "measured_object",
                                             sym_voltage_sensor_measured_object.data());

        input_dataset_2.add_buffer("sym_power_sensor", std::ssize(sym_power_sensor_id), std::ssize(sym_power_sensor_id),
                                   nullptr, nullptr);
        input_dataset_2.add_attribute_buffer("sym_power_sensor", "id", sym_power_sensor_id.data());
        input_dataset_2.add_attribute_buffer("sym_power_sensor", "measured_object",
                                             sym_power_sensor_measured_object.data());
        input_dataset_2.add_attribute_buffer("sym_power_sensor", "measured_terminal_type",
                                             sym_power_sensor_measured_terminal_type.data());

        auto construct_model = [&] { Model{50.0, input_dataset_2}; };

        SUBCASE("Correct type") { CHECK_NOTHROW(construct_model()); }
        SUBCASE("Wrong branch terminal node") {
            link_from_node[0] = 9;
            CHECK_THROWS_WITH_AS(construct_model(), "Wrong type for object with id 9\n", PowerGridRegularError);
        }
        SUBCASE("Wrong voltage sensor measured object") {
            sym_voltage_sensor_measured_object[0] = 5;
            CHECK_THROWS_WITH_AS(construct_model(), "Wrong type for object with id 5\n", PowerGridRegularError);
        }
        SUBCASE("Wrong power sensor measured object type") {
            using enum MeasuredTerminalType;
            std::vector<MeasuredTerminalType> const mt_types{branch_from, branch_to, generator, load, shunt, source};

            for (auto const& mt_type : mt_types) {
                CAPTURE(mt_type);
                sym_power_sensor_measured_terminal_type[0] = mt_type;
                CHECK_THROWS_WITH_AS(construct_model(), "Wrong type for object with id 3\n", PowerGridRegularError);
            }
        }
    }
}

namespace {
// node 1 (abc, source) --a-- node 2 (a); optionally with a three-phase fault at node 1
auto reduced_phase_model_json(bool with_fault) {
    using namespace std::string_literals;
    std::string const fault = with_fault ? R"json(,
    "fault": [
      {"id": 6, "status": 1, "fault_type": 0, "fault_phase": 0, "fault_object": 1}
    ])json"s
                                         : ""s;
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 10000},
      {"id": 2, "u_rated": 10000}
    ],
    "source": [
      {"id": 3, "node": 1, "status": 1, "u_ref": 1.0, "sk": 1e8}
    ],
    "asym_line": [
      {"id": 4, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "phases": 1,
       "r_aa": 0.5, "x_aa": 0.4, "c_aa": 1e-8}
    ],
    "asym_load": [
      {"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": [1000, 0, 0], "q_specified": [0, 0, 0]}
    ])json"s +
           fault +
           R"json(
  }
})json"s;
}
} // namespace

TEST_CASE("API Model - reduced phase") {
    using namespace std::string_literals;

    auto const owning_input_dataset = load_dataset(reduced_phase_model_json(true));
    Model model{50.0, owning_input_dataset.dataset};

    Options options{};
    options.set_experimental_features(PGM_experimental_features_enabled);

    Buffer node_output{PGM_def_asym_output_node, 2};
    node_output.set_nan();
    DatasetMutable asym_output{"asym_output", false, 1};
    asym_output.add_buffer("node", 2, 2, nullptr, node_output);

    Buffer sym_node_output{PGM_def_sym_output_node, 2};
    sym_node_output.set_nan();
    DatasetMutable sym_output{"sym_output", false, 1};
    sym_output.add_buffer("node", 2, 2, nullptr, sym_node_output);

    SUBCASE("Asymmetric power flow") {
        options.set_symmetric(0);
        options.set_calculation_method(PGM_newton_raphson);
        model.calculate(options, asym_output);
        std::array<std::array<double, 3>, 2> u{};
        node_output.get_value(PGM_def_asym_output_node_u, u.data(), -1);
        CHECK(u[1][0] > 0.0);
        CHECK(u[1][1] == 0.0);
        CHECK(u[1][2] == 0.0);
    }

    SUBCASE("Experimental features are required") {
        options.set_experimental_features(PGM_experimental_features_disabled);
        options.set_symmetric(0);
        check_throws_with([&model, &asym_output](auto const& opt) { model.calculate(opt, asym_output); },
                          PGM_regular_error, "is an experimental feature"s, options);
    }

    SUBCASE("Symmetric power flow is not supported") {
        check_throws_with([&model, &sym_output](auto const& opt) { model.calculate(opt, sym_output); },
                          PGM_regular_error,
                          "Symmetric power flow is not supported in a network with one- or two-phase nodes"s, options);
    }

    SUBCASE("State estimation is not supported") {
        options.set_calculation_type(PGM_state_estimation);
        for (Idx const symmetric : {0, 1}) {
            options.set_symmetric(symmetric);
            auto& output = symmetric == 0 ? asym_output : sym_output;
            check_throws_with([&model, &output](auto const& opt) { model.calculate(opt, output); }, PGM_regular_error,
                              "State estimation is not supported in a network with one- or two-phase nodes"s, options);
        }
    }

    SUBCASE("Three-phase fault at a three-phase node") {
        options.set_calculation_type(PGM_short_circuit);
        options.set_calculation_method(PGM_iec60909);
        Buffer fault_output{PGM_def_sc_output_fault, 1};
        fault_output.set_nan();
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("fault", 1, 1, nullptr, fault_output);
        model.calculate(options, sc_output);
        std::array<std::array<double, 3>, 1> i_f{};
        fault_output.get_value(PGM_def_sc_output_fault_i_f, i_f.data(), -1);
        CHECK(i_f[0][0] > 0.0);
        CHECK(i_f[0][1] == doctest::Approx(i_f[0][0]));
        CHECK(i_f[0][2] == doctest::Approx(i_f[0][0]));
    }

    SUBCASE("Short circuit without faults") {
        auto const owning_no_fault_dataset = load_dataset(reduced_phase_model_json(false));
        Model no_fault_model{50.0, owning_no_fault_dataset.dataset};
        options.set_calculation_type(PGM_short_circuit);
        options.set_calculation_method(PGM_iec60909);
        Buffer sc_node_output{PGM_def_sc_output_node, 2};
        sc_node_output.set_nan();
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("node", 2, 2, nullptr, sc_node_output);
        CHECK_NOTHROW(no_fault_model.calculate(options, sc_output));
        std::array<std::array<double, 3>, 2> u{};
        sc_node_output.get_value(PGM_def_sc_output_node_u, u.data(), -1);
        CHECK(u[1][0] > 0.0);
        CHECK(u[1][1] == 0.0);
    }
}

namespace {
// node 1 (abc, 22 kV, source) -- single_phase_transformer ab -> a (22 kV / 12.7 kV, earth electrode 5 ohm) -- node 2
// (SWER, u_rated 22 kV) with a load on phase a and a single-phase-to-earth fault
auto single_phase_transformer_model_json() {
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 22000},
      {"id": 2, "u_rated": 22000}
    ],
    "source": [
      {"id": 3, "node": 1, "status": 1, "u_ref": 1.0, "sk": 1e9}
    ],
    "single_phase_transformer": [
      {"id": 4, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "u1": 22000, "u2": 12700,
       "sn": 200000, "uk": 0.04, "pk": 2000, "i0": 0, "p0": 0, "from_connection": 4, "to_connection": 1,
       "r_grounding_to": 5}
    ],
    "asym_load": [
      {"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": [50000, 0, 0], "q_specified": [0, 0, 0]}
    ],
    "fault": [
      {"id": 6, "status": 1, "fault_type": 1, "fault_phase": 1, "fault_object": 2}
    ]
  }
})json";
}
} // namespace

TEST_CASE("API Model - single phase transformer") {
    using namespace std::string_literals;

    auto const owning_input_dataset = load_dataset(single_phase_transformer_model_json());
    Model model{50.0, owning_input_dataset.dataset};

    Options options{};
    options.set_experimental_features(PGM_experimental_features_enabled);
    options.set_symmetric(0);
    options.set_calculation_method(PGM_newton_raphson);

    Buffer node_output{PGM_def_asym_output_node, 2};
    node_output.set_nan();
    Buffer transformer_output{PGM_def_asym_output_single_phase_transformer, 1};
    transformer_output.set_nan();
    DatasetMutable asym_output{"asym_output", false, 1};
    asym_output.add_buffer("node", 2, 2, nullptr, node_output);
    asym_output.add_buffer("single_phase_transformer", 1, 1, nullptr, transformer_output);

    SUBCASE("Asymmetric power flow") {
        model.calculate(options, asym_output);
        std::array<std::array<double, 3>, 2> u{};
        std::array<std::array<double, 3>, 2> u_angle{};
        node_output.get_value(PGM_def_asym_output_node_u, u.data(), -1);
        node_output.get_value(PGM_def_asym_output_node_u_angle, u_angle.data(), -1);
        CHECK(u[1][0] > 0.0);
        CHECK(u[1][1] == 0.0);
        CHECK(u[1][2] == 0.0);
        // the SWER voltage follows U_a - U_b, which leads U_a by 30 degrees
        CHECK(std::abs(u_angle[1][0] - std::numbers::pi / 6.0) < 0.05);

        std::array<std::array<double, 3>, 1> p_to{};
        std::array<double, 1> loading{};
        transformer_output.get_value(PGM_def_asym_output_single_phase_transformer_p_to, p_to.data(), -1);
        transformer_output.get_value(PGM_def_asym_output_single_phase_transformer_loading, loading.data(), -1);
        CHECK(p_to[0][0] < 0.0);
        CHECK(p_to[0][1] == 0.0);
        CHECK(p_to[0][2] == 0.0);
        CHECK(loading[0] > 0.0);
    }

    SUBCASE("Experimental features are required") {
        options.set_experimental_features(PGM_experimental_features_disabled);
        check_throws_with([&model, &asym_output](auto const& opt) { model.calculate(opt, asym_output); },
                          PGM_regular_error, "is an experimental feature"s, options);
    }

    SUBCASE("Symmetric power flow is not supported") {
        Buffer sym_node_output{PGM_def_sym_output_node, 2};
        DatasetMutable sym_output{"sym_output", false, 1};
        sym_output.add_buffer("node", 2, 2, nullptr, sym_node_output);
        options.set_symmetric(1);
        check_throws_with([&model, &sym_output](auto const& opt) { model.calculate(opt, sym_output); },
                          PGM_regular_error,
                          "Symmetric power flow is not supported in a network with one- or two-phase nodes"s, options);
    }

    SUBCASE("Single-phase-to-earth fault at the SWER node") {
        options.set_calculation_type(PGM_short_circuit);
        options.set_calculation_method(PGM_iec60909);
        Buffer fault_output{PGM_def_sc_output_fault, 1};
        fault_output.set_nan();
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("fault", 1, 1, nullptr, fault_output);
        model.calculate(options, sc_output);
        std::array<std::array<double, 3>, 1> i_f{};
        fault_output.get_value(PGM_def_sc_output_fault_i_f, i_f.data(), -1);
        CHECK(i_f[0][0] > 0.0);
        CHECK(i_f[0][1] == 0.0);
        CHECK(i_f[0][2] == 0.0);
    }
}

namespace {
// node 1 (abc, 22 kV, source) -- single_phase_transformer a -> ab_center_tapped (12.7 kV / 240-0-240 V, solidly
// earthed centre tap, no magnetising) -- node 2 (split phase, u_rated 240 sqrt3); extra_components holds further
// components of node 2, each preceded by a comma
std::string split_phase_model_json(std::string_view extra_components) {
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 22000},
      {"id": 2, "u_rated": 415.69219381653056}
    ],
    "source": [
      {"id": 3, "node": 1, "status": 1, "u_ref": 1.0, "sk": 1e9}
    ],
    "single_phase_transformer": [
      {"id": 4, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "u1": 12701.705922171768, "u2": 480,
       "sn": 25000, "uk": 0.04, "pk": 400, "i0": 0, "p0": 0, "from_connection": 1, "to_connection": 7}
    ])json"s +
           std::string{extra_components} + R"json(
  }
})json"s;
}

constexpr std::string_view split_phase_unbalanced_loads = R"json(,
    "asym_load": [
      {"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": [6000, 2000, 0], "q_specified": [2000, 500, 0]},
      {"id": 6, "node": 2, "status": 1, "type": 2, "p_specified": [1000, 0, 0], "q_specified": [0, 0, 0]}
    ])json";

// FaultType and FaultPhase values
constexpr int fault_type_three_phase = 0;
constexpr int fault_type_single_phase_to_ground = 1;
constexpr int fault_type_two_phase = 2;
constexpr int fault_type_two_phase_to_ground = 3;
constexpr int fault_phase_abc = 0;
constexpr int fault_phase_b = 2;
constexpr int fault_phase_ab = 4;

// a load without power: the solution stays that of the unloaded network, but the iterative methods are used (with
// only constant-impedance loads, or none, the calculation switches to the linear method)
constexpr std::string_view split_phase_zero_load = R"json(,
    "asym_load": [
      {"id": 5, "node": 2, "status": 1, "type": 0, "p_specified": [0, 0, 0], "q_specified": [0, 0, 0]}
    ])json";

std::string split_phase_fault(int fault_type, int fault_phase) {
    return R"json(,
    "fault": [
      {"id": 7, "status": 1, "fault_type": )json"s +
           std::to_string(fault_type) + R"json(, "fault_phase": )json"s + std::to_string(fault_phase) +
           R"json(, "fault_object": 2}
    ])json"s;
}
} // namespace

TEST_CASE("API Model - centre-tapped single phase transformer") {
    using namespace std::string_literals;

    Options options{};
    options.set_experimental_features(PGM_experimental_features_enabled);
    options.set_symmetric(0);

    Buffer node_output{PGM_def_asym_output_node, 2};
    node_output.set_nan();
    DatasetMutable asym_output{"asym_output", false, 1};
    asym_output.add_buffer("node", 2, 2, nullptr, node_output);

    SUBCASE("The start values of the iterative methods carry the split-phase pattern") {
        // without load and losses the solution is U_1 at the source voltage and U_2 = (1, -1, 0) p.u., which is the
        // flat and the average-source start when slot b starts opposite to slot a; a start at the 120 degree pattern is
        // 60 degrees off at slot b and does not converge in one iteration
        auto const owning_input_dataset = load_dataset(split_phase_model_json(split_phase_zero_load));
        Model model{50.0, owning_input_dataset.dataset};
        options.set_max_iter(1);
        auto const check_one_iteration = [&model, &options, &node_output, &asym_output](Idx method,
                                                                                        Idx initialization) {
            CAPTURE(method);
            CAPTURE(initialization);
            options.set_calculation_method(method);
            options.set_calculation_initialization(initialization);
            node_output.set_nan();
            REQUIRE_NOTHROW(model.calculate(options, asym_output));
            std::array<std::array<double, 3>, 2> u_pu{};
            std::array<std::array<double, 3>, 2> u_angle{};
            node_output.get_value(PGM_def_asym_output_node_u_pu, u_pu.data(), -1);
            node_output.get_value(PGM_def_asym_output_node_u_angle, u_angle.data(), -1);
            CHECK(u_pu[1][0] == doctest::Approx(1.0));
            CHECK(u_pu[1][1] == doctest::Approx(1.0));
            CHECK(u_pu[1][2] == 0.0);
            CHECK(u_angle[1][0] == doctest::Approx(0.0).epsilon(1e-9));
            CHECK(std::abs(u_angle[1][1]) == doctest::Approx(std::numbers::pi));
        };
        SUBCASE("Newton-Raphson with a flat start") {
            check_one_iteration(PGM_newton_raphson, PGM_calculation_initialization_flat);
        }
        SUBCASE("Newton-Raphson with an average-source start") {
            check_one_iteration(PGM_newton_raphson, PGM_calculation_initialization_average_source);
        }
        SUBCASE("Iterative current") {
            check_one_iteration(PGM_iterative_current, PGM_calculation_initialization_default);
        }
    }

    SUBCASE("Unbalanced loads") {
        auto const owning_input_dataset = load_dataset(split_phase_model_json(split_phase_unbalanced_loads));
        Model model{50.0, owning_input_dataset.dataset};
        Buffer transformer_output{PGM_def_asym_output_single_phase_transformer, 1};
        transformer_output.set_nan();
        asym_output.add_buffer("single_phase_transformer", 1, 1, nullptr, transformer_output);
        for (auto const [method, initialization] :
             {std::pair{PGM_newton_raphson, PGM_calculation_initialization_default},
              std::pair{PGM_newton_raphson, PGM_calculation_initialization_flat},
              std::pair{PGM_newton_raphson, PGM_calculation_initialization_average_source},
              std::pair{PGM_iterative_current, PGM_calculation_initialization_default}}) {
            CAPTURE(static_cast<int>(method));
            CAPTURE(static_cast<int>(initialization));
            options.set_calculation_method(method);
            options.set_calculation_initialization(initialization);
            node_output.set_nan();
            model.calculate(options, asym_output);
            std::array<std::array<double, 3>, 2> u{};
            std::array<std::array<double, 3>, 2> u_angle{};
            node_output.get_value(PGM_def_asym_output_node_u, u.data(), -1);
            node_output.get_value(PGM_def_asym_output_node_u_angle, u_angle.data(), -1);
            CHECK(u[1][0] > 0.0);
            CHECK(u[1][1] > 0.0);
            CHECK(u[1][0] < u[1][1]); // leg a is loaded more
            CHECK(u[1][2] == 0.0);
            // the legs stay close to opposite
            double const leg_angle = std::remainder(u_angle[1][1] - u_angle[1][0], 2.0 * std::numbers::pi);
            CHECK(std::abs(std::abs(leg_angle) - std::numbers::pi) < 0.05);

            std::array<double, 1> loading{};
            transformer_output.get_value(PGM_def_asym_output_single_phase_transformer_loading, loading.data(), -1);
            CHECK(loading[0] > 0.0);
        }
    }

    SUBCASE("Short circuits at the split-phase node") {
        options.set_calculation_type(PGM_short_circuit);
        options.set_calculation_method(PGM_iec60909);
        Buffer fault_output{PGM_def_sc_output_fault, 1};
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("fault", 1, 1, nullptr, fault_output);
        auto const fault_current = [&options, &fault_output, &sc_output](int fault_type, int fault_phase) {
            auto const owning_input_dataset =
                load_dataset(split_phase_model_json(split_phase_fault(fault_type, fault_phase)));
            Model model{50.0, owning_input_dataset.dataset};
            fault_output.set_nan();
            model.calculate(options, sc_output);
            std::array<std::array<double, 3>, 1> i_f{};
            fault_output.get_value(PGM_def_sc_output_fault_i_f, i_f.data(), -1);
            return i_f[0];
        };

        SUBCASE("Two-phase ab") {
            auto const i_f = fault_current(fault_type_two_phase, fault_phase_ab);
            CHECK(i_f[0] > 0.0);
            CHECK(i_f[1] == doctest::Approx(i_f[0]));
            CHECK(i_f[2] == 0.0);
        }
        SUBCASE("Single-phase-to-ground b") {
            auto const i_f = fault_current(fault_type_single_phase_to_ground, fault_phase_b);
            CHECK(i_f[0] == 0.0);
            CHECK(i_f[1] > 0.0);
            CHECK(i_f[2] == 0.0);
        }
        SUBCASE("Three-phase is not possible") {
            check_throws_with([&fault_current]() { fault_current(fault_type_three_phase, fault_phase_abc); },
                              PGM_regular_error,
                              "The short circuit phases (abc) of fault 7 are not all present at its node 2 (node "
                              "phases: ab)"s);
        }
    }
}

namespace {
constexpr int fault_phase_a = 1;
constexpr int fault_phase_bc = 6;
constexpr int fault_phase_not_given = -128; // nan: the default of the fault type

// node 1 (abc, 10 kV, source) -- asym_line 11 (a) -- node 2; node 1 -- asym_line 12 (b) -- node 3; node 4 is not
// connected to anything (isolated); faults holds the fault components, preceded by a comma
std::string fault_phase_model_json(std::string_view faults) {
    using namespace std::string_literals;
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 10000},
      {"id": 2, "u_rated": 10000},
      {"id": 3, "u_rated": 10000},
      {"id": 4, "u_rated": 10000}
    ],
    "source": [
      {"id": 5, "node": 1, "status": 1, "u_ref": 1.0, "sk": 1e8}
    ],
    "asym_line": [
      {"id": 11, "from_node": 1, "to_node": 2, "from_status": 1, "to_status": 1, "phases": 1,
       "r_aa": 0.5, "x_aa": 0.4, "c_aa": 1e-8},
      {"id": 12, "from_node": 1, "to_node": 3, "from_status": 1, "to_status": 1, "phases": 2,
       "r_bb": 0.5, "x_bb": 0.4, "c_bb": 1e-8}
    ])json"s +
           std::string{faults} + R"json(
  }
})json"s;
}

// a fault component entry; fault_phase_not_given leaves fault_phase out; r_f in ohm (0: bolted)
std::string fault_json(int id, int status, int fault_type, int fault_phase, int fault_object, int r_f = 0) {
    using namespace std::string_literals;
    std::string const phase =
        fault_phase == fault_phase_not_given ? ""s : R"json(, "fault_phase": )json"s + std::to_string(fault_phase);
    return R"json({"id": )json"s + std::to_string(id) + R"json(, "status": )json"s + std::to_string(status) +
           R"json(, "fault_type": )json"s + std::to_string(fault_type) + phase + R"json(, "fault_object": )json"s +
           std::to_string(fault_object) + R"json(, "r_f": )json"s + std::to_string(r_f) + "}"s;
}

std::string faults_json(std::initializer_list<std::string> faults) {
    using namespace std::string_literals;
    std::string result = R"json(,
    "fault": [
      )json"s;
    bool first = true;
    for (auto const& fault : faults) {
        if (!first) {
            result += ",\n      "s;
        }
        result += fault;
        first = false;
    }
    return result + "\n    ]"s;
}
} // namespace

TEST_CASE("API Model - short circuit fault phases") {
    using namespace std::string_literals;

    Options options{};
    options.set_experimental_features(PGM_experimental_features_enabled);
    options.set_calculation_type(PGM_short_circuit);
    options.set_calculation_method(PGM_iec60909);

    // the fault currents of all faults, in input order
    auto const fault_currents = [&options](std::initializer_list<std::string> faults) {
        auto const owning_input_dataset = load_dataset(fault_phase_model_json(faults_json(faults)));
        Model model{50.0, owning_input_dataset.dataset};
        auto const n_faults = static_cast<Idx>(faults.size());
        Buffer fault_output{PGM_def_sc_output_fault, n_faults};
        fault_output.set_nan();
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("fault", n_faults, n_faults, nullptr, fault_output);
        model.calculate(options, sc_output);
        std::vector<std::array<double, 3>> i_f(faults.size());
        fault_output.get_value(PGM_def_sc_output_fault_i_f, i_f.data(), -1);
        return i_f;
    };
    auto const single_fault = [](int fault_type, int fault_phase, int node) {
        return fault_json(31, 1, fault_type, fault_phase, node);
    };

    SUBCASE("Faults on present phases") {
        SUBCASE("Single-phase-to-ground a at the a node") {
            auto const i_f = fault_currents({single_fault(fault_type_single_phase_to_ground, fault_phase_a, 2)});
            CHECK(i_f[0][0] > 0.0);
            CHECK(i_f[0][1] == 0.0);
            CHECK(i_f[0][2] == 0.0);
        }
        SUBCASE("Single-phase-to-ground with the default phase (a) at the a node") {
            auto const i_f =
                fault_currents({single_fault(fault_type_single_phase_to_ground, fault_phase_not_given, 2)});
            CHECK(i_f[0][0] > 0.0);
            CHECK(i_f[0][1] == 0.0);
            CHECK(i_f[0][2] == 0.0);
        }
        SUBCASE("Single-phase-to-ground b at the b node") {
            auto const i_f = fault_currents({single_fault(fault_type_single_phase_to_ground, fault_phase_b, 3)});
            CHECK(i_f[0][0] == 0.0);
            CHECK(i_f[0][1] > 0.0);
            CHECK(i_f[0][2] == 0.0);
        }
        SUBCASE("Two simultaneous faults at the abc node and the a node") {
            // the fault at node 1 has a fault resistance; a bolted one would leave no voltage for the fault at node 2
            auto const i_f =
                fault_currents({fault_json(31, 1, fault_type_single_phase_to_ground, fault_phase_a, 1, 100),
                                fault_json(32, 1, fault_type_single_phase_to_ground, fault_phase_a, 2)});
            CHECK(i_f[0][0] > 0.0);
            CHECK(i_f[1][0] > 0.0);
        }
    }

    SUBCASE("Faults on absent phases name the fault and its node") {
        SUBCASE("Single-phase-to-ground b at the a node") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_single_phase_to_ground, fault_phase_b, 2)});
                },
                PGM_regular_error,
                "The short circuit phases (b) of fault 31 are not all present at its node 2 (node phases: a); this "
                "short circuit type cannot be applied there\n"s);
        }
        SUBCASE("Three-phase at the a node") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_three_phase, fault_phase_abc, 2)});
                },
                PGM_regular_error, "The short circuit phases (abc) of fault 31 are not all present at its node 2"s);
        }
        SUBCASE("Two-phase bc at the a node") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_two_phase, fault_phase_bc, 2)});
                },
                PGM_regular_error, "The short circuit phases (bc) of fault 31 are not all present at its node 2"s);
        }
    }

    SUBCASE("Default phases are not adapted to the node") {
        SUBCASE("Single-phase-to-ground at the b node") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_single_phase_to_ground, fault_phase_not_given, 3)});
                },
                PGM_regular_error,
                "The short circuit phases (a) of fault 31 are not all present at its node 3 (node phases: b); this "
                "short circuit type cannot be applied there; the phases are the default of the fault type, so set "
                "fault_phase explicitly\n"s);
        }
        SUBCASE("Two-phase at the b node") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_two_phase, fault_phase_not_given, 3)});
                },
                PGM_regular_error,
                "The short circuit phases (bc) of fault 31 are not all present at its node 3 (node phases: b); this "
                "short circuit type cannot be applied there\n"s);
        }
    }

    SUBCASE("The hint is only given when an explicit phase can fit the node") {
        // the message ends directly after "there", so there is no hint
        SUBCASE("Default three-phase at the a node: no explicit phase can fit") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_three_phase, fault_phase_not_given, 2)});
                },
                PGM_regular_error,
                "The short circuit phases (abc) of fault 31 are not all present at its node 2 (node phases: a); this "
                "short circuit type cannot be applied there\n"s);
        }
        SUBCASE("Default two-phase at the a node: the node has fewer than two phases") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_two_phase, fault_phase_not_given, 2)});
                },
                PGM_regular_error,
                "The short circuit phases (bc) of fault 31 are not all present at its node 2 (node phases: a); this "
                "short circuit type cannot be applied there\n"s);
        }
        SUBCASE("Default two-phase-to-ground at the b node: the node has fewer than two phases") {
            check_throws_with(
                [&fault_currents, &single_fault]() {
                    fault_currents({single_fault(fault_type_two_phase_to_ground, fault_phase_not_given, 3)});
                },
                PGM_regular_error,
                "The short circuit phases (bc) of fault 31 are not all present at its node 3 (node phases: b); this "
                "short circuit type cannot be applied there\n"s);
        }
    }

    SUBCASE("Skipped faults") {
        SUBCASE("A fault with status off") {
            auto const i_f = fault_currents({fault_json(31, 1, fault_type_single_phase_to_ground, fault_phase_a, 1),
                                             fault_json(32, 0, fault_type_single_phase_to_ground, fault_phase_b, 2)});
            CHECK(i_f[0][0] > 0.0);
            CHECK(i_f[1][0] == 0.0);
            CHECK(i_f[1][1] == 0.0);
            CHECK(i_f[1][2] == 0.0);
        }
        SUBCASE("A fault at an isolated node") {
            auto const i_f = fault_currents({single_fault(fault_type_three_phase, fault_phase_abc, 4)});
            CHECK(i_f[0][0] == 0.0);
            CHECK(i_f[0][1] == 0.0);
            CHECK(i_f[0][2] == 0.0);
        }
    }
}

namespace {
// LoadGenType values
constexpr int load_gen_type_const_pq = 0;
constexpr int load_gen_type_const_y = 1;

// node 1 (abc, 400 V, source with the default sk), optionally node 4 (400 V), with the given components, each preceded
// by a comma
std::string delta_load_model_json(std::string_view components, bool with_node_4 = false) {
    using namespace std::string_literals;
    return R"json({
  "version": "1.0",
  "type": "input",
  "is_batch": false,
  "attributes": {},
  "data": {
    "node": [
      {"id": 1, "u_rated": 400})json"s +
           (with_node_4 ? R"json(,
      {"id": 4, "u_rated": 400})json"s
                        : ""s) +
           R"json(
    ],
    "source": [
      {"id": 2, "node": 1, "status": 1, "u_ref": 1.0}
    ])json"s +
           std::string{components} + R"json(
  }
})json"s;
}

// an asym_load with the given connection (1 = delta), type and specified active power at the given node, no reactive
// power
std::string delta_asym_load(int node, int type, std::array<double, 3> const& p_specified) {
    using namespace std::string_literals;
    return R"json(,
    "asym_load": [
      {"id": 3, "node": )json"s +
           std::to_string(node) + R"json(, "status": 1, "type": )json"s + std::to_string(type) +
           R"json(, "p_specified": [)json"s + std::to_string(p_specified[0]) + ", "s + std::to_string(p_specified[1]) +
           ", "s + std::to_string(p_specified[2]) + R"json(], "q_specified": [0, 0, 0], "connection": 1}
    ])json"s;
}

struct LoadGenPhaseResult {
    std::array<double, 3> p{};
    std::array<double, 3> q{};
    std::array<double, 3> i{};
};

LoadGenPhaseResult get_asym_load_result(Buffer const& load_output) {
    LoadGenPhaseResult result;
    load_output.get_value(PGM_def_asym_output_asym_load_p, result.p.data(), -1);
    load_output.get_value(PGM_def_asym_output_asym_load_q, result.q.data(), -1);
    load_output.get_value(PGM_def_asym_output_asym_load_i, result.i.data(), -1);
    return result;
}
} // namespace

TEST_CASE("API Model - delta load") {
    using namespace std::string_literals;

    // 40 kW between a and b at a 400 V node: 100 A in both conductors; at balanced voltages each phase takes half of
    // the power, with -P/(2 sqrt3) and +P/(2 sqrt3) reactive power on a and b
    constexpr double p_delta = 40e3;
    auto const const_pq_load = delta_load_model_json(delta_asym_load(1, load_gen_type_const_pq, {p_delta, 0.0, 0.0}));

    Options options{};
    options.set_experimental_features(PGM_experimental_features_enabled);
    options.set_symmetric(0);

    Buffer load_output{PGM_def_asym_output_asym_load, 1};
    load_output.set_nan();
    DatasetMutable asym_output{"asym_output", false, 1};
    asym_output.add_buffer("asym_load", 1, 1, nullptr, load_output);

    SUBCASE("Experimental features are required") {
        auto const owning_input_dataset = load_dataset(const_pq_load);
        Model model{50.0, owning_input_dataset.dataset};
        options.set_experimental_features(PGM_experimental_features_disabled);
        check_throws_with([&model, &asym_output](auto const& opt) { model.calculate(opt, asym_output); },
                          PGM_regular_error, "A delta-connected load or generator is an experimental feature"s,
                          options);
    }

    SUBCASE("Asymmetric power flow") {
        for (auto const [method, type] :
             {std::pair{PGM_newton_raphson, load_gen_type_const_pq},
              std::pair{PGM_iterative_current, load_gen_type_const_pq}, std::pair{PGM_linear, load_gen_type_const_y}}) {
            CAPTURE(static_cast<int>(method));
            CAPTURE(type);
            auto const owning_input_dataset =
                load_dataset(delta_load_model_json(delta_asym_load(1, type, {p_delta, 0.0, 0.0})));
            Model model{50.0, owning_input_dataset.dataset};
            options.set_calculation_method(method);
            load_output.set_nan();
            model.calculate(options, asym_output);
            auto const [p, q, i] = get_asym_load_result(load_output);

            double const q_ab = p_delta / (2.0 * std::numbers::sqrt3);
            CHECK(p[0] == doctest::Approx(p_delta / 2.0).epsilon(1e-4));
            CHECK(p[1] == doctest::Approx(p_delta / 2.0).epsilon(1e-4));
            CHECK(q[0] == doctest::Approx(-q_ab).epsilon(1e-4));
            CHECK(q[1] == doctest::Approx(q_ab).epsilon(1e-4));
            CHECK(i[0] == doctest::Approx(100.0).epsilon(1e-4));
            CHECK(i[1] == doctest::Approx(100.0).epsilon(1e-4));
            // the c conductor carries nothing
            CHECK(p[2] == 0.0);
            CHECK(q[2] == 0.0);
            CHECK(i[2] == 0.0);
            if (type == load_gen_type_const_pq) {
                // the branch power is exactly the specified power
                CHECK(p[0] + p[1] == doctest::Approx(p_delta).epsilon(1e-9));
                CHECK(q[0] + q[1] == doctest::Approx(0.0).epsilon(1e-9).scale(p_delta));
            }
        }
    }

    SUBCASE("Symmetric power flow is not supported") {
        auto const owning_input_dataset = load_dataset(const_pq_load);
        Model model{50.0, owning_input_dataset.dataset};
        Buffer sym_load_output{PGM_def_sym_output_asym_load, 1};
        DatasetMutable sym_output{"sym_output", false, 1};
        sym_output.add_buffer("asym_load", 1, 1, nullptr, sym_load_output);
        options.set_symmetric(1);
        check_throws_with(
            [&model, &sym_output](auto const& opt) { model.calculate(opt, sym_output); }, PGM_regular_error,
            "Symmetric power flow is not supported in a network with delta-connected loads or generators"s, options);
    }

    SUBCASE("State estimation is not supported") {
        auto const owning_input_dataset = load_dataset(const_pq_load);
        Model model{50.0, owning_input_dataset.dataset};
        Buffer sym_load_output{PGM_def_sym_output_asym_load, 1};
        DatasetMutable sym_output{"sym_output", false, 1};
        sym_output.add_buffer("asym_load", 1, 1, nullptr, sym_load_output);
        options.set_calculation_type(PGM_state_estimation);
        for (Idx const symmetric : {0, 1}) {
            options.set_symmetric(symmetric);
            auto& output = symmetric == 0 ? asym_output : sym_output;
            check_throws_with(
                [&model, &output](auto const& opt) { model.calculate(opt, output); }, PGM_regular_error,
                "State estimation is not supported in a network with delta-connected loads or generators"s, options);
        }
    }

    SUBCASE("A three-phase fault is calculated symmetrically and the load is ignored") {
        auto const owning_input_dataset = load_dataset(
            delta_load_model_json(delta_asym_load(1, load_gen_type_const_pq, {p_delta, 0.0, 0.0}) + R"json(,
    "fault": [
      {"id": 9, "status": 1, "fault_type": 0, "fault_phase": 0, "fault_object": 1}
    ])json"s));
        Model model{50.0, owning_input_dataset.dataset};
        options.set_calculation_type(PGM_short_circuit);
        options.set_calculation_method(PGM_iec60909);
        Buffer fault_output{PGM_def_sc_output_fault, 1};
        fault_output.set_nan();
        Buffer sc_load_output{PGM_def_sc_output_asym_load, 1};
        sc_load_output.set_nan();
        DatasetMutable sc_output{"sc_output", false, 1};
        sc_output.add_buffer("fault", 1, 1, nullptr, fault_output);
        sc_output.add_buffer("asym_load", 1, 1, nullptr, sc_load_output);
        REQUIRE_NOTHROW(model.calculate(options, sc_output));
        std::array<double, 3> i_f{};
        std::array<double, 3> i_load{};
        fault_output.get_value(PGM_def_sc_output_fault_i_f, i_f.data(), -1);
        sc_load_output.get_value(PGM_def_sc_output_asym_load_i, i_load.data(), -1);
        CHECK(i_f[0] > 0.0);
        CHECK(i_f[1] == doctest::Approx(i_f[0]));
        CHECK(i_f[2] == doctest::Approx(i_f[0]));
        CHECK(i_load == std::array{0.0, 0.0, 0.0});
    }

    SUBCASE("A sym_load cannot be delta") {
        auto const owning_input_dataset = load_dataset(delta_load_model_json(R"json(,
    "sym_load": [
      {"id": 6, "node": 1, "status": 1, "type": 0, "p_specified": 1000, "q_specified": 0, "connection": 1}
    ])json"));
        check_throws_with([&owning_input_dataset]() { Model{50.0, owning_input_dataset.dataset}; }, PGM_regular_error,
                          "Invalid phase connection for sym_load 6 with phases (delta): a delta connection is only "
                          "supported for asym_load and asym_gen"s);
    }

    SUBCASE("A delta generator cannot be voltage regulated") {
        auto const owning_input_dataset = load_dataset(delta_load_model_json(R"json(,
    "asym_gen": [
      {"id": 6, "node": 1, "status": 1, "type": 0, "p_specified": [1000, 0, 0], "q_specified": [0, 0, 0],
       "connection": 1}
    ],
    "voltage_regulator": [
      {"id": 7, "regulated_object": 6, "status": 1, "u_ref": 1.0}
    ])json"));
        Model model{50.0, owning_input_dataset.dataset};
        options.set_calculation_method(PGM_newton_raphson);
        Buffer gen_output{PGM_def_asym_output_asym_gen, 1};
        DatasetMutable gen_dataset{"asym_output", false, 1};
        gen_dataset.add_buffer("asym_gen", 1, 1, nullptr, gen_output);
        check_throws_with([&model, &gen_dataset](auto const& opt) { model.calculate(opt, gen_dataset); },
                          PGM_regular_error, "Unsupported load_gen type for voltage regulators 7."s, options);
    }

    SUBCASE("A 480 V delta load at a split-phase node") {
        // the network is antisymmetric (U_b = -U_a), so both legs take half of the power; the nominal branch voltage
        // of the delta entry is 2 u_rated / sqrt3 = 480 V
        constexpr double p_split = 4e3;
        Buffer node_output{PGM_def_asym_output_node, 2};
        asym_output.add_buffer("node", 2, 2, nullptr, node_output);
        for (auto const [method, type] : {std::pair{PGM_newton_raphson, load_gen_type_const_pq},
                                          std::pair{PGM_newton_raphson, load_gen_type_const_y},
                                          std::pair{PGM_iterative_current, load_gen_type_const_y}}) {
            CAPTURE(static_cast<int>(method));
            CAPTURE(type);
            auto const load = R"json(,
    "asym_load": [
      {"id": 5, "node": 2, "status": 1, "type": )json"s +
                              std::to_string(type) + R"json(, "p_specified": [4000, 0, 0], "q_specified": [0, 0, 0],
       "connection": 1}
    ])json"s;
            auto const owning_input_dataset = load_dataset(split_phase_model_json(load));
            Model model{50.0, owning_input_dataset.dataset};
            options.set_calculation_method(method);
            load_output.set_nan();
            node_output.set_nan();
            model.calculate(options, asym_output);
            auto const [p, q, i] = get_asym_load_result(load_output);
            std::array<std::array<double, 3>, 2> u{};
            std::array<std::array<double, 3>, 2> u_angle{};
            node_output.get_value(PGM_def_asym_output_node_u, u.data(), -1);
            node_output.get_value(PGM_def_asym_output_node_u_angle, u_angle.data(), -1);

            CHECK(p[0] == doctest::Approx(p[1]).epsilon(1e-9));
            CHECK(q[0] == doctest::Approx(q[1]).epsilon(1e-9).scale(p_split));
            CHECK(i[0] == doctest::Approx(i[1]).epsilon(1e-9));
            CHECK(p[2] == 0.0);
            CHECK(i[2] == 0.0);
            double const u_ab = std::abs(std::polar(u[1][0], u_angle[1][0]) - std::polar(u[1][1], u_angle[1][1]));
            CHECK(u_ab == doctest::Approx(480.0).epsilon(0.05));
            double const expected = type == load_gen_type_const_y ? p_split * (u_ab / 480.0) * (u_ab / 480.0) : p_split;
            CHECK(p[0] + p[1] == doctest::Approx(expected).epsilon(1e-9));
        }
    }

    SUBCASE("A delta entry needs both of its phases") {
        // node 4 behind an ab line
        auto const ab_line_load = [](std::array<double, 3> const& p_specified) {
            return delta_load_model_json(R"json(,
    "asym_line": [
      {"id": 5, "from_node": 1, "to_node": 4, "from_status": 1, "to_status": 1, "phases": 4,
       "r_aa": 0.05, "x_aa": 0.04, "r_ba": 0.01, "x_ba": 0.02, "r_bb": 0.05, "x_bb": 0.04,
       "c_aa": 1e-9, "c_ba": 0, "c_bb": 1e-9}
    ])json"s + delta_asym_load(4, load_gen_type_const_pq, p_specified),
                                         true);
        };
        options.set_calculation_method(PGM_newton_raphson);
        SUBCASE("An ab entry") {
            auto const owning_input_dataset = load_dataset(ab_line_load({1000.0, 0.0, 0.0}));
            Model model{50.0, owning_input_dataset.dataset};
            model.calculate(options, asym_output);
            auto const [p, q, i] = get_asym_load_result(load_output);
            CHECK(p[0] + p[1] == doctest::Approx(1000.0).epsilon(1e-9));
            CHECK(i[0] == doctest::Approx(i[1]).epsilon(1e-9));
            CHECK(p[2] == 0.0);
            CHECK(i[2] == 0.0);
        }
        SUBCASE("A bc entry") {
            auto const owning_input_dataset = load_dataset(ab_line_load({0.0, 1000.0, 0.0}));
            Model model{50.0, owning_input_dataset.dataset};
            check_throws_with([&model, &asym_output](auto const& opt) { model.calculate(opt, asym_output); },
                              PGM_regular_error,
                              "Load or generator 3 has a non-zero specified power on phases (c) that are not present "
                              "at its node (node phases: ab). A delta (phase-to-phase) entry needs both of its phases"s,
                              options);
        }
    }
}

} // namespace power_grid_model_cpp
