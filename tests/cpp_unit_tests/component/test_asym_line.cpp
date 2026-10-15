// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/component/asym_line.hpp>

#include <power_grid_model/auxiliary/input.hpp>
#include <power_grid_model/auxiliary/output.hpp>
#include <power_grid_model/auxiliary/update.hpp>
#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/exception.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>
#include <power_grid_model/component/base.hpp>
#include <power_grid_model/component/branch.hpp>
#include <power_grid_model/component/line_utils.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace power_grid_model {

using namespace std::complex_literals;

namespace {
// this function exceeds the reasonable threshold of 800 statements
// NOLINTNEXTLINE(readability-function-size,google-readability-function-size,hicpp-function-size)
void execute_subcases(const AsymLineInput& input, const ComplexTensor<asymmetric_t>& y_series,
                      const ComplexTensor<asymmetric_t>& c_matrix, const double base_i, const double base_y,
                      const double system_frequency, const double voltage_lvl) {
    CAPTURE(input.id);
    CAPTURE(y_series);
    AsymLine asym_line{input, system_frequency, voltage_lvl, voltage_lvl};
    Branch& branch = asym_line;

    ComplexTensor<asymmetric_t> const y_shunt = 1 / base_y * (2 * pi * system_frequency * c_matrix * 1.0i);

    DoubleComplex const y1_series =
        (y_series(0, 0) + y_series(1, 1) + y_series(2, 2)) / 3.0 -
        (y_series(0, 1) + y_series(0, 2) + y_series(1, 0) + y_series(1, 2) + y_series(2, 0) + y_series(2, 1)) / 6.0;
    DoubleComplex const y1_shunt =
        (y_shunt(0, 0) + y_shunt(1, 1) + y_shunt(2, 2)) / 3.0 -
        (y_shunt(0, 1) + y_shunt(0, 2) + y_shunt(1, 0) + y_shunt(1, 2) + y_shunt(2, 0) + y_shunt(2, 1)) / 6.0;

    // symmetric
    DoubleComplex const yff1 = y1_series + 0.5 * y1_shunt;
    DoubleComplex const yft1 = -y1_series;
    DoubleComplex const ys1 = 0.5 * y1_shunt + 1.0 / (1.0 / y1_series + 2.0 / y1_shunt);

    // asymmetric
    ComplexTensor<asymmetric_t> const ytt = y_series + 0.5 * y_shunt;
    ComplexTensor<asymmetric_t> const branch_shunt = 0.5 * y_shunt + inv(inv(y_series) + 2.0 * inv(y_shunt));

    // Short circuit results
    DoubleComplex const if_sc{1.0, 1.0};
    DoubleComplex const it_sc{2.0, 2.0 * sqrt3};
    ComplexValue<asymmetric_t> const if_sc_asym{1.0 + 1.0i};
    ComplexValue<asymmetric_t> const it_sc_asym{2.0 + (2.0i * sqrt3)};

    CHECK(asym_line.math_model_type() == ComponentType::branch);

    SUBCASE("Voltge error") { CHECK_THROWS_AS(AsymLine(input, 50.0, 10.0e3, 50.0e3), ConflictVoltage); }

    SUBCASE("General") {
        CHECK(branch.from_node() == 2);
        CHECK(branch.to_node() == 3);
        CHECK(branch.from_status() == true);
        CHECK(branch.to_status() == true);
        CHECK(branch.edge_status() == true);
        CHECK(branch.status(BranchSide::from) == branch.from_status());
        CHECK(branch.status(BranchSide::to) == branch.to_status());
        CHECK(branch.base_i_from() == doctest::Approx(base_i));
        CHECK(branch.base_i_to() == doctest::Approx(base_i));
        CHECK(branch.phase_shift() == 0.0);
        CHECK(branch.phase_incidence() == BranchPhaseIncidence{});
    }

    SUBCASE("Symmetric parameters") {
        // double connected
        BranchCalcParam<symmetric_t> param = branch.calc_param<symmetric_t>();
        CHECK(cabs(param.yff() - yff1) < numerical_tolerance);
        CHECK(cabs(param.ytt() - yff1) < numerical_tolerance);
        CHECK(cabs(param.ytf() - yft1) < numerical_tolerance);
        CHECK(cabs(param.yft() - yft1) < numerical_tolerance);
        // to connected
        CHECK(branch.update(BranchUpdate{1, false, na_IntS}).topo);
        param = branch.calc_param<symmetric_t>();
        CHECK(cabs(param.yff() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.ytt() - ys1) < numerical_tolerance);
        CHECK(cabs(param.ytf() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.yft() - 0.0) < numerical_tolerance);
        // not connected
        CHECK(branch.set_status(na_IntS, false));
        param = branch.calc_param<symmetric_t>();
        CHECK(cabs(param.yff() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.ytt() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.ytf() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.yft() - 0.0) < numerical_tolerance);
        // not changing
        CHECK(!branch.set_status(false, false));
        // from connected
        CHECK(branch.set_status(true, na_IntS));
        param = branch.calc_param<symmetric_t>();
        CHECK(cabs(param.yff() - ys1) < numerical_tolerance);
        CHECK(cabs(param.ytt() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.ytf() - 0.0) < numerical_tolerance);
        CHECK(cabs(param.yft() - 0.0) < numerical_tolerance);
    }

    SUBCASE("Asymmetric parameters") {
        // double connected
        BranchCalcParam<asymmetric_t> param = asym_line.calc_param<asymmetric_t>();
        CHECK((cabs(param.yff() - ytt) < numerical_tolerance).all());
        CHECK((cabs(param.ytt() - ytt) < numerical_tolerance).all());
        CHECK((cabs(param.ytf() - (-y_series)) < numerical_tolerance).all());
        CHECK((cabs(param.yft() - (-y_series)) < numerical_tolerance).all());
        // no source
        param = branch.calc_param<asymmetric_t>(false);
        CHECK((cabs(param.yff() - 0.0) < numerical_tolerance).all());
        CHECK((cabs(param.ytt() - 0.0) < numerical_tolerance).all());
        CHECK((cabs(param.ytf() - 0.0) < numerical_tolerance).all());
        CHECK((cabs(param.yft() - 0.0) < numerical_tolerance).all());
        // from connected
        CHECK(branch.set_status(na_IntS, false));
        param = asym_line.calc_param<asymmetric_t>();
        CHECK((cabs(param.yff() - branch_shunt) < numerical_tolerance).all()); // Fail
        CHECK((cabs(param.ytt() - 0.0) < numerical_tolerance).all());
        CHECK((cabs(param.ytf() - 0.0) < numerical_tolerance).all());
        CHECK((cabs(param.yft() - 0.0) < numerical_tolerance).all());
    }

    SUBCASE("Symmetric results with direct power and current output") {
        BranchSolverOutput<symmetric_t> branch_solver_output{};
        branch_solver_output.i_f = 1.0 - 2.0i;
        branch_solver_output.i_t = 2.0 - 1.0i;
        branch_solver_output.s_f = 1.0 - 1.5i;
        branch_solver_output.s_t = 1.5 - 1.5i;
        BranchOutput<symmetric_t> const output = branch.get_output<symmetric_t>(branch_solver_output);
        CHECK(output.id == 1);
        CHECK(output.energized);
        CHECK(output.loading == doctest::Approx(cabs(2.0 - 1.0i) * base_i / input.i_n));
        CHECK(output.i_from == doctest::Approx(cabs(1.0 - 2.0i) * base_i));
        CHECK(output.i_to == doctest::Approx(cabs(2.0 - 1.0i) * base_i));
        CHECK(output.s_from == doctest::Approx(cabs(1.0 - 1.5i) * base_power<symmetric_t>));
        CHECK(output.s_to == doctest::Approx(cabs(1.5 - 1.5i) * base_power<symmetric_t>));
        CHECK(output.p_from == doctest::Approx(1.0 * base_power<symmetric_t>));
        CHECK(output.p_to == doctest::Approx(1.5 * base_power<symmetric_t>));
        CHECK(output.q_from == doctest::Approx(-1.5 * base_power<symmetric_t>));
        CHECK(output.q_to == doctest::Approx(-1.5 * base_power<symmetric_t>));
    }

    SUBCASE("No source results") {
        BranchOutput<asymmetric_t> const output = branch.get_null_output<asymmetric_t>();
        CHECK(output.id == 1);
        CHECK(!output.energized);
        CHECK(output.loading == 0.0);
        CHECK(output.i_from(0) == 0.0);
        CHECK(output.i_to(1) == 0.0);
        CHECK(output.s_from(2) == 0.0);
        CHECK(output.s_to(0) == 0.0);
        CHECK(output.p_from(1) == 0.0);
        CHECK(output.p_to(2) == 0.0);
        CHECK(output.q_from(0) == 0.0);
        CHECK(output.q_to(1) == 0.0);
    }

    SUBCASE("No source short circuit results") {
        BranchShortCircuitOutput const output = branch.get_null_sc_output();
        CHECK(output.id == 1);
        CHECK(!output.energized);
        CHECK(output.i_from(0) == 0.0);
        CHECK(output.i_to(1) == 0.0);
        CHECK(output.i_from_angle(0) == 0.0);
        CHECK(output.i_to_angle(1) == 0.0);
    }

    SUBCASE("Asym short circuit results") {
        BranchShortCircuitOutput const asym_output = branch.get_sc_output(if_sc_asym, it_sc_asym);
        CHECK(asym_output.id == 1);
        CHECK(asym_output.energized);
        CHECK(asym_output.i_from(1) == doctest::Approx(cabs(if_sc) * base_i));
        CHECK(asym_output.i_from(2) == doctest::Approx(cabs(if_sc) * base_i));
        CHECK(asym_output.i_to(0) == doctest::Approx(cabs(it_sc) * base_i));
        CHECK(asym_output.i_to(1) == doctest::Approx(cabs(it_sc) * base_i));
        CHECK(asym_output.i_from_angle(0) == doctest::Approx(pi / 4));
        CHECK(asym_output.i_from_angle(2) == doctest::Approx(pi / 4 + deg_120));
        CHECK(asym_output.i_to_angle(1) == doctest::Approx(pi / 3 - deg_120));
        CHECK(asym_output.i_to_angle(2) == doctest::Approx(pi / 3 + deg_120));
        CHECK(asym_output.id == 1);
    }

    SUBCASE("Sym short circuit results") {
        BranchShortCircuitOutput const sym_output = branch.get_sc_output(if_sc, it_sc);
        BranchShortCircuitOutput const asym_output = branch.get_sc_output(if_sc_asym, it_sc_asym);
        CHECK(sym_output.energized == asym_output.energized);
        CHECK(sym_output.i_from(1) == doctest::Approx(asym_output.i_from(1)));
        CHECK(sym_output.i_from(2) == doctest::Approx(asym_output.i_from(2)));
        CHECK(sym_output.i_to(0) == doctest::Approx(asym_output.i_to(0)));
        CHECK(sym_output.i_to(1) == doctest::Approx(asym_output.i_to(1)));
        CHECK(sym_output.i_from_angle(0) == doctest::Approx(asym_output.i_from_angle(0)));
        CHECK(sym_output.i_from_angle(2) == doctest::Approx(asym_output.i_from_angle(2)));
        CHECK(sym_output.i_to_angle(1) == doctest::Approx(asym_output.i_to_angle(1)));
        CHECK(sym_output.i_to_angle(2) == doctest::Approx(asym_output.i_to_angle(2)));
    }

    SUBCASE("Update inverse") {
        BranchUpdate branch_update{.id = 1, .from_status = na_IntS, .to_status = na_IntS};
        auto expected = branch_update;

        SUBCASE("Identical") {
            // default values
        }

        SUBCASE("From status") {
            SUBCASE("same") { branch_update.from_status = status_to_int(asym_line.from_status()); }
            SUBCASE("different") { branch_update.from_status = IntS{0}; }
            expected.from_status = status_to_int(asym_line.from_status());
        }

        SUBCASE("To status") {
            SUBCASE("same") { branch_update.to_status = status_to_int(asym_line.to_status()); }
            SUBCASE("different") { branch_update.to_status = IntS{0}; }
            expected.to_status = status_to_int(asym_line.to_status());
        }

        SUBCASE("multiple") {
            branch_update.from_status = IntS{0};
            branch_update.to_status = IntS{0};
            expected.from_status = status_to_int(asym_line.from_status());
            expected.to_status = status_to_int(asym_line.to_status());
        }

        auto const inv = asym_line.inverse(branch_update);

        CHECK(inv.id == expected.id);
        CHECK(inv.from_status == expected.from_status);
        CHECK(inv.to_status == expected.to_status);
    }
}
} // namespace

TEST_CASE("Test asym line") {
    double const system_frequency = 50.0;
    double const voltage_lvl = 10.0e3;
    double const base_i = base_power_1p / (voltage_lvl / sqrt3);
    double const base_y = base_i * base_i / base_power_1p;

    SUBCASE("R and X matrix c0, c1 including neutral") {
        AsymLineInput const input = {.id = 1,
                                     .from_node = 2,
                                     .to_node = 3,
                                     .from_status = 1,
                                     .to_status = 1,
                                     .r_aa = 0.4369,
                                     .r_ba = 0.0496,
                                     .r_bb = 0.4369,
                                     .r_ca = 0.0485,
                                     .r_cb = 0.0496,
                                     .r_cc = 0.4369,
                                     .r_na = 0.0496,
                                     .r_nb = 0.0485,
                                     .r_nc = 0.0496,
                                     .r_nn = 0.4369,
                                     .x_aa = 0.8538,
                                     .x_ba = 0.7886,
                                     .x_bb = 0.8538,
                                     .x_ca = 0.7663,
                                     .x_cb = 0.7886,
                                     .x_cc = 0.8538,
                                     .x_na = 0.7886,
                                     .x_nb = 0.7663,
                                     .x_nc = 0.7886,
                                     .x_nn = 0.8538,
                                     .c0 = 0.18,
                                     .c1 = 0.308,
                                     .i_n = 216.0};

        ComplexTensor4 const r_matrix{input.r_aa, input.r_bb, input.r_cc, input.r_nn, input.r_ba,
                                      input.r_ca, input.r_na, input.r_cb, input.r_nb, input.r_nc};
        ComplexTensor4 const x_matrix{input.x_aa, input.x_bb, input.x_cc, input.x_nn, input.x_ba,
                                      input.x_ca, input.x_na, input.x_cb, input.x_nb, input.x_nc};
        auto const c_matrix =
            ComplexTensor<asymmetric_t>{(2.0 * input.c1 + input.c0) / 3.0, (input.c0 - input.c1) / 3.0};
        ComplexTensor4 const z = r_matrix + 1.0i * x_matrix;

        ComplexTensor<asymmetric_t> const y_series = 1 / base_y * inv(kron_reduction(z));
        execute_subcases(input, y_series, c_matrix, base_i, base_y, system_frequency, voltage_lvl);
    }

    SUBCASE("R and X matrix, c0, c1 excluding neutral") {
        AsymLineInput const input = {.id = 1,
                                     .from_node = 2,
                                     .to_node = 3,
                                     .from_status = 1,
                                     .to_status = 1,
                                     .r_aa = 0.4369,
                                     .r_ba = 0.0496,
                                     .r_bb = 0.4369,
                                     .r_ca = 0.0485,
                                     .r_cb = 0.0496,
                                     .r_cc = 0.4369,
                                     .x_aa = 0.8538,
                                     .x_ba = 0.7886,
                                     .x_bb = 0.8538,
                                     .x_ca = 0.7663,
                                     .x_cb = 0.7886,
                                     .x_cc = 0.8538,
                                     .c0 = 0.18,
                                     .c1 = 0.308,
                                     .i_n = 216.0};
        ComplexTensor<asymmetric_t> const r_matrix{input.r_aa, input.r_bb, input.r_cc,
                                                   input.r_ba, input.r_ca, input.r_cb};
        ComplexTensor<asymmetric_t> const x_matrix{input.x_aa, input.x_bb, input.x_cc,
                                                   input.x_ba, input.x_ca, input.x_cb};
        ComplexTensor<asymmetric_t> const y_series = 1 / base_y * inv(r_matrix + 1.0i * x_matrix);
        auto const c_matrix =
            ComplexTensor<asymmetric_t>{(2.0 * input.c1 + input.c0) / 3.0, (input.c0 - input.c1) / 3.0};
        execute_subcases(input, y_series, c_matrix, base_i, base_y, system_frequency, voltage_lvl);
    }

    SUBCASE("R, X excluding neutral C matrix") {
        AsymLineInput const input = {.id = 1,
                                     .from_node = 2,
                                     .to_node = 3,
                                     .from_status = 1,
                                     .to_status = 1,
                                     .r_aa = 0.4369,
                                     .r_ba = 0.0496,
                                     .r_bb = 0.4369,
                                     .r_ca = 0.0485,
                                     .r_cb = 0.0496,
                                     .r_cc = 0.4369,
                                     .x_aa = 0.8538,
                                     .x_ba = 0.7886,
                                     .x_bb = 0.8538,
                                     .x_ca = 0.7663,
                                     .x_cb = 0.7886,
                                     .x_cc = 0.8538,
                                     .c_aa = 0.3200,
                                     .c_ba = 0.5400,
                                     .c_bb = 0.3200,
                                     .c_ca = 0.7600,
                                     .c_cb = 0.5400,
                                     .c_cc = 0.3200,
                                     .i_n = 216.0};
        ComplexTensor<asymmetric_t> const r_matrix{input.r_aa, input.r_bb, input.r_cc,
                                                   input.r_ba, input.r_ca, input.r_cb};
        ComplexTensor<asymmetric_t> const x_matrix{input.x_aa, input.x_bb, input.x_cc,
                                                   input.x_ba, input.x_ca, input.x_cb};
        ComplexTensor<asymmetric_t> const y_series = 1 / base_y * inv(r_matrix + 1.0i * x_matrix);
        auto const c_matrix =
            ComplexTensor<asymmetric_t>(input.c_aa, input.c_bb, input.c_cc, input.c_ba, input.c_ca, input.c_cb);
        execute_subcases(input, y_series, c_matrix, base_i, base_y, system_frequency, voltage_lvl);
    }
}

namespace {
// conductor data a, b, c, n of the existing test cases, symmetric
using ConductorMatrix = std::array<std::array<double, 4>, 4>;
constexpr ConductorMatrix r_conductors{{{0.4369, 0.0496, 0.0485, 0.0496},
                                        {0.0496, 0.4369, 0.0496, 0.0485},
                                        {0.0485, 0.0496, 0.4369, 0.0496},
                                        {0.0496, 0.0485, 0.0496, 0.4569}}};
constexpr ConductorMatrix x_conductors{{{0.8538, 0.7886, 0.7663, 0.7886},
                                        {0.7886, 0.8538, 0.7886, 0.7663},
                                        {0.7663, 0.7886, 0.8538, 0.7886},
                                        {0.7886, 0.7663, 0.7886, 0.8738}}};
constexpr ConductorMatrix c_conductors{{{0.32e-6, -0.05e-6, -0.04e-6, 0.0},
                                        {-0.05e-6, 0.33e-6, -0.05e-6, 0.0},
                                        {-0.04e-6, -0.05e-6, 0.31e-6, 0.0},
                                        {0.0, 0.0, 0.0, 0.0}}};
constexpr double reduced_frequency = 50.0;
constexpr double reduced_voltage = 10.0e3;
constexpr double reduced_tolerance = 1e-12;

struct ReducedPhaseSet {
    PhaseConnection connection;
    PhaseMask mask;
    std::vector<Idx> phases;
};

std::vector<ReducedPhaseSet> reduced_phase_sets() {
    using enum PhaseConnection;
    return {{.connection = a, .mask = phase_mask_a, .phases = {0}},
            {.connection = b, .mask = phase_mask_b, .phases = {1}},
            {.connection = c, .mask = phase_mask_c, .phases = {2}},
            {.connection = ab, .mask = phase_mask_ab, .phases = {0, 1}},
            {.connection = ac, .mask = phase_mask_ac, .phases = {0, 2}},
            {.connection = bc, .mask = phase_mask_bc, .phases = {1, 2}}};
}

// input with values only for the present phases (and the neutral); the absent entries are left empty
AsymLineInput reduced_phase_input(ReducedPhaseSet const& phase_set, bool neutral) {
    AsymLineInput input{
        .id = 1, .from_node = 2, .to_node = 3, .from_status = 1, .to_status = 1, .phases = phase_set.connection};
    input.i_n = 216.0;
    auto const is_present = [&phase_set, neutral](Idx conductor) {
        return conductor == 3 ? neutral : phase_set.mask.has(conductor);
    };
    auto const set_rx = [&is_present](double& r, double& x, Idx row, Idx col) {
        if (is_present(row) && is_present(col)) {
            r = r_conductors[static_cast<size_t>(row)][static_cast<size_t>(col)];
            x = x_conductors[static_cast<size_t>(row)][static_cast<size_t>(col)];
        }
    };
    set_rx(input.r_aa, input.x_aa, 0, 0);
    set_rx(input.r_ba, input.x_ba, 1, 0);
    set_rx(input.r_bb, input.x_bb, 1, 1);
    set_rx(input.r_ca, input.x_ca, 2, 0);
    set_rx(input.r_cb, input.x_cb, 2, 1);
    set_rx(input.r_cc, input.x_cc, 2, 2);
    set_rx(input.r_na, input.x_na, 3, 0);
    set_rx(input.r_nb, input.x_nb, 3, 1);
    set_rx(input.r_nc, input.x_nc, 3, 2);
    set_rx(input.r_nn, input.x_nn, 3, 3);
    auto const set_c = [&phase_set](double& c, Idx row, Idx col) {
        if (phase_set.mask.has(row) && phase_set.mask.has(col)) {
            c = c_conductors[static_cast<size_t>(row)][static_cast<size_t>(col)];
        }
    };
    set_c(input.c_aa, 0, 0);
    set_c(input.c_ba, 1, 0);
    set_c(input.c_bb, 1, 1);
    set_c(input.c_ca, 2, 0);
    set_c(input.c_cb, 2, 1);
    set_c(input.c_cc, 2, 2);
    return input;
}

// a 1x1 or 2x2 block of the present phases, with closed-form operations, independent of the code under test
struct Block {
    size_t n{};
    std::array<std::array<DoubleComplex, 2>, 2> value{};

    friend Block operator+(Block const& x, Block const& y) {
        Block result{.n = x.n};
        for (size_t i = 0; i != x.n; ++i) {
            for (size_t j = 0; j != x.n; ++j) {
                result.value[i][j] = x.value[i][j] + y.value[i][j];
            }
        }
        return result;
    }
    friend Block operator*(DoubleComplex factor, Block const& x) {
        Block result{.n = x.n};
        for (size_t i = 0; i != x.n; ++i) {
            for (size_t j = 0; j != x.n; ++j) {
                result.value[i][j] = factor * x.value[i][j];
            }
        }
        return result;
    }
};

Block inverse(Block const& x) {
    if (x.n == 1) {
        return {.n = 1, .value = {{{1.0 / x.value[0][0], 0.0}, {0.0, 0.0}}}};
    }
    DoubleComplex const det = x.value[0][0] * x.value[1][1] - x.value[0][1] * x.value[1][0];
    return {.n = 2,
            .value = {{{x.value[1][1] / det, -x.value[0][1] / det}, {-x.value[1][0] / det, x.value[0][0] / det}}}};
}

// series impedance of the present phases, with the neutral eliminated (Z_PP - Z_Pn Z_nP / Z_nn)
Block expected_z_series(std::vector<Idx> const& phases, bool neutral) {
    auto const z = [](Idx row, Idx col) {
        auto const r = static_cast<size_t>(row);
        auto const c = static_cast<size_t>(col);
        return DoubleComplex{r_conductors[r][c], x_conductors[r][c]};
    };
    Block result{.n = phases.size()};
    for (size_t i = 0; i != phases.size(); ++i) {
        for (size_t j = 0; j != phases.size(); ++j) {
            result.value[i][j] = z(phases[i], phases[j]);
            if (neutral) {
                result.value[i][j] -= z(phases[i], 3) * z(3, phases[j]) / z(3, 3);
            }
        }
    }
    return result;
}

Block expected_c(std::vector<Idx> const& phases) {
    Block result{.n = phases.size()};
    for (size_t i = 0; i != phases.size(); ++i) {
        for (size_t j = 0; j != phases.size(); ++j) {
            result.value[i][j] = c_conductors[static_cast<size_t>(phases[i])][static_cast<size_t>(phases[j])];
        }
    }
    return result;
}

void check_block(ComplexTensor<asymmetric_t> const& actual, Block const& expected, std::vector<Idx> const& phases) {
    CHECK(!actual.isNaN().any());
    for (Idx row = 0; row != 3; ++row) {
        for (Idx col = 0; col != 3; ++col) {
            CAPTURE(row);
            CAPTURE(col);
            auto const row_it = std::ranges::find(phases, row);
            auto const col_it = std::ranges::find(phases, col);
            if (row_it == phases.end() || col_it == phases.end()) {
                CHECK(actual(row, col) == 0.0);
                continue;
            }
            DoubleComplex const value =
                expected
                    .value[static_cast<size_t>(row_it - phases.begin())][static_cast<size_t>(col_it - phases.begin())];
            CHECK(cabs(actual(row, col) - value) < reduced_tolerance * std::max(1.0, cabs(value)));
        }
    }
}

void check_zero(BranchCalcParam<asymmetric_t> const& param) {
    for (auto const& value : param.value) {
        CHECK((value == ComplexTensor<asymmetric_t>{}).all());
    }
}

void check_reduced_phase_line(ReducedPhaseSet const& phase_set, bool neutral) {
    double const base_i = base_power_3p / reduced_voltage / sqrt3;
    double const base_y = base_i / (reduced_voltage / sqrt3);
    AsymLine asym_line{reduced_phase_input(phase_set, neutral), reduced_frequency, reduced_voltage, reduced_voltage};
    Branch& branch = asym_line;
    auto const& phases = phase_set.phases;

    Block const y_series = (1.0 / base_y) * inverse(expected_z_series(phases, neutral));
    Block const y_shunt = (1.0i * 2.0 * pi * reduced_frequency / base_y) * expected_c(phases);
    Block const y_tt = y_series + 0.5 * y_shunt;
    Block const y_ft = -1.0 * y_series;
    Block const branch_shunt = 0.5 * y_shunt + inverse(inverse(y_series) + 2.0 * inverse(y_shunt));

    CHECK(branch.phase_incidence() ==
          BranchPhaseIncidence{.from = phase_set.mask, .to = phase_set.mask, .coupling = PhaseCoupling::per_phase});
    CHECK_THROWS_AS(branch.calc_param<symmetric_t>(), ReducedPhaseNotSupported);

    // both sides connected
    BranchCalcParam<asymmetric_t> param = branch.calc_param<asymmetric_t>();
    check_block(param.yff(), y_tt, phases);
    check_block(param.ytt(), y_tt, phases);
    check_block(param.yft(), y_ft, phases);
    check_block(param.ytf(), y_ft, phases);
    // not energized
    check_zero(branch.calc_param<asymmetric_t>(false));
    // from side connected: open-end shunt 0.5 Ysh + (Zs + 2 Ysh^-1)^-1
    CHECK(branch.set_status(na_IntS, false));
    param = branch.calc_param<asymmetric_t>();
    check_block(param.yff(), branch_shunt, phases);
    check_zero({.value = {param.ytt(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
    // to side connected
    CHECK(branch.set_status(false, true));
    param = branch.calc_param<asymmetric_t>();
    check_block(param.ytt(), branch_shunt, phases);
    check_zero({.value = {param.yff(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
    // not connected
    CHECK(branch.set_status(false, false));
    check_zero(branch.calc_param<asymmetric_t>());
}

ComplexTensor<asymmetric_t> test_impedance() {
    ComplexTensor<asymmetric_t> result;
    for (Idx row = 0; row != 3; ++row) {
        for (Idx col = 0; col != 3; ++col) {
            auto const r = static_cast<size_t>(row);
            auto const c = static_cast<size_t>(col);
            result(row, col) = DoubleComplex{r_conductors[r][c], x_conductors[r][c]};
        }
    }
    return result;
}
} // namespace

TEST_CASE("Test masked_inv") {
    ComplexTensor<asymmetric_t> const z = test_impedance();

    SUBCASE("All phases is the plain inverse") { CHECK((masked_inv(z, phase_mask_abc) == inv(z)).all()); }

    SUBCASE("Single phase with nan on the absent entries") {
        ComplexTensor<asymmetric_t> z_b;
        z_b.setConstant(DoubleComplex{nan, nan});
        z_b(1, 1) = z(1, 1);
        ComplexTensor<asymmetric_t> const y = masked_inv(z_b, phase_mask_b);
        CHECK(!y.isNaN().any());
        CHECK(cabs(y(1, 1) - 1.0 / z(1, 1)) < reduced_tolerance);
        for (Idx row = 0; row != 3; ++row) {
            for (Idx col = 0; col != 3; ++col) {
                if (row != 1 || col != 1) {
                    CHECK(y(row, col) == 0.0);
                }
            }
        }
    }

    SUBCASE("Two phases against the closed-form inverse") {
        ComplexTensor<asymmetric_t> z_ab = z;
        z_ab(2, 2) = nan;
        z_ab(0, 2) = nan;
        z_ab(2, 1) = nan;
        Block const z_block{.n = 2, .value = {{{z(0, 0), z(0, 1)}, {z(1, 0), z(1, 1)}}}};
        check_block(masked_inv(z_ab, phase_mask_ab), inverse(z_block), {0, 1});
    }
}

TEST_CASE("Test asym line - reduced phases") {
    SUBCASE("All one- and two-phase sets, with and without a neutral") {
        for (auto const& phase_set : reduced_phase_sets()) {
            for (bool const neutral : {false, true}) {
                CAPTURE(phase_letters(phase_set.mask));
                CAPTURE(neutral);
                check_reduced_phase_line(phase_set, neutral);
            }
        }
    }

    SUBCASE("Not provided and abc give identical parameters") {
        AsymLineInput input = reduced_phase_input(
            {.connection = PhaseConnection::abc, .mask = phase_mask_abc, .phases = {0, 1, 2}}, true);
        AsymLine const line_abc{input, reduced_frequency, reduced_voltage, reduced_voltage};
        input.phases = AsymLineInput{}.phases; // not provided
        AsymLine const line_na{input, reduced_frequency, reduced_voltage, reduced_voltage};
        auto const param_na = line_na.calc_param<asymmetric_t>();
        auto const param_abc = line_abc.calc_param<asymmetric_t>();
        for (size_t k = 0; k != param_na.value.size(); ++k) {
            CHECK((param_na.value[k] == param_abc.value[k]).all());
        }
        CHECK(line_abc.phase_incidence() == BranchPhaseIncidence{});
        CHECK(line_na.phase_incidence() == BranchPhaseIncidence{});
    }

    SUBCASE("SWER line") {
        double const r = 1.40 + (pi * pi * reduced_frequency * 1e-4);
        double const x = 4.0 * pi * reduced_frequency * 1e-4 * std::log(658.5 * std::sqrt(100.0 / 50.0) / 0.00127);
        AsymLineInput const input{.id = 1,
                                  .from_node = 2,
                                  .to_node = 3,
                                  .from_status = 1,
                                  .to_status = 1,
                                  .phases = PhaseConnection::a,
                                  .r_aa = 20.0 * r,
                                  .x_aa = 20.0 * x,
                                  .c_aa = 20.0 * 6e-9};
        double const base_y = base_power_3p / (22.0e3 * 22.0e3);
        AsymLine const swer{input, reduced_frequency, 22.0e3, 22.0e3};
        auto const param = swer.calc_param<asymmetric_t>();
        DoubleComplex const y_series = 1.0 / DoubleComplex{20.0 * r, 20.0 * x} / base_y;
        DoubleComplex const y_shunt = 1.0i * 2.0 * pi * reduced_frequency * 20.0 * 6e-9 / base_y;
        check_block(param.yft(), {.n = 1, .value = {{{-y_series, 0.0}, {0.0, 0.0}}}}, {0});
        check_block(param.yff(), {.n = 1, .value = {{{y_series + 0.5 * y_shunt, 0.0}, {0.0, 0.0}}}}, {0});
    }

    SUBCASE("Invalid input") {
        ReducedPhaseSet const phase_a = reduced_phase_sets().front();
        auto const construct = [](AsymLineInput const& input) {
            return AsymLine{input, reduced_frequency, reduced_voltage, reduced_voltage};
        };
        AsymLineInput input = reduced_phase_input(phase_a, false);
        REQUIRE_NOTHROW(construct(input));

        SUBCASE("Unknown phases") {
            // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
            input.phases = static_cast<PhaseConnection>(IntS{9});
            CHECK_THROWS_AS(construct(input), MissingCaseForEnumError);
        }
        SUBCASE("c0 and c1") {
            input.c_aa = nan;
            input.c0 = 0.2e-6;
            input.c1 = 0.3e-6;
            CHECK_THROWS_WITH_AS(construct(input),
                                 "Invalid phase connection for asym_line 1 with phases (a): c0 and c1 are only "
                                 "supported with all three phases; use the c matrix\n",
                                 InvalidPhaseConnection);
        }
        SUBCASE("Missing entry of a present phase") {
            input.r_aa = nan;
            CHECK_THROWS_WITH_AS(construct(input),
                                 "Invalid phase connection for asym_line 1 with phases (a): r_aa is required\n",
                                 InvalidPhaseConnection);
        }
        SUBCASE("Entry of an absent phase") {
            input.r_bb = 0.4;
            CHECK_THROWS_WITH_AS(construct(input),
                                 "Invalid phase connection for asym_line 1 with phases (a): r_bb refers to an absent "
                                 "phase or a missing neutral and must be left empty\n",
                                 InvalidPhaseConnection);
        }
        SUBCASE("Neutral entry of an absent phase") {
            input = reduced_phase_input(phase_a, true);
            REQUIRE_NOTHROW(construct(input));
            input.r_nb = 0.05;
            CHECK_THROWS_AS(construct(input), InvalidPhaseConnection);
        }
        SUBCASE("Neutral entries without the neutral self impedance") {
            input = reduced_phase_input(phase_a, true);
            input.r_nn = nan;
            input.x_nn = nan;
            CHECK_THROWS_WITH_AS(construct(input),
                                 "Invalid phase connection for asym_line 1 with phases (a): r_na refers to an absent "
                                 "phase or a missing neutral and must be left empty\n",
                                 InvalidPhaseConnection);
        }
        SUBCASE("Zero capacitance is allowed") {
            input.c_aa = 0.0;
            CHECK_NOTHROW(construct(input));
        }
    }
}

} // namespace power_grid_model
