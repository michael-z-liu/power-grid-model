// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

// Delta (phase-to-phase) load_gen math for the power-flow solvers
//
// The expected values come from closed forms that are derived by hand (balanced three-phase and split-phase voltages)
// and from central finite differences of phase_power, which computes U o conj(D^T I_delta) and does not share the
// per-phase power and derivative terms of power_and_jacobian.

#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/exception.hpp>
#include <power_grid_model/common/grouped_index_vector.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>
#include <power_grid_model/math_solver/delta_load_gen.hpp>

#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>

namespace power_grid_model::math_solver::delta_load_gen {
namespace {
using CValue = ComplexValue<asymmetric_t>;
using RTensor = RealTensor<asymmetric_t>;
using CTensor = ComplexTensor<asymmetric_t>;

constexpr double tol = 1e-14;
constexpr double fd_step = 1e-6;
constexpr double fd_tol = 1e-8;
constexpr std::array all_types{LoadGenType::const_pq, LoadGenType::const_i, LoadGenType::const_y};

double deg(double degrees) { return degrees * std::numbers::pi / 180.0; }

CValue polar_value(std::array<double, 3> const& magnitude, std::array<double, 3> const& angle_deg) {
    CValue result{};
    for (Idx phase = 0; phase != 3; ++phase) {
        auto const p = static_cast<std::size_t>(phase);
        result(phase) = std::polar(magnitude[p], deg(angle_deg[p]));
    }
    return result;
}

// 2x2 block in rows/columns a and b, zero in row c and column c
RTensor ab_block(double aa, double ab, double ba, double bb) {
    RTensor result{};
    result(0, 0) = aa;
    result(0, 1) = ab;
    result(1, 0) = ba;
    result(1, 1) = bb;
    return result;
}

void check_close(CValue const& actual, CValue const& expected, double tolerance) {
    for (Idx phase = 0; phase != 3; ++phase) {
        CAPTURE(phase);
        CHECK(std::abs(actual(phase) - expected(phase)) < tolerance);
    }
}

void check_close(RTensor const& actual, RTensor const& expected, double tolerance) {
    for (Idx row = 0; row != 3; ++row) {
        for (Idx col = 0; col != 3; ++col) {
            CAPTURE(row);
            CAPTURE(col);
            CHECK(std::abs(actual(row, col) - expected(row, col)) < tolerance);
        }
    }
}

// exactly +0, not -0 and not NaN
bool is_positive_zero(double x) { return x == 0.0 && !std::signbit(x); }
bool is_positive_zero(DoubleComplex const& x) { return is_positive_zero(x.real()) && is_positive_zero(x.imag()); }

void check_phase_zero(CValue const& value, Idx phase) { CHECK(is_positive_zero(value(phase))); }

template <typename TensorType> void check_row_col_zero(TensorType const& value, Idx phase) {
    for (Idx other = 0; other != 3; ++other) {
        CAPTURE(other);
        CHECK(is_positive_zero(value(phase, other)));
        CHECK(is_positive_zero(value(other, phase)));
    }
}

void check_jacobian_row_col_zero(PowerAndJacobian const& result, Idx phase) {
    check_phase_zero(result.s, phase);
    check_row_col_zero(result.dp_dtheta, phase);
    check_row_col_zero(result.dp_dv, phase);
    check_row_col_zero(result.dq_dtheta, phase);
    check_row_col_zero(result.dq_dv, phase);
}

// the derivatives of section 1 of the design for a single branch ab, evaluated from given alpha_a, alpha_b and S(U)
PowerAndJacobian expected_ab_jacobian(DoubleComplex const& s_u, DoubleComplex const& alpha_a,
                                      DoubleComplex const& alpha_b, int n) {
    PowerAndJacobian result{};
    std::array const alpha{alpha_a, alpha_b};
    DoubleComplex const beta = s_u * alpha_a * alpha_b;
    auto const n_real = static_cast<double>(n);
    for (Idx r = 0; r != 2; ++r) {
        DoubleComplex const g = s_u * alpha[static_cast<std::size_t>(r)];
        result.s(r) = g;
        for (Idx c = 0; c != 2; ++c) {
            DoubleComplex const alpha_c = alpha[static_cast<std::size_t>(c)];
            double const sigma = r == c ? 1.0 : -1.0;
            result.dp_dtheta(r, c) = -sigma * beta.imag() - n_real * g.real() * alpha_c.imag();
            result.dq_dtheta(r, c) = sigma * beta.real() - n_real * g.imag() * alpha_c.imag();
            result.dp_dv(r, c) = sigma * beta.real() + n_real * g.real() * alpha_c.real();
            result.dq_dv(r, c) = sigma * beta.imag() + n_real * g.imag() * alpha_c.real();
        }
    }
    return result;
}

void check_jacobian_close(PowerAndJacobian const& actual, PowerAndJacobian const& expected, double tolerance) {
    check_close(actual.s, expected.s, tolerance);
    check_close(actual.dp_dtheta, expected.dp_dtheta, tolerance);
    check_close(actual.dp_dv, expected.dp_dv, tolerance);
    check_close(actual.dq_dtheta, expected.dq_dtheta, tolerance);
    check_close(actual.dq_dv, expected.dq_dv, tolerance);
}

// central differences of phase_power with respect to theta_s and V_s (scaled by V_s), compared with the analytic
// Jacobian; returns the largest absolute deviation
double finite_difference_deviation(CValue const& s_delta, LoadGenType type, CValue const& u, PhaseMask phases,
                                   double u_nom) {
    PowerAndJacobian const analytic = power_and_jacobian(s_delta, type, u, phases, u_nom);
    double max_deviation = 0.0;
    for (Idx col = 0; col != 3; ++col) {
        CValue u_theta_plus = u;
        CValue u_theta_minus = u;
        u_theta_plus(col) *= std::polar(1.0, fd_step);
        u_theta_minus(col) *= std::polar(1.0, -fd_step);
        CValue u_v_plus = u;
        CValue u_v_minus = u;
        u_v_plus(col) *= 1.0 + fd_step;
        u_v_minus(col) *= 1.0 - fd_step;
        CValue const ds_dtheta = (phase_power(s_delta, type, u_theta_plus, phases, u_nom) -
                                  phase_power(s_delta, type, u_theta_minus, phases, u_nom)) /
                                 (2.0 * fd_step);
        CValue const v_ds_dv = (phase_power(s_delta, type, u_v_plus, phases, u_nom) -
                                phase_power(s_delta, type, u_v_minus, phases, u_nom)) /
                               (2.0 * fd_step);
        for (Idx row = 0; row != 3; ++row) {
            max_deviation = std::max({max_deviation, std::abs(ds_dtheta(row).real() - analytic.dp_dtheta(row, col)),
                                      std::abs(ds_dtheta(row).imag() - analytic.dq_dtheta(row, col)),
                                      std::abs(v_ds_dv(row).real() - analytic.dp_dv(row, col)),
                                      std::abs(v_ds_dv(row).imag() - analytic.dq_dv(row, col))});
        }
    }
    return max_deviation;
}

// rotation invariance (row sums over theta are zero) and Euler homogeneity (row sums of V d/dV are n times the power)
void check_invariants(CValue const& s_delta, LoadGenType type, CValue const& u, PhaseMask phases, double u_nom) {
    PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phases, u_nom);
    auto const n_real = static_cast<double>(voltage_exponent(type));
    for (Idx row = 0; row != 3; ++row) {
        CAPTURE(row);
        CHECK(std::abs(result.dp_dtheta.row(row).sum()) < 1e-13);
        CHECK(std::abs(result.dq_dtheta.row(row).sum()) < 1e-13);
        CHECK(std::abs(result.dp_dv.row(row).sum() - n_real * result.s(row).real()) < 1e-13);
        CHECK(std::abs(result.dq_dv.row(row).sum() - n_real * result.s(row).imag()) < 1e-13);
    }
}

// unbalanced test point of the design
CValue unbalanced_u() { return polar_value({1.02, 0.97, 1.01}, {2.0, -118.0, 121.0}); }
CValue unbalanced_s_delta() { return CValue{{0.3, 0.1}, {-0.2, 0.05}, {0.15, -0.07}}; }
} // namespace

TEST_CASE("Delta load_gen - nominal voltage, voltage exponent and active entries") {
    CHECK(nominal_voltage(SlotPattern::three_phase) == sqrt3);
    CHECK(nominal_voltage(SlotPattern::split_phase) == 2.0);
    CHECK_THROWS_AS(nominal_voltage(static_cast<SlotPattern>(7)), MissingCaseForEnumError);

    CHECK(voltage_exponent(LoadGenType::const_pq) == 0);
    CHECK(voltage_exponent(LoadGenType::const_i) == 1);
    CHECK(voltage_exponent(LoadGenType::const_y) == 2);
    CHECK_THROWS_AS(voltage_exponent(static_cast<LoadGenType>(7)), MissingCaseForEnumError);

    CHECK(delta_branches[0].from == 0);
    CHECK(delta_branches[0].to == 1);
    CHECK(delta_branches[0].phases == phase_mask_ab);
    CHECK(delta_branches[1].from == 1);
    CHECK(delta_branches[1].to == 2);
    CHECK(delta_branches[1].phases == phase_mask_bc);
    CHECK(delta_branches[2].from == 2);
    CHECK(delta_branches[2].to == 0);
    CHECK(delta_branches[2].phases == phase_mask_ac);

    DoubleComplex const s{0.1, 0.2};
    CHECK(is_active(s, phase_mask_abc, delta_branches[0]));
    CHECK(is_active(s, phase_mask_ab, delta_branches[0]));
    CHECK_FALSE(is_active(s, phase_mask_ab, delta_branches[1]));
    CHECK_FALSE(is_active(s, phase_mask_ab, delta_branches[2]));
    CHECK(is_active(s, phase_mask_ac, delta_branches[2]));
    CHECK_FALSE(is_active(s, phase_mask_a, delta_branches[0]));
    CHECK_FALSE(is_active(DoubleComplex{}, phase_mask_abc, delta_branches[0]));
    CHECK(is_active(DoubleComplex{nan, 0.0}, phase_mask_abc, delta_branches[0])); // NaN counts as specified
    CHECK(is_active(DoubleComplex{0.0, 0.1}, phase_mask_abc, delta_branches[0]));
}

TEST_CASE("Delta load_gen - MathModelTopology connection accessors") {
    MathModelTopology topo{};
    topo.load_gens_per_bus = SparseGroupedIdxVector{from_sparse, {0, 0, 3}};
    REQUIRE(topo.n_load_gen() == 3);
    CHECK_FALSE(topo.has_delta_load_gen());
    CHECK_FALSE(topo.is_delta_load_gen(0));
    CHECK_FALSE(topo.is_delta_load_gen(2));

    topo.load_gen_connection = {LoadGenConnection::wye, LoadGenConnection::delta, LoadGenConnection::wye};
    CHECK(topo.has_delta_load_gen());
    CHECK_FALSE(topo.is_delta_load_gen(0));
    CHECK(topo.is_delta_load_gen(1));
    CHECK_FALSE(topo.is_delta_load_gen(2));
}

TEST_CASE("Delta load_gen - closed form, balanced three-phase voltages") {
    // U = (1, a^2, a), only the ab entry, U_nom = sqrt(3) = |U_ab|, so S(U) = S for every type
    // alpha_a = 1/2 - j c, alpha_b = 1/2 + j c with c = 1 / (2 sqrt(3)); beta = S / 3
    CValue const u{DoubleComplex{1.0}, a2, a};
    double const u_nom = nominal_voltage(SlotPattern::three_phase);
    double const c = 1.0 / (2.0 * sqrt3);
    DoubleComplex const alpha_a{0.5, -c};
    DoubleComplex const alpha_b{0.5, c};

    SUBCASE("Real power") {
        double const p = 1.0;
        CValue const s_delta{DoubleComplex{p}, DoubleComplex{}, DoubleComplex{}};
        for (LoadGenType const type : all_types) {
            int const n = voltage_exponent(type);
            auto const n_real = static_cast<double>(n);
            CAPTURE(n);
            PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phase_mask_abc, u_nom);

            CValue const expected_s{DoubleComplex{p / 2.0, -p / (2.0 * sqrt3)},
                                    DoubleComplex{p / 2.0, p / (2.0 * sqrt3)}, DoubleComplex{}};
            check_close(result.s, expected_s, tol);
            check_close(phase_power(s_delta, type, u, phase_mask_abc, u_nom), expected_s, tol);

            CValue const i = phase_current(s_delta, type, u, phase_mask_abc, u_nom);
            DoubleComplex const i_a = std::polar(p / sqrt3, std::numbers::pi / 6.0);
            check_close(i, CValue{i_a, -i_a, DoubleComplex{}}, tol);
            check_phase_zero(i, 2);

            double const k = n_real * p / (4.0 * sqrt3);
            double const q_theta = (p / 3.0) - (n_real * p / 12.0);
            double const quarter = n_real * p / 4.0;
            check_close(result.dp_dtheta, ab_block(k, -k, k, -k), tol);
            check_close(result.dq_dtheta, ab_block(q_theta, -q_theta, -q_theta, q_theta), tol);
            check_close(result.dp_dv,
                        ab_block((p / 3.0) + quarter, (-p / 3.0) + quarter, (-p / 3.0) + quarter, (p / 3.0) + quarter),
                        tol);
            check_close(result.dq_dv, ab_block(-k, -k, k, k), tol);

            if (n == 1) { // numbers from the design
                CHECK(result.dp_dtheta(0, 0) == doctest::Approx(0.144338).epsilon(1e-6));
                CHECK(result.dq_dtheta(0, 0) == doctest::Approx(0.25).epsilon(1e-12));
                CHECK(result.dp_dv(0, 0) == doctest::Approx(0.583333).epsilon(1e-6));
                CHECK(result.dp_dv(0, 1) == doctest::Approx(-0.083333).epsilon(1e-5));
            }
            check_jacobian_row_col_zero(result, 2);
        }
    }

    SUBCASE("Complex power") {
        DoubleComplex const s{1.0, 0.5};
        CValue const s_delta{s, DoubleComplex{}, DoubleComplex{}};
        for (LoadGenType const type : all_types) {
            int const n = voltage_exponent(type);
            CAPTURE(n);
            PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phase_mask_abc, u_nom);
            check_jacobian_close(result, expected_ab_jacobian(s, alpha_a, alpha_b, n), tol);

            CValue const i = phase_current(s_delta, type, u, phase_mask_abc, u_nom);
            DoubleComplex const i_a = std::conj(s) * std::polar(1.0 / sqrt3, std::numbers::pi / 6.0);
            check_close(i, CValue{i_a, -i_a, DoubleComplex{}}, tol);
            check_phase_zero(i, 2);
            check_jacobian_row_col_zero(result, 2);
        }
    }
}

TEST_CASE("Delta load_gen - closed form, split-phase voltages") {
    // U = (1, -1, 0) at a split-phase bus (phases ab), U_nom = 2 = |U_ab|; alpha_a = alpha_b = 1/2, beta = S / 4
    CValue const u{DoubleComplex{1.0}, DoubleComplex{-1.0}, DoubleComplex{}};
    double const u_nom = nominal_voltage(SlotPattern::split_phase);
    double const p = 0.8;
    CValue const s_delta{DoubleComplex{p}, DoubleComplex{}, DoubleComplex{}};

    for (LoadGenType const type : all_types) {
        auto const n_real = static_cast<double>(voltage_exponent(type));
        CAPTURE(n_real);
        PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phase_mask_ab, u_nom);
        CValue const expected_s{DoubleComplex{p / 2.0}, DoubleComplex{p / 2.0}, DoubleComplex{}};
        check_close(result.s, expected_s, tol);
        check_close(phase_power(s_delta, type, u, phase_mask_ab, u_nom), expected_s, tol);

        double const quarter = p / 4.0;
        check_close(result.dp_dtheta, RTensor{}, tol);
        check_close(result.dq_dtheta, ab_block(quarter, -quarter, -quarter, quarter), tol);
        check_close(result.dp_dv,
                    ab_block(quarter + n_real * quarter, -quarter + n_real * quarter, -quarter + n_real * quarter,
                             quarter + n_real * quarter),
                    tol);
        check_close(result.dq_dv, RTensor{}, tol);
        check_jacobian_row_col_zero(result, 2);
    }

    SUBCASE("Voltage dependence") {
        // |U_ab| = 2.2 = 1.1 U_nom
        CValue const u_high{DoubleComplex{1.1}, DoubleComplex{-1.1}, DoubleComplex{}};
        auto const total = [&](LoadGenType type) {
            CValue const s = phase_power(s_delta, type, u_high, phase_mask_ab, u_nom);
            return s(0) + s(1) + s(2);
        };
        CHECK(std::abs(total(LoadGenType::const_pq) - p) < tol);
        CHECK(std::abs(total(LoadGenType::const_i) - 1.1 * p) < tol);
        CHECK(std::abs(total(LoadGenType::const_y) - 1.21 * p) < tol);
        CHECK(std::abs(branch_power(DoubleComplex{p}, 1, DoubleComplex{2.2}, u_nom) - 1.1 * p) < tol);
        CHECK(std::abs(branch_power(DoubleComplex{p}, 2, DoubleComplex{2.2}, u_nom) - 1.21 * p) < tol);
    }
}

TEST_CASE("Delta load_gen - power conservation") {
    CValue const u = unbalanced_u();
    CValue const s_delta = unbalanced_s_delta();
    for (double const u_nom : {sqrt3, 2.0}) {
        for (LoadGenType const type : all_types) {
            int const n = voltage_exponent(type);
            CAPTURE(u_nom);
            CAPTURE(n);
            CValue const s = phase_power(s_delta, type, u, phase_mask_abc, u_nom);
            // the branch voltages written out: ab, bc, ca
            DoubleComplex const expected_total = branch_power(s_delta(0), n, u(0) - u(1), u_nom) +
                                                 branch_power(s_delta(1), n, u(1) - u(2), u_nom) +
                                                 branch_power(s_delta(2), n, u(2) - u(0), u_nom);
            CHECK(std::abs(s(0) + s(1) + s(2) - expected_total) < tol);
            check_close(power_and_jacobian(s_delta, type, u, phase_mask_abc, u_nom).s, s, tol);
        }
    }
}

TEST_CASE("Delta load_gen - finite-difference check of the Jacobian") {
    SUBCASE("Three-phase bus, all entries") {
        CValue const u = unbalanced_u();
        CValue const s_delta = unbalanced_s_delta();
        for (double const u_nom : {sqrt3, 2.0}) {
            for (LoadGenType const type : all_types) {
                CAPTURE(u_nom);
                CAPTURE(voltage_exponent(type));
                CHECK(finite_difference_deviation(s_delta, type, u, phase_mask_abc, u_nom) < fd_tol);
                check_invariants(s_delta, type, u, phase_mask_abc, u_nom);
            }
        }
    }
    SUBCASE("Split-phase bus, ab entry") {
        CValue const u = polar_value({1.01, 0.99, 0.0}, {3.0, -178.0, 0.0});
        CValue const s_delta{DoubleComplex{0.4, 0.1}, DoubleComplex{}, DoubleComplex{}};
        double const u_nom = nominal_voltage(SlotPattern::split_phase);
        for (LoadGenType const type : all_types) {
            CAPTURE(voltage_exponent(type));
            CHECK(finite_difference_deviation(s_delta, type, u, phase_mask_ab, u_nom) < fd_tol);
            check_invariants(s_delta, type, u, phase_mask_ab, u_nom);
            check_jacobian_row_col_zero(power_and_jacobian(s_delta, type, u, phase_mask_ab, u_nom), 2);
        }
    }
}

TEST_CASE("Delta load_gen - inactive entries contribute exactly zero") {
    double const u_nom = sqrt3;

    SUBCASE("Two-phase bus: the bc and ca entries are ignored") {
        CValue const u = polar_value({1.02, 0.97, 0.0}, {2.0, -118.0, 0.0});
        CValue const s_all = unbalanced_s_delta();
        CValue const s_ab_only{s_all(0), DoubleComplex{}, DoubleComplex{}};
        for (LoadGenType const type : all_types) {
            CAPTURE(voltage_exponent(type));
            PowerAndJacobian const all = power_and_jacobian(s_all, type, u, phase_mask_ab, u_nom);
            PowerAndJacobian const ab_only = power_and_jacobian(s_ab_only, type, u, phase_mask_ab, u_nom);
            CHECK((all.s == ab_only.s).all());
            CHECK((all.dp_dtheta == ab_only.dp_dtheta).all());
            CHECK((all.dp_dv == ab_only.dp_dv).all());
            CHECK((all.dq_dtheta == ab_only.dq_dtheta).all());
            CHECK((all.dq_dv == ab_only.dq_dv).all());
            check_jacobian_row_col_zero(all, 2);

            CValue const i = phase_current(s_all, type, u, phase_mask_ab, u_nom);
            CHECK((i == phase_current(s_ab_only, type, u, phase_mask_ab, u_nom)).all());
            check_phase_zero(i, 2);
            check_phase_zero(phase_power(s_all, type, u, phase_mask_ab, u_nom), 2);
        }
        CTensor const y = linear_admittance(s_all, phase_mask_ab, u_nom);
        CHECK((y == linear_admittance(s_ab_only, phase_mask_ab, u_nom)).all());
        check_row_col_zero(y, 2);
    }

    SUBCASE("One-phase bus: every entry is ignored") {
        CValue const u{DoubleComplex{1.0}, DoubleComplex{}, DoubleComplex{}};
        DoubleComplex const s{0.3, 0.1};
        CValue const s_delta{s, s, s};
        for (LoadGenType const type : all_types) {
            CAPTURE(voltage_exponent(type));
            PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phase_mask_a, u_nom);
            CValue const i = phase_current(s_delta, type, u, phase_mask_a, u_nom);
            CValue const s_phase = phase_power(s_delta, type, u, phase_mask_a, u_nom);
            for (Idx phase = 0; phase != 3; ++phase) {
                CAPTURE(phase);
                check_jacobian_row_col_zero(result, phase);
                check_phase_zero(i, phase);
                check_phase_zero(s_phase, phase);
            }
        }
        CTensor const y = linear_admittance(s_delta, phase_mask_a, u_nom);
        for (Idx phase = 0; phase != 3; ++phase) {
            check_row_col_zero(y, phase);
        }
    }

    SUBCASE("An entry of exactly zero is skipped") {
        // U_b == U_c, so U_bc = 0; the bc entry is zero and must not be evaluated
        CValue const u{DoubleComplex{1.0}, DoubleComplex{-0.5, -0.2}, DoubleComplex{-0.5, -0.2}};
        CValue const s_delta{DoubleComplex{0.3, 0.1}, DoubleComplex{}, DoubleComplex{0.15, -0.07}};
        for (LoadGenType const type : all_types) {
            CAPTURE(voltage_exponent(type));
            PowerAndJacobian const result = power_and_jacobian(s_delta, type, u, phase_mask_abc, u_nom);
            CHECK(result.s.isFinite().all());
            CHECK(result.dp_dtheta.isFinite().all());
            CHECK(result.dp_dv.isFinite().all());
            CHECK(result.dq_dtheta.isFinite().all());
            CHECK(result.dq_dv.isFinite().all());
            CHECK(phase_current(s_delta, type, u, phase_mask_abc, u_nom).isFinite().all());
            CHECK(phase_power(s_delta, type, u, phase_mask_abc, u_nom).isFinite().all());
        }
    }
}

TEST_CASE("Delta load_gen - linear admittance") {
    SUBCASE("Closed form") {
        DoubleComplex const s{0.6, -0.2};
        CValue const s_delta{s, DoubleComplex{}, DoubleComplex{}};
        for (double const u_nom : {sqrt3, 2.0}) {
            CAPTURE(u_nom);
            DoubleComplex const y_k = -std::conj(s) / (u_nom * u_nom);
            CTensor expected{};
            expected(0, 0) = y_k;
            expected(0, 1) = -y_k;
            expected(1, 0) = -y_k;
            expected(1, 1) = y_k;
            CTensor const y = linear_admittance(s_delta, phase_mask_abc, u_nom);
            for (Idx row = 0; row != 3; ++row) {
                for (Idx col = 0; col != 3; ++col) {
                    CAPTURE(row);
                    CAPTURE(col);
                    CHECK(std::abs(y(row, col) - expected(row, col)) < tol);
                }
            }
            check_row_col_zero(y, 2);
        }
    }

    SUBCASE("Constant admittance equals the const_y current") {
        CValue const u = unbalanced_u();
        CValue const s_delta = unbalanced_s_delta();
        for (double const u_nom : {sqrt3, 2.0}) {
            CAPTURE(u_nom);
            CTensor const y = linear_admittance(s_delta, phase_mask_abc, u_nom);
            CValue const i = phase_current(s_delta, LoadGenType::const_y, u, phase_mask_abc, u_nom);
            for (Idx row = 0; row != 3; ++row) {
                DoubleComplex y_u{};
                for (Idx col = 0; col != 3; ++col) {
                    y_u += y(row, col) * u(col);
                }
                CAPTURE(row);
                CHECK(std::abs(y_u + i(row)) < tol);
            }
        }
    }
}

} // namespace power_grid_model::math_solver::delta_load_gen
