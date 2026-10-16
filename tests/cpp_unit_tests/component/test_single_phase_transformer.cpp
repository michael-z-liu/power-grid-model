// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/component/single_phase_transformer.hpp>

#include <power_grid_model/auxiliary/input.hpp>
#include <power_grid_model/auxiliary/output.hpp>
#include <power_grid_model/auxiliary/update.hpp>
#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/exception.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>
#include <power_grid_model/component/component.hpp>
#include <power_grid_model/component/line_utils.hpp>

#include <doctest/doctest.h>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace power_grid_model {

using namespace std::complex_literals;

namespace {
// the expected values below are computed from the closed forms of the design (C^T Y_prim C with 2x2 algebra only),
// independently of the nodal stamp and the Kron elimination of the code under test
constexpr double test_tolerance = 1e-12;
constexpr double test_u_rated_from = 11.0e3;
constexpr double test_u_rated_to = 415.0;
constexpr double test_u1 = 11.0e3;
constexpr double test_u2 = 250.0;
constexpr double test_sn = 50.0e3;
constexpr double test_uk = 0.04;
constexpr double test_pk = 600.0;
constexpr double test_i0 = 0.015;
constexpr double test_p0 = 120.0;
constexpr double test_base_power_1p = 1.0e6 / 3.0;

struct ConnectionCase {
    WindingConnection1P connection;
    PhaseMask mask;
    // the winding voltage is the sum of incidence[k] * U_k over the slots k = a, b, c
    std::array<double, 3> incidence;

    bool phase_to_phase() const {
        return std::ranges::count_if(incidence, [](double value) { return value != 0.0; }) == 2;
    }
    // angle of the winding voltage for a balanced positive-sequence set U_a = 1, U_b = a^2, U_c = a
    double angle() const {
        std::array<DoubleComplex, 3> const balanced{1.0, std::polar(1.0, -2.0 * pi / 3.0),
                                                    std::polar(1.0, 2.0 * pi / 3.0)};
        DoubleComplex winding_voltage{};
        for (size_t k = 0; k != 3; ++k) {
            winding_voltage += incidence[k] * balanced[k];
        }
        return std::arg(winding_voltage);
    }
};

std::vector<ConnectionCase> connection_cases() {
    using enum WindingConnection1P;
    return {{.connection = a, .mask = phase_mask_a, .incidence = {1.0, 0.0, 0.0}},
            {.connection = b, .mask = phase_mask_b, .incidence = {0.0, 1.0, 0.0}},
            {.connection = c, .mask = phase_mask_c, .incidence = {0.0, 0.0, 1.0}},
            {.connection = ab, .mask = phase_mask_ab, .incidence = {1.0, -1.0, 0.0}},
            {.connection = ac, .mask = phase_mask_ac, .incidence = {1.0, 0.0, -1.0}},
            {.connection = bc, .mask = phase_mask_bc, .incidence = {0.0, 1.0, -1.0}}};
}

ConnectionCase connection_case(WindingConnection1P connection) {
    auto const cases = connection_cases();
    return *std::ranges::find(cases, connection, &ConnectionCase::connection);
}

SinglePhaseTransformerInput test_input(WindingConnection1P from_connection, WindingConnection1P to_connection) {
    return SinglePhaseTransformerInput{.id = 1,
                                       .from_node = 2,
                                       .to_node = 3,
                                       .from_status = 1,
                                       .to_status = 1,
                                       .u1 = test_u1,
                                       .u2 = test_u2,
                                       .sn = test_sn,
                                       .uk = test_uk,
                                       .pk = test_pk,
                                       .i0 = test_i0,
                                       .p0 = test_p0,
                                       .from_connection = from_connection,
                                       .to_connection = to_connection};
}

SinglePhaseTransformer make_transformer(SinglePhaseTransformerInput const& input) {
    return SinglePhaseTransformer{input, test_u_rated_from, test_u_rated_to};
}

// 2x2 complex block with closed-form operations
using Block2 = std::array<std::array<DoubleComplex, 2>, 2>;

Block2 inverse(Block2 const& x) {
    DoubleComplex const det = x[0][0] * x[1][1] - x[0][1] * x[1][0];
    return {{{x[1][1] / det, -x[0][1] / det}, {-x[1][0] / det, x[0][0] / det}}};
}

Block2 add_diagonal(Block2 x, DoubleComplex d0, DoubleComplex d1) {
    x[0][0] += d0;
    x[1][1] += d1;
    return x;
}

struct WindingData {
    double u1{test_u1};
    double u2{test_u2};
    double uk{test_uk};
    double pk{test_pk};
    double i0{test_i0};
    double p0{test_p0};
};

DoubleComplex series_admittance(WindingData const& data) {
    double const z_abs = data.uk * data.u2 * data.u2 / test_sn;
    double const r = data.pk * data.u2 * data.u2 / test_sn / test_sn;
    return 1.0 / DoubleComplex{r, std::sqrt(z_abs * z_abs - r * r)};
}

DoubleComplex magnetising_admittance(WindingData const& data) {
    double const y_abs = data.i0 * test_sn / data.u2 / data.u2;
    double const g = data.p0 / data.u2 / data.u2;
    return {g, -std::sqrt(y_abs * y_abs - g * g)};
}

// primitive admittance in siemens: y_s [[1/n^2, -1/n], [-1/n, 1]] + y_m / 2 diag(1/n^2, 1)
Block2 primitive(WindingData const& data) {
    DoubleComplex const y_s = series_admittance(data);
    DoubleComplex const y_m = magnetising_admittance(data);
    double const n = data.u1 / data.u2;
    return {{{(y_s + 0.5 * y_m) / (n * n), -y_s / n}, {-y_s / n, y_s + 0.5 * y_m}}};
}

// open-end shunt of the connected winding, referred to side 2
DoubleComplex open_end_shunt(WindingData const& data) {
    DoubleComplex const y_s = series_admittance(data);
    DoubleComplex const y_m = magnetising_admittance(data);
    return 0.5 * y_m + 1.0 / (1.0 / y_s + 2.0 / y_m);
}

struct ExpectedCase {
    ConnectionCase from;
    ConnectionCase to;
    Block2 y{};           // admittance of the windings in siemens, earthing included
    double polarity{1.0}; // -1 for clock 6
    bool from_connected{true};
    bool to_connected{true};
    double u_rated_from{test_u_rated_from};
    double u_rated_to{test_u_rated_to};
};

// C^T Y C scaled to per unit, plus the low susceptance y_low / 2 on each connected phase-to-phase side
BranchCalcParam<asymmetric_t> expected_param(ExpectedCase const& e) {
    double const base_u_from = e.u_rated_from / sqrt3;
    double const base_u_to = e.u_rated_to / sqrt3;
    double const s_ff = base_u_from * base_u_from / test_base_power_1p;
    double const s_ft = base_u_from * base_u_to / test_base_power_1p;
    double const s_tt = base_u_to * base_u_to / test_base_power_1p;
    DoubleComplex const half_low{0.0, -0.5e-8 * test_sn / test_base_power_1p / test_uk};
    auto const low_susceptance = [&half_low](ConnectionCase const& side, bool connected, size_t r, size_t c) {
        if (!connected || !side.phase_to_phase() || side.incidence[r] == 0.0 || side.incidence[c] == 0.0) {
            return DoubleComplex{};
        }
        return half_low;
    };
    auto const& c_f = e.from.incidence;
    std::array<double, 3> c_t = e.to.incidence;
    for (auto& value : c_t) {
        value *= e.polarity;
    }

    BranchCalcParam<asymmetric_t> param{};
    for (size_t r = 0; r != 3; ++r) {
        for (size_t c = 0; c != 3; ++c) {
            auto const row = static_cast<Idx>(r);
            auto const col = static_cast<Idx>(c);
            param.yff()(row, col) =
                e.y[0][0] * s_ff * c_f[r] * c_f[c] + low_susceptance(e.from, e.from_connected, r, c);
            param.yft()(row, col) = e.y[0][1] * s_ft * c_f[r] * c_t[c];
            param.ytf()(row, col) = e.y[1][0] * s_ft * c_t[r] * c_f[c];
            param.ytt()(row, col) = e.y[1][1] * s_tt * c_t[r] * c_t[c] + low_susceptance(e.to, e.to_connected, r, c);
        }
    }
    return param;
}

// exactly +0 (not -0, not nan)
bool is_positive_zero(DoubleComplex const& value) {
    return value == 0.0 && !std::signbit(value.real()) && !std::signbit(value.imag());
}

void check_block(ComplexTensor<asymmetric_t> const& actual, ComplexTensor<asymmetric_t> const& expected,
                 PhaseMask row_mask, PhaseMask col_mask) {
    CHECK(!actual.isNaN().any());
    for (Idx row = 0; row != 3; ++row) {
        for (Idx col = 0; col != 3; ++col) {
            CAPTURE(row);
            CAPTURE(col);
            if (!row_mask.has(row) || !col_mask.has(col)) {
                CHECK(is_positive_zero(actual(row, col)));
                continue;
            }
            CAPTURE(actual(row, col));
            CAPTURE(expected(row, col));
            CHECK(cabs(actual(row, col) - expected(row, col)) <
                  test_tolerance * std::max(1.0, cabs(expected(row, col))));
        }
    }
}

void check_param(BranchCalcParam<asymmetric_t> const& actual, ExpectedCase const& e) {
    BranchCalcParam<asymmetric_t> const expected = expected_param(e);
    check_block(actual.yff(), expected.yff(), e.from.mask, e.from.mask);
    check_block(actual.yft(), expected.yft(), e.from.mask, e.to.mask);
    check_block(actual.ytf(), expected.ytf(), e.to.mask, e.from.mask);
    check_block(actual.ytt(), expected.ytt(), e.to.mask, e.to.mask);
}

void check_all_zero(BranchCalcParam<asymmetric_t> const& param) {
    for (auto const& value : param.value) {
        for (Idx row = 0; row != 3; ++row) {
            for (Idx col = 0; col != 3; ++col) {
                CHECK(is_positive_zero(value(row, col)));
            }
        }
    }
}

bool params_equal(BranchCalcParam<asymmetric_t> const& x, BranchCalcParam<asymmetric_t> const& y) {
    return std::ranges::all_of(std::array<size_t, 4>{0, 1, 2, 3},
                               [&x, &y](size_t k) { return (x.value[k] == y.value[k]).all(); });
}

// equal modulo 2 pi
bool same_angle(double x, double y) { return cabs(std::polar(1.0, x) - std::polar(1.0, y)) < test_tolerance; }

void check_connection_pair(ConnectionCase const& from, ConnectionCase const& to) {
    SinglePhaseTransformer const transformer = make_transformer(test_input(from.connection, to.connection));
    CHECK(transformer.phase_incidence() ==
          BranchPhaseIncidence{.from = from.mask, .to = to.mask, .coupling = PhaseCoupling::winding});
    CHECK(same_angle(transformer.phase_shift(), to.angle() - from.angle()));
    check_param(transformer.calc_param<asymmetric_t>(), {.from = from, .to = to, .y = primitive({})});
}
} // namespace

TEST_CASE("Test single phase transformer - every connection pair, solidly earthed") {
    for (auto const& from : connection_cases()) {
        for (auto const& to : connection_cases()) {
            CAPTURE(phase_letters(from.mask));
            CAPTURE(phase_letters(to.mask));
            check_connection_pair(from, to);
        }
    }
}

TEST_CASE("Test single phase transformer - phase shift table") {
    using enum WindingConnection1P;
    auto const shift = [](WindingConnection1P from, WindingConnection1P to, IntS clock = 0) {
        SinglePhaseTransformerInput input = test_input(from, to);
        input.clock = clock;
        return make_transformer(input).phase_shift();
    };
    CHECK(same_angle(shift(ab, a), -pi / 6.0));
    CHECK(same_angle(shift(bc, a), pi / 2.0));
    CHECK(same_angle(shift(ac, a, 6), -5.0 * pi / 6.0)); // the design's ca winding
    CHECK(same_angle(shift(a, b), -2.0 * pi / 3.0));
    CHECK(same_angle(shift(a, a, 6), pi));
    CHECK(same_angle(shift(ab, ab), 0.0));
}

TEST_CASE("Test single phase transformer - earthing impedance") {
    ConnectionCase const phase_a = connection_case(WindingConnection1P::a);
    ConnectionCase const phase_ab = connection_case(WindingConnection1P::ab);
    DoubleComplex const z_from{3.0, 0.5};
    DoubleComplex const z_to{0.2, 0.05};
    Block2 const y = primitive({});

    // Sherman-Morrison: (Z_prim + diag(0, z))^-1 = Y - z / (1 + z Y22) Y[:, 1] Y[1, :]
    auto const with_to_electrode = [&y, &z_to]() {
        Block2 result = y;
        for (size_t i = 0; i != 2; ++i) {
            for (size_t j = 0; j != 2; ++j) {
                result[i][j] -= z_to / (1.0 + z_to * y[1][1]) * y[i][1] * y[1][j];
            }
        }
        return result;
    };

    SUBCASE("Phase-to-earth to side") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::a, WindingConnection1P::a);
        input.r_grounding_to = z_to.real();
        input.x_grounding_to = z_to.imag();
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a, .to = phase_a, .y = with_to_electrode()});
    }
    SUBCASE("Phase-to-earth from side") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::a, WindingConnection1P::a);
        input.r_grounding_from = z_from.real();
        input.x_grounding_from = z_from.imag();
        Block2 expected = y;
        for (size_t i = 0; i != 2; ++i) {
            for (size_t j = 0; j != 2; ++j) {
                expected[i][j] -= z_from / (1.0 + z_from * y[0][0]) * y[i][0] * y[0][j];
            }
        }
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a, .to = phase_a, .y = expected});
    }
    SUBCASE("Both sides, each with its own electrode") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::a, WindingConnection1P::a);
        input.r_grounding_from = z_from.real();
        input.x_grounding_from = z_from.imag();
        input.r_grounding_to = z_to.real();
        input.x_grounding_to = z_to.imag();
        // (Y_prim^-1 + diag(z_from, z_to))^-1, which needs a magnetising admittance
        Block2 const expected = inverse(add_diagonal(inverse(y), z_from, z_to));
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a, .to = phase_a, .y = expected});
    }
    SUBCASE("Phase-to-phase primary, earthed secondary") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::ab, WindingConnection1P::a);
        input.r_grounding_to = z_to.real();
        input.x_grounding_to = z_to.imag();
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_ab, .to = phase_a, .y = with_to_electrode()});
    }
    SUBCASE("Phase-to-phase primary, earthed secondary, clock 6") {
        // clock 6 negates the to-side incidence; the electrode is in series with the to winding whatever its
        // polarity, so the winding admittance is the same Sherman-Morrison form
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::ab, WindingConnection1P::a);
        input.r_grounding_to = z_to.real();
        input.x_grounding_to = z_to.imag();
        auto const param_0 = make_transformer(input).calc_param<asymmetric_t>();
        input.clock = 6;
        auto const param_6 = make_transformer(input).calc_param<asymmetric_t>();
        check_param(param_6, {.from = phase_ab, .to = phase_a, .y = with_to_electrode(), .polarity = -1.0});
        CHECK((param_6.yff() == param_0.yff()).all());
        CHECK((param_6.ytt() == param_0.ytt()).all());
        CHECK((param_6.yft() == -param_0.yft()).all());
        CHECK((param_6.ytf() == -param_0.ytf()).all());
    }
    SUBCASE("Only a resistance or only a reactance") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::a, WindingConnection1P::a);
        input.x_grounding_to = 0.4;
        Block2 expected = y;
        DoubleComplex const z{0.0, 0.4};
        for (size_t i = 0; i != 2; ++i) {
            for (size_t j = 0; j != 2; ++j) {
                expected[i][j] -= z / (1.0 + z * y[1][1]) * y[i][1] * y[1][j];
            }
        }
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a, .to = phase_a, .y = expected});
    }
    SUBCASE("Zero earthing impedance is solidly earthed") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::a, WindingConnection1P::b);
        input.r_grounding_from = 0.0;
        input.x_grounding_to = 0.0;
        CHECK(params_equal(
            make_transformer(input).calc_param<asymmetric_t>(),
            make_transformer(test_input(WindingConnection1P::a, WindingConnection1P::b)).calc_param<asymmetric_t>()));
    }
}

TEST_CASE("Test single phase transformer - low susceptance") {
    DoubleComplex const low_susceptance{0.0, -1e-8 * test_sn / test_base_power_1p / test_uk};
    // the differential part of a phase-to-phase winding has zero row sums, so a row sum of its block is the common
    // mode admittance y_low
    auto const row_sum = [](ComplexTensor<asymmetric_t> const& y_block, Idx row) {
        return y_block(row, 0) + y_block(row, 1) + y_block(row, 2);
    };

    SUBCASE("Present on every connected phase-to-phase side") {
        auto const param =
            make_transformer(test_input(WindingConnection1P::ab, WindingConnection1P::bc)).calc_param<asymmetric_t>();
        CHECK(cabs(row_sum(param.yff(), 0) - low_susceptance) < 1e-6 * cabs(low_susceptance));
        CHECK(cabs(row_sum(param.yff(), 1) - low_susceptance) < 1e-6 * cabs(low_susceptance));
        CHECK(cabs(row_sum(param.ytt(), 1) - low_susceptance) < 1e-6 * cabs(low_susceptance));
        CHECK(cabs(row_sum(param.ytt(), 2) - low_susceptance) < 1e-6 * cabs(low_susceptance));
    }
    SUBCASE("Present on the connected side only") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::ab, WindingConnection1P::ab);
        input.i0 = 0.0;
        input.p0 = 0.0;
        input.to_status = 0;
        auto const param = make_transformer(input).calc_param<asymmetric_t>();
        for (Idx row = 0; row != 2; ++row) {
            for (Idx col = 0; col != 2; ++col) {
                CHECK(cabs(param.yff()(row, col) - 0.5 * low_susceptance) < test_tolerance * cabs(low_susceptance));
            }
        }
        CHECK(is_positive_zero(param.yff()(2, 2)));
        check_all_zero({.value = {param.ytt(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
    }
    SUBCASE("Absent on a phase-to-earth side") {
        SinglePhaseTransformerInput input = test_input(WindingConnection1P::ab, WindingConnection1P::c);
        input.i0 = 0.0;
        input.p0 = 0.0;
        input.from_status = 0;
        // to side only, without a magnetising admittance: nothing at all
        check_all_zero(make_transformer(input).calc_param<asymmetric_t>());
    }
}

TEST_CASE("Test single phase transformer - clock") {
    using enum WindingConnection1P;
    auto const with_clock = [](IntS clock) {
        SinglePhaseTransformerInput input = test_input(ab, a);
        input.clock = clock;
        return make_transformer(input);
    };
    SinglePhaseTransformer const clock_0 = with_clock(0);
    SinglePhaseTransformer const clock_6 = with_clock(6);
    auto const param_0 = clock_0.calc_param<asymmetric_t>();
    auto const param_6 = clock_6.calc_param<asymmetric_t>();

    SUBCASE("Clock 6 reverses the polarity of the to winding") {
        check_param(param_6,
                    {.from = connection_case(ab), .to = connection_case(a), .y = primitive({}), .polarity = -1.0});
        CHECK((param_6.yff() == param_0.yff()).all());
        CHECK((param_6.ytt() == param_0.ytt()).all());
        CHECK((param_6.yft() == -param_0.yft()).all());
        CHECK((param_6.ytf() == -param_0.ytf()).all());
        CHECK(same_angle(clock_6.phase_shift(), clock_0.phase_shift() + pi));
    }
    SUBCASE("Reversed phase-to-phase and phase-to-earth windings") {
        SinglePhaseTransformerInput input = test_input(bc, ac);
        input.clock = 6;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = connection_case(bc), .to = connection_case(ac), .y = primitive({}), .polarity = -1.0});
        input = test_input(c, b);
        input.clock = 6;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = connection_case(c), .to = connection_case(b), .y = primitive({}), .polarity = -1.0});
    }
    SUBCASE("Periodic clock") {
        for (IntS const clock : {IntS{18}, IntS{-6}}) {
            CAPTURE(clock);
            SinglePhaseTransformer const transformer = with_clock(clock);
            CHECK(transformer.clock() == 6);
            CHECK(params_equal(transformer.calc_param<asymmetric_t>(), param_6));
        }
        CHECK(with_clock(12).clock() == 0);
        CHECK(params_equal(with_clock(-12).calc_param<asymmetric_t>(), param_0));
    }
    SUBCASE("Not provided is 0") {
        SinglePhaseTransformer const transformer = with_clock(na_IntS);
        CHECK(transformer.clock() == 0);
        CHECK(params_equal(transformer.calc_param<asymmetric_t>(), param_0));
    }
    SUBCASE("Other values are invalid") {
        for (IntS const clock : {IntS{1}, IntS{3}, IntS{5}, IntS{11}, IntS{-1}}) {
            CAPTURE(clock);
            CHECK_THROWS_AS(with_clock(clock), InvalidTransformerClock);
        }
    }
}

TEST_CASE("Test single phase transformer - status") {
    using enum WindingConnection1P;
    ConnectionCase const phase_a = connection_case(a);
    WindingData const nominal{};
    double const n = test_u1 / test_u2;

    SUBCASE("From side connected") {
        SinglePhaseTransformerInput input = test_input(a, a);
        input.to_status = 0;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a,
                     .to = phase_a,
                     .y = {{{open_end_shunt(nominal) / (n * n), 0.0}, {0.0, 0.0}}},
                     .to_connected = false});
    }
    SUBCASE("To side connected") {
        SinglePhaseTransformerInput input = test_input(a, a);
        input.from_status = 0;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a,
                     .to = phase_a,
                     .y = {{{0.0, 0.0}, {0.0, open_end_shunt(nominal)}}},
                     .from_connected = false});
    }
    SUBCASE("To side connected with an electrode") {
        SinglePhaseTransformerInput input = test_input(a, a);
        input.from_status = 0;
        input.r_grounding_to = 2.0;
        DoubleComplex const shunt = open_end_shunt(nominal);
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_a,
                     .to = phase_a,
                     .y = {{{0.0, 0.0}, {0.0, 1.0 / (1.0 / shunt + 2.0)}}},
                     .from_connected = false});
    }
    SUBCASE("Without magnetising admittance a single connected side is exactly zero") {
        for (auto const [from_status, to_status] : {std::array<IntS, 2>{1, 0}, std::array<IntS, 2>{0, 1}}) {
            SinglePhaseTransformerInput input = test_input(a, a);
            input.i0 = 0.0;
            input.p0 = 0.0;
            input.r_grounding_to = 5.0;
            input.from_status = from_status;
            input.to_status = to_status;
            check_all_zero(make_transformer(input).calc_param<asymmetric_t>());
        }
    }
    SUBCASE("Both sides open or not energized") {
        SinglePhaseTransformerInput input = test_input(ab, a);
        check_all_zero(make_transformer(input).calc_param<asymmetric_t>(false));
        input.from_status = 0;
        input.to_status = 0;
        check_all_zero(make_transformer(input).calc_param<asymmetric_t>());
    }
}

TEST_CASE("Test single phase transformer - taps") {
    using enum WindingConnection1P;
    ConnectionCase const phase_a = connection_case(a);
    ConnectionCase const phase_bc = connection_case(bc);
    auto const tapped_input = [](BranchSide tap_side, IntS tap_pos, double tap_size) {
        SinglePhaseTransformerInput input = test_input(bc, a);
        input.tap_side = tap_side;
        input.tap_pos = tap_pos;
        input.tap_min = -5;
        input.tap_max = 5;
        input.tap_nom = 0;
        input.tap_size = tap_size;
        return input;
    };

    SUBCASE("From side tap scales the ratio, with the same impedance on side 2") {
        SinglePhaseTransformer const transformer = make_transformer(tapped_input(BranchSide::from, 2, 200.0));
        CHECK(transformer.tap_pos() == 2);
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 + 400.0})});
    }
    SUBCASE("To side tap") {
        SinglePhaseTransformer const transformer = make_transformer(tapped_input(BranchSide::to, -3, 5.0));
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u2 = test_u2 - 15.0})});
    }
    SUBCASE("Reversed tap direction") {
        SinglePhaseTransformerInput input = tapped_input(BranchSide::from, 2, 200.0);
        input.tap_min = 5;
        input.tap_max = -5;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 - 400.0})});
    }
    SUBCASE("Tap dependent short circuit parameters") {
        SinglePhaseTransformerInput input = tapped_input(BranchSide::from, 2, 200.0);
        input.uk_min = 0.03;
        input.uk_max = 0.05;
        input.pk_min = 500.0;
        input.pk_max = 800.0;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 + 400.0,
                                     .uk = test_uk + 2.0 * (0.05 - test_uk) / 5.0,
                                     .pk = test_pk + 2.0 * (800.0 - test_pk) / 5.0})});
        input.tap_pos = -3;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 - 600.0,
                                     .uk = test_uk + 3.0 * (0.03 - test_uk) / 5.0,
                                     .pk = test_pk + 3.0 * (500.0 - test_pk) / 5.0})});
    }
    SUBCASE("Tap fields not provided give the nominal ratio") {
        SinglePhaseTransformer const transformer = make_transformer(test_input(bc, a));
        CHECK(transformer.tap_pos() == 0);
        CHECK(transformer.tap_nom() == 0);
        CHECK(transformer.tap_min() == 0);
        CHECK(transformer.tap_max() == 0);
        CHECK(transformer.tap_side() == BranchSide::from);
        check_param(transformer.calc_param<asymmetric_t>(), {.from = phase_bc, .to = phase_a, .y = primitive({})});
    }
    SUBCASE("Tap position not provided is the nominal tap; the tap position is limited") {
        SinglePhaseTransformerInput input = tapped_input(BranchSide::from, na_IntS, 200.0);
        input.tap_nom = 1;
        CHECK(make_transformer(input).tap_pos() == 1);
        input.tap_pos = 9;
        CHECK(make_transformer(input).tap_pos() == 5);
    }
    SUBCASE("Update and inverse") {
        SinglePhaseTransformer transformer = make_transformer(tapped_input(BranchSide::from, 2, 200.0));
        TransformerUpdate const update{.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = 4};
        TransformerUpdate const inverse = transformer.inverse(update);
        CHECK(inverse.tap_pos == 2);
        CHECK(inverse.from_status == na_IntS);

        auto const changed = transformer.update(update);
        CHECK_FALSE(changed.topo);
        CHECK(changed.param);
        CHECK(transformer.tap_pos() == 4);
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 + 800.0})});

        auto const unchanged = transformer.update(update);
        CHECK_FALSE(unchanged.topo);
        CHECK_FALSE(unchanged.param);

        transformer.update(inverse);
        CHECK(transformer.tap_pos() == 2);

        auto const status_changed =
            transformer.update({.id = 1, .from_status = 0, .to_status = na_IntS, .tap_pos = na_IntS});
        CHECK(status_changed.topo);
        CHECK(status_changed.param);
    }
}

TEST_CASE("Test single phase transformer - tap limits not provided") {
    using enum WindingConnection1P;
    ConnectionCase const phase_a = connection_case(a);
    ConnectionCase const phase_bc = connection_case(bc);
    auto const unlimited_input = [](IntS tap_pos) {
        SinglePhaseTransformerInput input = test_input(bc, a);
        input.tap_side = BranchSide::from;
        input.tap_pos = tap_pos;
        input.tap_nom = 0;
        input.tap_size = 200.0;
        return input;
    };
    auto const check_tap = [&](SinglePhaseTransformer const& transformer, IntS tap_pos) {
        CHECK(transformer.tap_pos() == tap_pos);
        CHECK(transformer.tap_min() == std::min(tap_pos, IntS{0}));
        CHECK(transformer.tap_max() == std::max(tap_pos, IntS{0}));
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 + tap_pos * 200.0})});
    };

    SUBCASE("The tap position is not limited") {
        check_tap(make_transformer(unlimited_input(3)), 3);
        check_tap(make_transformer(unlimited_input(-2)), -2);
        check_tap(make_transformer(unlimited_input(na_IntS)), 0);
    }
    SUBCASE("An update outside the limits at construction is not limited either") {
        SinglePhaseTransformer transformer = make_transformer(unlimited_input(3));
        TransformerUpdate const update{.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = 7};
        TransformerUpdate const inverse = transformer.inverse(update);
        CHECK(transformer.update(update).param);
        check_tap(transformer, 7);
        transformer.update({.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = -4});
        check_tap(transformer, -4);
        transformer.update(inverse);
        check_tap(transformer, 3);
    }
    SUBCASE("From the nominal tap, the tap direction is increasing") {
        SinglePhaseTransformer transformer = make_transformer(unlimited_input(na_IntS));
        transformer.update({.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = 2});
        check_tap(transformer, 2);
    }
    SUBCASE("A provided limit still limits the tap position") {
        SinglePhaseTransformerInput input = unlimited_input(3);
        input.tap_min = -1;
        SinglePhaseTransformer transformer = make_transformer(input);
        CHECK(transformer.tap_pos() == 3);
        CHECK(transformer.tap_min() == -1);
        CHECK(transformer.tap_max() == 3);
        transformer.update({.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = -4});
        CHECK(transformer.tap_pos() == -1);
        CHECK(transformer.tap_max() == 0);
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 - 200.0})});
    }
    SUBCASE("Tap dependent short circuit parameters only apply at a provided limit") {
        SinglePhaseTransformerInput input = unlimited_input(2);
        input.tap_min = -5;
        input.uk_min = 0.03;
        input.uk_max = 0.05;
        // no tap_max: uk_max does not apply
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 + 400.0})});
        input.tap_pos = -3;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 - 600.0, .uk = test_uk + 3.0 * (0.03 - test_uk) / 5.0})});
    }
}

TEST_CASE("Test single phase transformer - tap direction and limits") {
    using enum WindingConnection1P;
    ConnectionCase const phase_a = connection_case(a);
    ConnectionCase const phase_bc = connection_case(bc);
    constexpr IntS tap_nom = 5;
    constexpr double tap_size = 200.0;
    auto const limited_input = [](IntS tap_min, IntS tap_max, IntS tap_pos) {
        SinglePhaseTransformerInput input = test_input(bc, a);
        input.tap_side = BranchSide::from;
        input.tap_pos = tap_pos;
        input.tap_min = tap_min;
        input.tap_max = tap_max;
        input.tap_nom = tap_nom;
        input.tap_size = tap_size;
        return input;
    };
    // the tap position and the ratio: u1 + direction * (tap_pos - tap_nom) * tap_size
    auto const check_tap = [&](SinglePhaseTransformer const& transformer, IntS tap_pos, double direction) {
        CHECK(transformer.tap_pos() == tap_pos);
        check_param(transformer.calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 + direction * (tap_pos - tap_nom) * tap_size})});
    };
    auto const update_tap = [](SinglePhaseTransformer& transformer, IntS tap_pos) {
        transformer.update({.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = tap_pos});
    };

    struct LimitCase {
        IntS tap_min;
        IntS tap_max;
        double direction;
        // tap positions with the position they are limited to
        std::vector<std::array<IntS, 2>> positions;
    };
    auto const check_case = [&](LimitCase const& limit_case) {
        CAPTURE(limit_case.tap_min);
        CAPTURE(limit_case.tap_max);
        for (auto const [tap_pos, limited] : limit_case.positions) {
            CAPTURE(tap_pos);
            // at construction
            check_tap(make_transformer(limited_input(limit_case.tap_min, limit_case.tap_max, tap_pos)), limited,
                      limit_case.direction);
            // in an update, and back with the inverse
            SinglePhaseTransformer transformer =
                make_transformer(limited_input(limit_case.tap_min, limit_case.tap_max, tap_nom));
            TransformerUpdate const update{.id = 1, .from_status = na_IntS, .to_status = na_IntS, .tap_pos = tap_pos};
            TransformerUpdate const inverse = transformer.inverse(update);
            CHECK(inverse.tap_pos == tap_nom);
            transformer.update(update);
            check_tap(transformer, limited, limit_case.direction);
            transformer.update(inverse);
            check_tap(transformer, tap_nom, limit_case.direction);
        }
    };

    SUBCASE("Only tap_min, above tap_nom: reversed, limited from above, unbounded below") {
        check_case({.tap_min = 9,
                    .tap_max = na_IntS,
                    .direction = -1.0,
                    .positions = {{10, 9}, {12, 9}, {9, 9}, {7, 7}, {2, 2}, {-3, -3}}});
        SinglePhaseTransformer transformer = make_transformer(limited_input(9, na_IntS, 7));
        CHECK(transformer.tap_min() == 9);
        CHECK(transformer.tap_max() == tap_nom);
        update_tap(transformer, 2);
        CHECK(transformer.tap_min() == 9);
        CHECK(transformer.tap_max() == 2);
    }
    SUBCASE("Only tap_max, below tap_nom: reversed, limited from below, unbounded above") {
        check_case({.tap_min = na_IntS,
                    .tap_max = 1,
                    .direction = -1.0,
                    .positions = {{0, 1}, {-1, 1}, {1, 1}, {3, 3}, {8, 8}, {12, 12}}});
        SinglePhaseTransformer transformer = make_transformer(limited_input(na_IntS, 1, 3));
        CHECK(transformer.tap_min() == tap_nom);
        CHECK(transformer.tap_max() == 1);
        update_tap(transformer, 8);
        CHECK(transformer.tap_min() == 8);
        CHECK(transformer.tap_max() == 1);
    }
    SUBCASE("Only tap_min, at or below tap_nom: normal, limited from below, unbounded above") {
        check_case({.tap_min = 1,
                    .tap_max = na_IntS,
                    .direction = 1.0,
                    .positions = {{-2, 1}, {0, 1}, {3, 3}, {9, 9}, {12, 12}}});
        check_case({.tap_min = tap_nom, .tap_max = na_IntS, .direction = 1.0, .positions = {{3, tap_nom}, {9, 9}}});
        SinglePhaseTransformer transformer = make_transformer(limited_input(1, na_IntS, 3));
        CHECK(transformer.tap_min() == 1);
        CHECK(transformer.tap_max() == tap_nom);
        update_tap(transformer, 9);
        CHECK(transformer.tap_max() == 9);
    }
    SUBCASE("Only tap_max, at or above tap_nom: normal, limited from above, unbounded below") {
        check_case({.tap_min = na_IntS,
                    .tap_max = 9,
                    .direction = 1.0,
                    .positions = {{12, 9}, {10, 9}, {7, 7}, {0, 0}, {-3, -3}}});
        check_case({.tap_min = na_IntS, .tap_max = tap_nom, .direction = 1.0, .positions = {{7, tap_nom}, {0, 0}}});
        SinglePhaseTransformer transformer = make_transformer(limited_input(na_IntS, 9, 7));
        CHECK(transformer.tap_min() == tap_nom);
        CHECK(transformer.tap_max() == 9);
        update_tap(transformer, 0);
        CHECK(transformer.tap_min() == 0);
    }
    SUBCASE("Both limits, normal") {
        check_case({.tap_min = 1,
                    .tap_max = 9,
                    .direction = 1.0,
                    .positions = {{12, 9}, {9, 9}, {7, 7}, {3, 3}, {1, 1}, {-1, 1}}});
    }
    SUBCASE("Both limits, reversed") {
        check_case({.tap_min = 9,
                    .tap_max = 1,
                    .direction = -1.0,
                    .positions = {{12, 9}, {9, 9}, {7, 7}, {3, 3}, {1, 1}, {-1, 1}}});
        SinglePhaseTransformer const transformer = make_transformer(limited_input(9, 1, 7));
        CHECK(transformer.tap_min() == 9);
        CHECK(transformer.tap_max() == 1);
    }
    SUBCASE("No limits: normal and unbounded") {
        check_case({.tap_min = na_IntS,
                    .tap_max = na_IntS,
                    .direction = 1.0,
                    .positions = {{12, 12}, {7, 7}, {3, 3}, {-1, -1}, {-20, -20}}});
    }
    SUBCASE("Tap dependent short circuit parameters with one limit on the reversed side") {
        // only tap_min above tap_nom: uk_min applies between tap_nom and tap_min, uk_max does not apply
        SinglePhaseTransformerInput input = limited_input(9, na_IntS, 7);
        input.uk_min = 0.03;
        input.uk_max = 0.05;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 - 2.0 * tap_size, .uk = test_uk + 2.0 * (0.03 - test_uk) / 4.0})});
        input.tap_pos = 3;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 + 2.0 * tap_size})});
        // only tap_max below tap_nom: uk_max applies between tap_max and tap_nom, uk_min does not apply
        input = limited_input(na_IntS, 1, 3);
        input.uk_min = 0.03;
        input.uk_max = 0.05;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc,
                     .to = phase_a,
                     .y = primitive({.u1 = test_u1 + 2.0 * tap_size, .uk = test_uk + 2.0 * (0.05 - test_uk) / 4.0})});
        input.tap_pos = 7;
        check_param(make_transformer(input).calc_param<asymmetric_t>(),
                    {.from = phase_bc, .to = phase_a, .y = primitive({.u1 = test_u1 - 2.0 * tap_size})});
    }
}

TEST_CASE("Test single phase transformer - invalid input") {
    using enum WindingConnection1P;

    SUBCASE("Connection not provided or unknown") {
        SinglePhaseTransformerInput input = test_input(a, a);
        input.from_connection = SinglePhaseTransformerInput{}.from_connection;
        CHECK_THROWS_AS(make_transformer(input), MissingCaseForEnumError);
        input = test_input(a, a);
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
        input.to_connection = static_cast<WindingConnection1P>(IntS{0});
        CHECK_THROWS_AS(make_transformer(input), MissingCaseForEnumError);
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
        input.to_connection = static_cast<WindingConnection1P>(IntS{9});
        CHECK_THROWS_AS(make_transformer(input), MissingCaseForEnumError);
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
        input.to_connection = static_cast<WindingConnection1P>(IntS{8});
        CHECK_THROWS_AS(make_transformer(input), MissingCaseForEnumError);
        input = test_input(a, a);
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
        input.from_connection = static_cast<WindingConnection1P>(IntS{8});
        CHECK_THROWS_AS(make_transformer(input), MissingCaseForEnumError);
    }
    SUBCASE("Earthing impedance on a phase-to-phase winding") {
        SinglePhaseTransformerInput input = test_input(ab, a);
        input.r_grounding_from = 1.0;
        CHECK_THROWS_WITH_AS(make_transformer(input),
                             "Invalid phase connection for single_phase_transformer 1 with phases (ab): "
                             "r_grounding_from refers to a phase-to-phase winding and must be left empty\n",
                             InvalidPhaseConnection);
        input = test_input(a, bc);
        input.x_grounding_to = 0.0;
        CHECK_THROWS_WITH_AS(make_transformer(input),
                             "Invalid phase connection for single_phase_transformer 1 with phases (bc): "
                             "x_grounding_to refers to a phase-to-phase winding and must be left empty\n",
                             InvalidPhaseConnection);
    }
    SUBCASE("Symmetric parameters are not supported") {
        SinglePhaseTransformer const transformer = make_transformer(test_input(a, a));
        CHECK_THROWS_WITH_AS(transformer.calc_param<symmetric_t>(),
                             "Symmetric calculation is not supported in a network with one- or two-phase nodes.\n",
                             ReducedPhaseNotSupported);
    }
}

TEST_CASE("Test single phase transformer - output") {
    using enum WindingConnection1P;
    SinglePhaseTransformer const transformer = make_transformer(test_input(ab, a));
    CHECK(transformer.math_model_type() == ComponentType::branch);
    CHECK(transformer.base_i_from() == doctest::Approx(1.0e6 / test_u_rated_from / sqrt3));
    CHECK(transformer.base_i_to() == doctest::Approx(1.0e6 / test_u_rated_to / sqrt3));

    SUBCASE("Asymmetric: loading from the complex sum of the per-phase powers of each side") {
        DoubleComplex const s_a{0.06, 0.02};
        DoubleComplex const s_b{0.05, -0.01};
        BranchSolverOutput<asymmetric_t> const solver_output{
            .s_f = ComplexValue<asymmetric_t>{s_a, s_b, 0.0},
            .s_t = ComplexValue<asymmetric_t>{DoubleComplex{-0.1, -0.005}, 0.0, 0.0},
            .i_f = ComplexValue<asymmetric_t>{0.1, 0.1, 0.0},
            .i_t = ComplexValue<asymmetric_t>{2.0, 0.0, 0.0}};
        auto const output = transformer.get_output<asymmetric_t>(solver_output);
        CHECK(output.loading == doctest::Approx(cabs(s_a + s_b) * test_base_power_1p / test_sn));
        CHECK(output.p_from(0) == doctest::Approx(0.06 * test_base_power_1p));
        CHECK(output.p_from(2) == 0.0);
        CHECK(output.s_from(1) == doctest::Approx(cabs(s_b) * test_base_power_1p));
        CHECK(output.i_to(0) == doctest::Approx(2.0 * transformer.base_i_to()));
        // the sum of the per-phase apparent powers would be larger
        CHECK(output.loading < (cabs(s_a) + cabs(s_b)) * test_base_power_1p / test_sn);
    }
    SUBCASE("Symmetric") {
        BranchSolverOutput<symmetric_t> const solver_output{
            .s_f = DoubleComplex{0.01, 0.002}, .s_t = DoubleComplex{-0.02, 0.0}, .i_f = 0.01, .i_t = 0.02};
        auto const output = transformer.get_output<symmetric_t>(solver_output);
        CHECK(output.loading == doctest::Approx(0.02 * 1.0e6 / test_sn));
    }
}

TEST_CASE("Test single phase transformer - per unit with phase-to-phase primary") {
    // 22 kV ab -> 12.7 kV a, both nodes at u_rated 22 kV: nominal ratio in per unit, with a 30 degree shift
    double const u_rated = 22.0e3;
    SinglePhaseTransformerInput input = test_input(WindingConnection1P::ab, WindingConnection1P::a);
    input.u1 = u_rated;
    input.u2 = u_rated / sqrt3;
    input.i0 = 0.0;
    input.p0 = 0.0;
    SinglePhaseTransformer const transformer{input, u_rated, u_rated};
    auto const param = transformer.calc_param<asymmetric_t>();

    // series admittance in per unit of the to side
    WindingData const data{.u1 = input.u1, .u2 = input.u2, .i0 = 0.0, .p0 = 0.0};
    DoubleComplex const y_series = series_admittance(data) * (u_rated * u_rated / 1.0e6);
    CHECK(cabs(param.ytt()(0, 0) - y_series) < test_tolerance * cabs(y_series));
    CHECK(cabs(param.yft()(0, 0) + y_series / sqrt3) < test_tolerance * cabs(y_series));
    CHECK(cabs(param.yft()(1, 0) - y_series / sqrt3) < test_tolerance * cabs(y_series));

    // no load: balanced from voltages and U_a,to = (U_a - U_b) / sqrt3 = 1 at +30 degrees give no to current
    ComplexValue<asymmetric_t> const u_from{DoubleComplex{1.0, 0.0}, std::polar(1.0, -2.0 * pi / 3.0),
                                            std::polar(1.0, 2.0 * pi / 3.0)};
    ComplexValue<asymmetric_t> const u_to{std::polar(1.0, pi / 6.0), 0.0, 0.0};
    ComplexValue<asymmetric_t> const i_to = dot(param.ytf(), u_from) + dot(param.ytt(), u_to);
    ComplexValue<asymmetric_t> const i_from = dot(param.yff(), u_from) + dot(param.yft(), u_to);
    CHECK(cabs(i_to(0)) < test_tolerance * cabs(y_series));
    // only the common-mode low susceptance remains at the from side
    CHECK(cabs(i_from(0)) < 1e-6);
    CHECK(same_angle(transformer.phase_shift(), -pi / 6.0));
}

namespace {
// centre-tapped (split-phase) secondary: the expected values below are computed from the closed forms of the M6
// design in per unit on sn (the primitive admittance of the primary and the two halves, entry by entry, then
// C^T Y_prim C over the from slots, the to slots and the centre tap, with a scalar Schur complement for an
// electrode), with std::array scalars, independently of the N-winding product, the nodal stamp and the Kron
// elimination of the code under test
constexpr double ct_u2 = 480.0;                 // full winding, 240-0-240 V
constexpr double ct_u_rated_to = 240.0 * sqrt3; // each leg at 240 V to the earthed centre tap

using Block3 = std::array<std::array<DoubleComplex, 3>, 3>;

SinglePhaseTransformer make_ct_transformer(SinglePhaseTransformerInput const& input) {
    return SinglePhaseTransformer{input, test_u_rated_from, ct_u_rated_to};
}

SinglePhaseTransformerInput ct_input(WindingConnection1P from_connection) {
    SinglePhaseTransformerInput input = test_input(from_connection, WindingConnection1P::ab_center_tapped);
    input.u2 = ct_u2;
    return input;
}

struct CenterTappedData {
    double u1{test_u1};
    double u2{ct_u2};
    // tap-adjusted short circuit parameters, for the full winding z_f
    double uk{test_uk};
    double pk{test_pk};
    // nominal short circuit parameters, for the default half-winding impedance and the half branch z_b
    double uk_nominal{test_uk};
    double pk_nominal{test_pk};
    double uk_half{nan};
    double pk_half{nan};
    double i0{test_i0};
    double p0{test_p0};
};

// short circuit impedance in per unit on sn
DoubleComplex short_circuit_pu(double uk, double pk) {
    double const r = pk / test_sn;
    return {r, std::sqrt(uk * uk - r * r)};
}

struct CenterTappedImpedances {
    DoubleComplex z_full;           // z_f: primary to both halves in series
    DoubleComplex z_half;           // z_hs: primary to one half, the other half open
    DoubleComplex z_primary_branch; // z_0
    DoubleComplex z_half_branch;    // z_b
    DoubleComplex y_m;              // magnetising admittance
};

CenterTappedImpedances ct_impedances(CenterTappedData const& data) {
    DoubleComplex const z_full = short_circuit_pu(data.uk, data.pk);
    DoubleComplex const z_nominal = short_circuit_pu(data.uk_nominal, data.pk_nominal);
    // Kersting: primary 0.5 r + j0.8 x, each half r + j0.4 x
    DoubleComplex const kersting{1.5 * z_nominal.real(), 1.2 * z_nominal.imag()};
    DoubleComplex z_half = kersting;
    if (!is_nan(data.uk_half) || !is_nan(data.pk_half)) {
        // on the pair rating sn / 2 and the half voltage: r = 4 pk_half / sn, |z| = 2 uk_half
        double const r_half = is_nan(data.pk_half) ? kersting.real() : 4.0 * data.pk_half / test_sn;
        double const z_half_abs = is_nan(data.uk_half) ? cabs(kersting) : 2.0 * data.uk_half;
        z_half = {r_half, std::sqrt(z_half_abs * z_half_abs - r_half * r_half)};
    }
    DoubleComplex const z_half_branch = 2.0 * (z_half - z_nominal);
    double const g_m = data.p0 / test_sn;
    return {.z_full = z_full,
            .z_half = z_half,
            .z_primary_branch = z_full - 0.5 * z_half_branch,
            .z_half_branch = z_half_branch,
            .y_m = {g_m, -std::sqrt(data.i0 * data.i0 - g_m * g_m)}};
}

// primitive admittance of the primary and the two halves in siemens: Y_ij(pu) sn / (V_i V_j), with
// Y_00 = 1/z_f + y_m/2, Y_01 = Y_02 = -1/(2 z_f), Y_11 = Y_22 = (z_0 + z_b)/(2 z_b z_f) + y_m/4,
// Y_12 = -z_0/(2 z_b z_f)
Block3 primitive_ct(CenterTappedData const& data) {
    auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances(data);
    DoubleComplex const self_primary = 1.0 / z_f + 0.5 * y_m;
    DoubleComplex const mutual = -1.0 / (2.0 * z_f);
    DoubleComplex const self_half = (z_0 + z_b) / (2.0 * z_b * z_f) + 0.25 * y_m;
    DoubleComplex const between_halves = -z_0 / (2.0 * z_b * z_f);
    Block3 const y_pu{
        {{self_primary, mutual, mutual}, {mutual, self_half, between_halves}, {mutual, between_halves, self_half}}};
    std::array<double, 3> const voltage{data.u1, 0.5 * data.u2, 0.5 * data.u2};
    Block3 y{};
    for (size_t i = 0; i != 3; ++i) {
        for (size_t j = 0; j != 3; ++j) {
            y[i][j] = y_pu[i][j] * test_sn / (voltage[i] * voltage[j]);
        }
    }
    return y;
}

struct ExpectedCenterTapped {
    ConnectionCase from;
    Block3 y{};           // admittance of the primary and the two halves in siemens
    double polarity{1.0}; // -1 for clock 6
    DoubleComplex z_grounding{};
    bool from_connected{true};
};

// C^T Y C over the from slots (0-2), the to slots (3-5) and the centre tap (6): half 1 is U_a - U_N, half 2 is
// U_N - U_b; a solid centre tap is earth (its row and column are dropped), an electrode is a scalar Schur complement;
// scaled to per unit, plus the low susceptance y_low / 2 on a connected phase-to-phase from side
BranchCalcParam<asymmetric_t> expected_param_ct(ExpectedCenterTapped const& e) {
    constexpr size_t n_node = 7;
    constexpr size_t centre_tap = 6;
    std::array<std::array<double, n_node>, 3> incidence{};
    for (size_t k = 0; k != 3; ++k) {
        incidence[0][k] = e.from.incidence[k];
    }
    incidence[1][3] = e.polarity;
    incidence[1][centre_tap] = -e.polarity;
    incidence[2][centre_tap] = e.polarity;
    incidence[2][4] = -e.polarity;

    std::array<std::array<DoubleComplex, n_node>, n_node> y_node{};
    for (size_t w = 0; w != 3; ++w) {
        for (size_t v = 0; v != 3; ++v) {
            for (size_t r = 0; r != n_node; ++r) {
                for (size_t c = 0; c != n_node; ++c) {
                    y_node[r][c] += e.y[w][v] * incidence[w][r] * incidence[v][c];
                }
            }
        }
    }
    if (e.z_grounding != DoubleComplex{}) {
        y_node[centre_tap][centre_tap] += 1.0 / e.z_grounding;
        for (size_t r = 0; r != centre_tap; ++r) {
            for (size_t c = 0; c != centre_tap; ++c) {
                y_node[r][c] -= y_node[r][centre_tap] * y_node[centre_tap][c] / y_node[centre_tap][centre_tap];
            }
        }
    }

    double const base_u_from = test_u_rated_from / sqrt3;
    double const base_u_to = ct_u_rated_to / sqrt3;
    double const s_ff = base_u_from * base_u_from / test_base_power_1p;
    double const s_ft = base_u_from * base_u_to / test_base_power_1p;
    double const s_tt = base_u_to * base_u_to / test_base_power_1p;
    DoubleComplex const half_low{0.0, -0.5e-8 * test_sn / test_base_power_1p / test_uk};
    bool const low_susceptance = e.from_connected && e.from.phase_to_phase();

    BranchCalcParam<asymmetric_t> param{};
    for (size_t r = 0; r != 3; ++r) {
        for (size_t c = 0; c != 3; ++c) {
            auto const row = static_cast<Idx>(r);
            auto const col = static_cast<Idx>(c);
            bool const both_present = e.from.incidence[r] != 0.0 && e.from.incidence[c] != 0.0;
            param.yff()(row, col) =
                y_node[r][c] * s_ff + (low_susceptance && both_present ? half_low : DoubleComplex{});
            param.yft()(row, col) = y_node[r][c + 3] * s_ft;
            param.ytf()(row, col) = y_node[r + 3][c] * s_ft;
            param.ytt()(row, col) = y_node[r + 3][c + 3] * s_tt;
        }
    }
    return param;
}

void check_param_ct(BranchCalcParam<asymmetric_t> const& actual, ExpectedCenterTapped const& e) {
    BranchCalcParam<asymmetric_t> const expected = expected_param_ct(e);
    check_block(actual.yff(), expected.yff(), e.from.mask, e.from.mask);
    check_block(actual.yft(), expected.yft(), e.from.mask, phase_mask_ab);
    check_block(actual.ytf(), expected.ytf(), phase_mask_ab, e.from.mask);
    check_block(actual.ytt(), expected.ytt(), phase_mask_ab, phase_mask_ab);
}

bool params_close(BranchCalcParam<asymmetric_t> const& x, BranchCalcParam<asymmetric_t> const& y) {
    for (size_t k = 0; k != 4; ++k) {
        for (Idx row = 0; row != 3; ++row) {
            for (Idx col = 0; col != 3; ++col) {
                if (cabs(x.value[k](row, col) - y.value[k](row, col)) >
                    test_tolerance * std::max(1.0, cabs(y.value[k](row, col)))) {
                    return false;
                }
            }
        }
    }
    return true;
}

// a 2x2 block of the to slots a and b
Block2 to_block(ComplexTensor<asymmetric_t> const& y) { return {{{y(0, 0), y(0, 1)}, {y(1, 0), y(1, 1)}}}; }

// impedance in ohm seen at slot a with the from side shorted and slot b open, from ytt in per unit
DoubleComplex leg_impedance(BranchCalcParam<asymmetric_t> const& param) {
    Block2 const y = to_block(param.ytt());
    double const base_u_to = ct_u_rated_to / sqrt3;
    return 1.0 / (y[0][0] - y[0][1] * y[1][0] / y[1][1]) * base_u_to * base_u_to / test_base_power_1p;
}

// impedance in ohm between slots a and b with the from side shorted: [1, -1] ytt^-1 [1, -1]^T
DoubleComplex a_to_b_impedance(BranchCalcParam<asymmetric_t> const& param) {
    Block2 const z = inverse(to_block(param.ytt()));
    double const base_u_to = ct_u_rated_to / sqrt3;
    return (z[0][0] - z[0][1] - z[1][0] + z[1][1]) * base_u_to * base_u_to / test_base_power_1p;
}

bool close_to(DoubleComplex const& x, DoubleComplex const& y) { return cabs(x - y) < test_tolerance * cabs(y); }
} // namespace

TEST_CASE("Test single phase transformer - centre-tapped secondary, every primary, solidly earthed") {
    for (auto const& from : connection_cases()) {
        CAPTURE(phase_letters(from.mask));
        SinglePhaseTransformer const transformer = make_ct_transformer(ct_input(from.connection));
        CHECK(transformer.phase_incidence() == BranchPhaseIncidence{.from = from.mask,
                                                                    .to = phase_mask_ab,
                                                                    .coupling = PhaseCoupling::winding,
                                                                    .to_pattern = SlotPattern::split_phase});
        // the reference slot a is the positive terminal of half 1, at angle 0
        CHECK(same_angle(transformer.phase_shift(), -from.angle()));
        // the expected param has the low susceptance on a phase-to-phase from side only, never on the earthed
        // centre-tapped side
        auto const param = transformer.calc_param<asymmetric_t>();
        check_param_ct(param, {.from = from, .y = primitive_ct({})});
        if (from.phase_to_phase()) {
            // the differential part of the from block has zero row sums, so a row sum is the common mode y_low
            DoubleComplex const low_susceptance{0.0, -1e-8 * test_sn / test_base_power_1p / test_uk};
            for (Idx row = 0; row != 3; ++row) {
                if (from.mask.has(row)) {
                    CAPTURE(row);
                    DoubleComplex const row_sum = param.yff()(row, 0) + param.yff()(row, 1) + param.yff()(row, 2);
                    CHECK(cabs(row_sum - low_susceptance) < 1e-6 * cabs(low_susceptance));
                }
            }
        }
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, half-winding impedance") {
    using enum WindingConnection1P;
    DoubleComplex const z_nominal = short_circuit_pu(test_uk, test_pk);
    double const uk_half_kersting = 0.5 * cabs(DoubleComplex{1.5 * z_nominal.real(), 1.2 * z_nominal.imag()});
    double const pk_half_kersting = 0.375 * test_pk;

    SUBCASE("The Kersting defaults equal the explicit overrides") {
        SinglePhaseTransformerInput input = ct_input(a);
        auto const defaults = make_ct_transformer(input).calc_param<asymmetric_t>();
        input.uk_half = uk_half_kersting;
        input.pk_half = pk_half_kersting;
        auto const overrides = make_ct_transformer(input).calc_param<asymmetric_t>();
        CHECK(params_close(overrides, defaults));
        check_param_ct(defaults, {.from = connection_case(a), .y = primitive_ct({})});
        check_param_ct(overrides, {.from = connection_case(a),
                                   .y = primitive_ct({.uk_half = uk_half_kersting, .pk_half = pk_half_kersting})});
    }
    SUBCASE("Overrides, each on its own and together") {
        constexpr double uk_half = 0.03;
        constexpr double pk_half = 200.0;
        SinglePhaseTransformerInput input = ct_input(bc);
        input.uk_half = uk_half;
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = connection_case(bc), .y = primitive_ct({.uk_half = uk_half})});
        input.uk_half = nan;
        input.pk_half = pk_half;
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = connection_case(bc), .y = primitive_ct({.pk_half = pk_half})});
        input.uk_half = uk_half;
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = connection_case(bc), .y = primitive_ct({.uk_half = uk_half, .pk_half = pk_half})});
    }
    SUBCASE("Leg and a-to-b impedances with the primary shorted") {
        constexpr double uk_half = 0.03;
        constexpr double pk_half = 200.0;
        SinglePhaseTransformerInput input = ct_input(a);
        input.uk_half = uk_half;
        input.pk_half = pk_half;
        input.i0 = 0.0;
        input.p0 = 0.0;
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        double const r_half = 4.0 * pk_half / test_sn;
        DoubleComplex const z_half{r_half, std::sqrt(4.0 * uk_half * uk_half - r_half * r_half)};
        CHECK(close_to(leg_impedance(param), z_half * (0.5 * ct_u2) * (0.5 * ct_u2) / test_sn));
        CHECK(close_to(a_to_b_impedance(param), z_nominal * ct_u2 * ct_u2 / test_sn));
        // the Kersting defaults: 1.5 r + j1.2 x
        input.uk_half = nan;
        input.pk_half = nan;
        auto const param_default = make_ct_transformer(input).calc_param<asymmetric_t>();
        CHECK(close_to(leg_impedance(param_default), DoubleComplex{1.5 * z_nominal.real(), 1.2 * z_nominal.imag()} *
                                                         (0.5 * ct_u2) * (0.5 * ct_u2) / test_sn));
        CHECK(close_to(a_to_b_impedance(param_default), z_nominal * ct_u2 * ct_u2 / test_sn));
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary equals a phase-to-phase secondary in "
          "differential mode") {
    using enum WindingConnection1P;
    // a -> ct and a -> ab, both with u2 = 480 V at the same node: with U_to = (v, -v, 0) both halves see v
    SinglePhaseTransformer const center_tapped = make_ct_transformer(ct_input(a));
    SinglePhaseTransformerInput input = test_input(a, ab);
    input.u2 = ct_u2;
    SinglePhaseTransformer const phase_to_phase{input, test_u_rated_from, ct_u_rated_to};
    auto const param_ct = center_tapped.calc_param<asymmetric_t>();
    auto const param_ab = phase_to_phase.calc_param<asymmetric_t>();

    ComplexValue<asymmetric_t> const u_from{DoubleComplex{0.98, 0.03}, DoubleComplex{-0.4, -0.9},
                                            DoubleComplex{-0.6, 0.8}};
    DoubleComplex const v{0.97, -0.05};
    ComplexValue<asymmetric_t> const u_to{v, -v, 0.0};
    auto const currents = [&u_from, &u_to](BranchCalcParam<asymmetric_t> const& param) {
        return std::array<ComplexValue<asymmetric_t>, 2>{dot(param.yff(), u_from) + dot(param.yft(), u_to),
                                                         dot(param.ytf(), u_from) + dot(param.ytt(), u_to)};
    };
    auto const [i_from_ct, i_to_ct] = currents(param_ct);
    auto const [i_from_ab, i_to_ab] = currents(param_ab);
    double const scale = std::max(cabs(i_from_ab(0)), cabs(i_to_ab(0)));
    for (Idx k = 0; k != 3; ++k) {
        CAPTURE(k);
        CHECK(cabs(i_from_ct(k) - i_from_ab(k)) < test_tolerance * scale);
        CHECK(cabs(i_to_ct(k) - i_to_ab(k)) < test_tolerance * scale);
    }
    CHECK(cabs(i_to_ct(0)) > 1e-3 * scale);
    CHECK(same_angle(center_tapped.phase_shift(), 0.0));
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, earthing, clock and magnetising") {
    using enum WindingConnection1P;
    DoubleComplex const z_grounding{0.5, 0.05};

    SUBCASE("Electrode on the centre tap") {
        SinglePhaseTransformerInput input = ct_input(ab);
        input.r_grounding_to = z_grounding.real();
        input.x_grounding_to = z_grounding.imag();
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        check_param_ct(param, {.from = connection_case(ab), .y = primitive_ct({}), .z_grounding = z_grounding});

        // differential excitation: no current in the centre tap, so the currents equal those of the solid case
        auto const solid = make_ct_transformer(ct_input(ab)).calc_param<asymmetric_t>();
        ComplexValue<asymmetric_t> const u_from{DoubleComplex{1.0, 0.02}, DoubleComplex{-0.5, -0.85}, 0.0};
        DoubleComplex const v{0.96, 0.31};
        ComplexValue<asymmetric_t> const u_to{v, -v, 0.0};
        ComplexValue<asymmetric_t> const i_to = dot(param.ytf(), u_from) + dot(param.ytt(), u_to);
        ComplexValue<asymmetric_t> const i_to_solid = dot(solid.ytf(), u_from) + dot(solid.ytt(), u_to);
        ComplexValue<asymmetric_t> const i_from = dot(param.yff(), u_from) + dot(param.yft(), u_to);
        ComplexValue<asymmetric_t> const i_from_solid = dot(solid.yff(), u_from) + dot(solid.yft(), u_to);
        double const scale = std::max(cabs(i_to_solid(0)), cabs(i_from_solid(0)));
        for (Idx k = 0; k != 3; ++k) {
            CAPTURE(k);
            CHECK(cabs(i_to(k) - i_to_solid(k)) < test_tolerance * scale);
            CHECK(cabs(i_from(k) - i_from_solid(k)) < test_tolerance * scale);
        }

        // the leg impedance includes the electrode
        input.i0 = 0.0;
        input.p0 = 0.0;
        auto const no_shunt = make_ct_transformer(input).calc_param<asymmetric_t>();
        DoubleComplex const z_nominal = short_circuit_pu(test_uk, test_pk);
        DoubleComplex const z_half{1.5 * z_nominal.real(), 1.2 * z_nominal.imag()};
        CHECK(close_to(leg_impedance(no_shunt), z_half * (0.5 * ct_u2) * (0.5 * ct_u2) / test_sn + z_grounding));
        // a to b does not pass the centre tap
        CHECK(close_to(a_to_b_impedance(no_shunt), z_nominal * ct_u2 * ct_u2 / test_sn));
    }
    SUBCASE("Only a reactance, and zero is solid") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.x_grounding_to = 0.3;
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = connection_case(a), .y = primitive_ct({}), .z_grounding = {0.0, 0.3}});
        input.x_grounding_to = 0.0;
        input.r_grounding_to = 0.0;
        CHECK(params_equal(make_ct_transformer(input).calc_param<asymmetric_t>(),
                           make_ct_transformer(ct_input(a)).calc_param<asymmetric_t>()));
    }
    SUBCASE("Clock 6 reverses both halves") {
        for (auto const& from : connection_cases()) {
            CAPTURE(phase_letters(from.mask));
            SinglePhaseTransformerInput input = ct_input(from.connection);
            auto const param_0 = make_ct_transformer(input).calc_param<asymmetric_t>();
            input.clock = 6;
            SinglePhaseTransformer const transformer = make_ct_transformer(input);
            auto const param_6 = transformer.calc_param<asymmetric_t>();
            check_param_ct(param_6, {.from = from, .y = primitive_ct({}), .polarity = -1.0});
            CHECK((param_6.yff() == param_0.yff()).all());
            CHECK((param_6.ytt() == param_0.ytt()).all());
            CHECK((param_6.yft() == -param_0.yft()).all());
            CHECK((param_6.ytf() == -param_0.ytf()).all());
            CHECK(same_angle(transformer.phase_shift(), -from.angle() + pi));
        }
        // with an electrode on the centre tap
        SinglePhaseTransformerInput input = ct_input(c);
        input.clock = 6;
        input.r_grounding_to = z_grounding.real();
        input.x_grounding_to = z_grounding.imag();
        check_param_ct(
            make_ct_transformer(input).calc_param<asymmetric_t>(),
            {.from = connection_case(c), .y = primitive_ct({}), .polarity = -1.0, .z_grounding = z_grounding});
    }
    SUBCASE("Magnetising admittance diag(y_m / 2, y_m / 4, y_m / 4)") {
        SinglePhaseTransformerInput input = ct_input(a);
        auto const with_shunt = make_ct_transformer(input).calc_param<asymmetric_t>();
        input.i0 = 0.0;
        input.p0 = 0.0;
        auto const without_shunt = make_ct_transformer(input).calc_param<asymmetric_t>();
        double const g_m = test_p0 / test_sn;
        DoubleComplex const y_m{g_m, -std::sqrt(test_i0 * test_i0 - g_m * g_m)};
        double const base_u_from = test_u_rated_from / sqrt3;
        double const base_u_to = ct_u_rated_to / sqrt3;
        DoubleComplex const from_shunt =
            0.5 * y_m * test_sn / (test_u1 * test_u1) * base_u_from * base_u_from / test_base_power_1p;
        DoubleComplex const half_shunt =
            0.25 * y_m * test_sn / (0.25 * ct_u2 * ct_u2) * base_u_to * base_u_to / test_base_power_1p;
        CHECK(close_to(with_shunt.yff()(0, 0) - without_shunt.yff()(0, 0), from_shunt));
        CHECK(close_to(with_shunt.ytt()(0, 0) - without_shunt.ytt()(0, 0), half_shunt));
        CHECK(close_to(with_shunt.ytt()(1, 1) - without_shunt.ytt()(1, 1), half_shunt));
        CHECK(cabs(with_shunt.ytt()(0, 1) - without_shunt.ytt()(0, 1)) < test_tolerance * cabs(half_shunt));
        CHECK(cabs(with_shunt.yft()(0, 0) - without_shunt.yft()(0, 0)) < test_tolerance * cabs(half_shunt));
        CHECK(cabs(with_shunt.yft()(0, 1) - without_shunt.yft()(0, 1)) < test_tolerance * cabs(half_shunt));
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, status") {
    using enum WindingConnection1P;
    ConnectionCase const phase_a = connection_case(a);
    CenterTappedData const no_shunt{.i0 = 0.0, .p0 = 0.0};

    SUBCASE("From side open: the halves stay coupled as a balancer") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.i0 = 0.0;
        input.p0 = 0.0;
        input.from_status = 0;
        auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances(no_shunt);
        // 1 / (2 z_b) [[1, -1], [-1, 1]] in half-winding per unit
        DoubleComplex const balancer = 1.0 / (2.0 * z_b) * test_sn / (0.25 * ct_u2 * ct_u2);
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        check_param_ct(param, {.from = phase_a,
                               .y = {{{0.0, 0.0, 0.0}, {0.0, balancer, -balancer}, {0.0, -balancer, balancer}}},
                               .from_connected = false});
        CHECK(cabs(param.ytt()(0, 0)) > 0.0);
        check_all_zero({.value = {param.yff(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
    }
    SUBCASE("From side open with magnetising admittance: the primary shunt stays in the balancer") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.from_status = 0;
        auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances({});
        // Schur complement of the primary: Y_hh - Y_h0 Y_00^-1 Y_0h with Y_00 = 1 / z_f + y_m / 2 and
        // Y_01 = Y_02 = -1 / (2 z_f), so every entry of the 2x2 block of the halves drops by
        // Y_01^2 / Y_00 = 1 / (2 z_f (2 + y_m z_f))
        DoubleComplex const drop = 1.0 / (2.0 * z_f * (2.0 + y_m * z_f));
        double const to_scale = test_sn / (0.25 * ct_u2 * ct_u2);
        DoubleComplex const self_half = ((z_0 + z_b) / (2.0 * z_b * z_f) + 0.25 * y_m - drop) * to_scale;
        DoubleComplex const between_halves = (-z_0 / (2.0 * z_b * z_f) - drop) * to_scale;
        Block3 const y{{{0.0, 0.0, 0.0}, {0.0, self_half, between_halves}, {0.0, between_halves, self_half}}};
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        check_param_ct(param, {.from = phase_a, .y = y, .from_connected = false});
        check_all_zero({.value = {param.yff(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
        // the shunt makes it differ from the balancer without magnetising admittance
        input.i0 = 0.0;
        input.p0 = 0.0;
        CHECK_FALSE(params_close(param, make_ct_transformer(input).calc_param<asymmetric_t>()));

        // with a centre-tap electrode: the same primitive block, the electrode as a scalar Schur complement
        DoubleComplex const z_grounding{0.5, 0.05};
        input = ct_input(a);
        input.from_status = 0;
        input.r_grounding_to = z_grounding.real();
        input.x_grounding_to = z_grounding.imag();
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = phase_a, .y = y, .z_grounding = z_grounding, .from_connected = false});
    }
    SUBCASE("From side open with a centre-tap electrode: the balancer through the electrode") {
        DoubleComplex const z_grounding{0.5, 0.05};
        SinglePhaseTransformerInput input = ct_input(a);
        input.i0 = 0.0;
        input.p0 = 0.0;
        input.from_status = 0;
        input.r_grounding_to = z_grounding.real();
        input.x_grounding_to = z_grounding.imag();
        auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances(no_shunt);
        // the halves are the balancer g w w^T with g = 1 / (2 z_b) and w = (1, -1), so the nodal admittance over
        // (a, N, b) is g d d^T with d = (1, -2, 1); the electrode y_e = 1 / z_g adds to Y_NN = 4 g, and eliminating
        // N leaves g - (2 g)^2 / (4 g + y_e) = g y_e / (4 g + y_e) in all four entries of the (a, b) block
        DoubleComplex const g = 1.0 / (2.0 * z_b) * test_sn / (0.25 * ct_u2 * ct_u2);
        DoubleComplex const y_e = 1.0 / z_grounding;
        DoubleComplex const y_ab =
            g * y_e / (4.0 * g + y_e) * (ct_u_rated_to / sqrt3) * (ct_u_rated_to / sqrt3) / test_base_power_1p;
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        for (Idx row = 0; row != 2; ++row) {
            for (Idx col = 0; col != 2; ++col) {
                CAPTURE(row);
                CAPTURE(col);
                CHECK(close_to(param.ytt()(row, col), y_ab));
            }
        }
        check_all_zero({.value = {param.yff(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
        // the generic expected parameter agrees
        double const half = 0.5 * ct_u2;
        DoubleComplex const balancer = 1.0 / (2.0 * z_b) * test_sn / (half * half);
        check_param_ct(param, {.from = phase_a,
                               .y = {{{0.0, 0.0, 0.0}, {0.0, balancer, -balancer}, {0.0, -balancer, balancer}}},
                               .z_grounding = z_grounding,
                               .from_connected = false});
    }
    SUBCASE("To side open without magnetising admittance is exactly zero") {
        SinglePhaseTransformerInput input = ct_input(ab);
        input.i0 = 0.0;
        input.p0 = 0.0;
        input.to_status = 0;
        input.r_grounding_to = 0.5;
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        // only the low susceptance of the phase-to-phase primary
        check_all_zero({.value = {param.ytt(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
        input.from_connection = a;
        check_all_zero(make_ct_transformer(input).calc_param<asymmetric_t>());
    }
    SUBCASE("To side open with magnetising admittance: the open-end shunt of the primary") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.to_status = 0;
        auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances({});
        // Y_00 - [Y_01 Y_02] [[Y_11, Y_12], [Y_12, Y_22]]^-1 [Y_10 Y_20]^T with Y_01 = Y_02 and the eigenvalue
        // Y_11 + Y_12 = 1 / (2 z_f) + y_m / 4 of (1, 1): y_m / 2 + (z_f + 2 / y_m)^-1, the open-end shunt of the
        // two-winding transformer
        DoubleComplex const open_end = (0.5 * y_m + 1.0 / (z_f + 2.0 / y_m)) * test_sn / (test_u1 * test_u1);
        check_param_ct(make_ct_transformer(input).calc_param<asymmetric_t>(),
                       {.from = phase_a, .y = {{{open_end, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}}});
    }
    SUBCASE("To side open with an electrode and magnetising admittance: the electrode is floating") {
        SinglePhaseTransformerInput input = ct_input(ab);
        input.to_status = 0;
        input.r_grounding_to = 0.5;
        input.x_grounding_to = 0.05;
        auto const [z_f, z_hs, z_0, z_b, y_m] = ct_impedances({});
        // the open halves leave the centre tap with only the electrode (Y_NN = y_e, no coupling), so the Schur
        // complement of the electrode is zero and the open-end shunt of the primary is unchanged
        DoubleComplex const open_end = (0.5 * y_m + 1.0 / (z_f + 2.0 / y_m)) * test_sn / (test_u1 * test_u1);
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        check_param_ct(param, {.from = connection_case(ab),
                               .y = {{{open_end, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}},
                               .z_grounding = {0.5, 0.05}});
        check_all_zero({.value = {param.ytt(), param.yft(), param.ytf(), ComplexTensor<asymmetric_t>{}}});
    }
    SUBCASE("Both sides open") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.from_status = 0;
        input.to_status = 0;
        check_all_zero(make_ct_transformer(input).calc_param<asymmetric_t>());
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, taps") {
    using enum WindingConnection1P;
    ConnectionCase const phase_bc = connection_case(bc);
    auto const tapped_input = [](BranchSide tap_side, IntS tap_pos, double tap_size) {
        SinglePhaseTransformerInput input = ct_input(bc);
        input.tap_side = tap_side;
        input.tap_pos = tap_pos;
        input.tap_min = -5;
        input.tap_max = 5;
        input.tap_nom = 0;
        input.tap_size = tap_size;
        return input;
    };

    SUBCASE("From side tap scales the primary voltage") {
        check_param_ct(make_ct_transformer(tapped_input(BranchSide::from, 2, 200.0)).calc_param<asymmetric_t>(),
                       {.from = phase_bc, .y = primitive_ct({.u1 = test_u1 + 400.0})});
    }
    SUBCASE("To side tap scales both halves") {
        check_param_ct(make_ct_transformer(tapped_input(BranchSide::to, -3, 5.0)).calc_param<asymmetric_t>(),
                       {.from = phase_bc, .y = primitive_ct({.u2 = ct_u2 - 15.0})});
    }
    SUBCASE("Tap dependent short circuit parameters change the primary branch only") {
        SinglePhaseTransformerInput input = tapped_input(BranchSide::from, 2, 200.0);
        input.uk_max = 0.05;
        input.pk_max = 800.0;
        input.uk_half = 0.03;
        double const uk = test_uk + 2.0 * (0.05 - test_uk) / 5.0;
        double const pk = test_pk + 2.0 * (800.0 - test_pk) / 5.0;
        // z_b from the nominal uk and pk, z_0 = z_f(tap) - z_b / 2
        CenterTappedData const data{.u1 = test_u1 + 400.0, .uk = uk, .pk = pk, .uk_half = 0.03};
        auto const param = make_ct_transformer(input).calc_param<asymmetric_t>();
        check_param_ct(param, {.from = phase_bc, .y = primitive_ct(data)});
        // the check discriminates: a half branch from the tap-adjusted uk and pk would differ
        CenterTappedData const tapped_half_branch{
            .u1 = test_u1 + 400.0, .uk = uk, .pk = pk, .uk_nominal = uk, .pk_nominal = pk, .uk_half = 0.03};
        CHECK_FALSE(params_close(param, expected_param_ct({.from = phase_bc, .y = primitive_ct(tapped_half_branch)})));
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, loading") {
    using enum WindingConnection1P;
    SinglePhaseTransformer const transformer = make_ct_transformer(ct_input(a));
    // one leg at sn / 2 is full load: each half is rated sn / 2
    double const half_rating = 0.5 * test_sn / test_base_power_1p;
    DoubleComplex const s_leg = std::polar(half_rating, 0.3);

    SUBCASE("A single leg at sn / 2") {
        BranchSolverOutput<asymmetric_t> const solver_output{.s_f = ComplexValue<asymmetric_t>{0.5 * s_leg, 0.0, 0.0},
                                                             .s_t = ComplexValue<asymmetric_t>{-s_leg, 0.0, 0.0},
                                                             .i_f = ComplexValue<asymmetric_t>{0.1, 0.0, 0.0},
                                                             .i_t = ComplexValue<asymmetric_t>{1.0, 0.0, 0.0}};
        CHECK(transformer.get_output<asymmetric_t>(solver_output).loading == doctest::Approx(1.0));
    }
    SUBCASE("Balanced full load") {
        BranchSolverOutput<asymmetric_t> const solver_output{.s_f = ComplexValue<asymmetric_t>{2.0 * s_leg, 0.0, 0.0},
                                                             .s_t = ComplexValue<asymmetric_t>{-s_leg, -s_leg, 0.0},
                                                             .i_f = ComplexValue<asymmetric_t>{0.1, 0.0, 0.0},
                                                             .i_t = ComplexValue<asymmetric_t>{1.0, 1.0, 0.0}};
        CHECK(transformer.get_output<asymmetric_t>(solver_output).loading == doctest::Approx(1.0));
    }
    SUBCASE("The primary side uses the complex sum") {
        SinglePhaseTransformer const ab_primary = make_ct_transformer(ct_input(ab));
        DoubleComplex const s_a{0.06, 0.02};
        DoubleComplex const s_b{0.05, -0.01};
        BranchSolverOutput<asymmetric_t> const solver_output{.s_f = ComplexValue<asymmetric_t>{s_a, s_b, 0.0},
                                                             .s_t = ComplexValue<asymmetric_t>{0.01, -0.01, 0.0},
                                                             .i_f = ComplexValue<asymmetric_t>{0.1, 0.1, 0.0},
                                                             .i_t = ComplexValue<asymmetric_t>{1.0, 1.0, 0.0}};
        CHECK(ab_primary.get_output<asymmetric_t>(solver_output).loading ==
              doctest::Approx(cabs(s_a + s_b) * test_base_power_1p / test_sn));
    }
    SUBCASE("Symmetric") {
        BranchSolverOutput<symmetric_t> const solver_output{
            .s_f = DoubleComplex{0.01, 0.002}, .s_t = DoubleComplex{-0.02, 0.0}, .i_f = 0.01, .i_t = 0.02};
        CHECK(transformer.get_output<symmetric_t>(solver_output).loading == doctest::Approx(0.02 * 1.0e6 / test_sn));
    }
}

TEST_CASE("Test single phase transformer - centre-tapped secondary, invalid input") {
    using enum WindingConnection1P;

    SUBCASE("Centre-tapped from side") {
        CHECK_THROWS_WITH_AS(make_transformer(test_input(ab_center_tapped, a)),
                             "Invalid phase connection for single_phase_transformer 1 with phases (ab): "
                             "a centre-tapped winding is only allowed on the to side\n",
                             InvalidPhaseConnection);
        CHECK_THROWS_AS(make_transformer(test_input(ab_center_tapped, ab_center_tapped)), InvalidPhaseConnection);
    }
    SUBCASE("Half-winding parameters on a winding that is not centre-tapped") {
        SinglePhaseTransformerInput input = test_input(a, a);
        input.uk_half = 0.03;
        CHECK_THROWS_WITH_AS(make_transformer(input),
                             "Invalid phase connection for single_phase_transformer 1 with phases (a): "
                             "uk_half refers to a centre-tapped winding and must be left empty\n",
                             InvalidPhaseConnection);
        input = test_input(ab, bc);
        input.pk_half = 0.0;
        CHECK_THROWS_WITH_AS(make_transformer(input),
                             "Invalid phase connection for single_phase_transformer 1 with phases (bc): "
                             "pk_half refers to a centre-tapped winding and must be left empty\n",
                             InvalidPhaseConnection);
    }
    SUBCASE("An earthing impedance on the centre tap is accepted") {
        SinglePhaseTransformerInput input = ct_input(a);
        input.r_grounding_to = 1.0;
        input.x_grounding_to = 0.1;
        CHECK_NOTHROW(make_ct_transformer(input));
    }
}

TEST_CASE("Test kron_eliminate") {
    using Matrix4 = Eigen::Matrix<DoubleComplex, 4, 4>;
    using Matrix5 = Eigen::Matrix<DoubleComplex, 5, 5>;
    ComplexTensor4 const y4{DoubleComplex{4.0, -8.0},  DoubleComplex{5.0, -9.0}, DoubleComplex{3.5, -7.0},
                            DoubleComplex{6.0, -10.0}, DoubleComplex{-1.0, 2.0}, DoubleComplex{-0.5, 1.5},
                            DoubleComplex{-2.0, 3.0},  DoubleComplex{-1.2, 2.2}, DoubleComplex{-1.5, 2.5},
                            DoubleComplex{-1.8, 2.8}};

    SUBCASE("Eliminating the last node of a 4x4 equals kron_reduction") {
        Matrix4 y = y4.matrix();
        kron_eliminate(y, 3);
        ComplexTensor<asymmetric_t> const reduced = kron_reduction(y4);
        for (Idx row = 0; row != 4; ++row) {
            for (Idx col = 0; col != 4; ++col) {
                CAPTURE(row);
                CAPTURE(col);
                if (row == 3 || col == 3) {
                    CHECK(y(row, col) == 0.0);
                } else {
                    CHECK(cabs(y(row, col) - reduced(row, col)) < test_tolerance);
                }
            }
        }
    }

    SUBCASE("Eliminating two nodes one after the other equals the Schur complement") {
        Matrix5 y_full;
        for (Idx row = 0; row != 5; ++row) {
            for (Idx col = 0; col != 5; ++col) {
                y_full(row, col) = row == col ? DoubleComplex{3.0 + static_cast<double>(row), -6.0}
                                              : DoubleComplex{-0.3 * static_cast<double>(row + 1),
                                                              0.7 + 0.1 * static_cast<double>(col)};
            }
        }
        Matrix5 y = y_full;
        kron_eliminate(y, 1);
        kron_eliminate(y, 3);

        // A - B D^-1 C with the kept nodes {0, 2, 4} and the eliminated nodes {1, 3}
        std::array<Idx, 3> const kept{0, 2, 4};
        Block2 const d_inv = inverse({{{y_full(1, 1), y_full(1, 3)}, {y_full(3, 1), y_full(3, 3)}}});
        std::array<Idx, 2> const eliminated{1, 3};
        for (Idx const row : kept) {
            for (Idx const col : kept) {
                DoubleComplex expected = y_full(row, col);
                for (size_t p = 0; p != 2; ++p) {
                    for (size_t q = 0; q != 2; ++q) {
                        expected -= y_full(row, eliminated[p]) * d_inv[p][q] * y_full(eliminated[q], col);
                    }
                }
                CAPTURE(row);
                CAPTURE(col);
                CHECK(cabs(y(row, col) - expected) < test_tolerance);
            }
        }
        for (Idx const k : eliminated) {
            CHECK(y.row(k).isZero(0.0));
            CHECK(y.col(k).isZero(0.0));
        }
    }
}

} // namespace power_grid_model
