// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

// Power flow on networks with one- and two-phase buses (absent phase slots)
//
// The expected values are derived by hand with Kirchhoff's current law: the voltages are chosen, the branch, source
// and shunt currents follow from the admittances, and the loads are whatever injection closes the current balance.
// The solvers must reproduce the chosen voltages, with the absent phase slots at exactly zero.

#include "test_math_solver_common.hpp"
#include "test_math_solver_pf.hpp" // NOLINT(misc-include-cleaner)

#include <power_grid_model/calculation_parameters.hpp>
#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/dummy_logging.hpp>
#include <power_grid_model/common/enum.hpp>
#include <power_grid_model/common/exception.hpp>
#include <power_grid_model/common/grouped_index_vector.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>
#include <power_grid_model/math_solver/common_solver_functions.hpp>
#include <power_grid_model/math_solver/iterative_current_pf_solver.hpp>
#include <power_grid_model/math_solver/linear_pf_solver.hpp>
#include <power_grid_model/math_solver/math_solver.hpp>
#include <power_grid_model/math_solver/newton_raphson_pf_solver.hpp>
#include <power_grid_model/math_solver/y_bus.hpp>

#include <doctest/doctest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <memory>
#include <optional>
#include <ranges>
#include <vector>

namespace power_grid_model::math_solver {
namespace {
// solver type and the initialization of the power flow
template <typename SolverType, CalculationInitialization init> struct SolverCase {
    using Solver = SolverType;
    static constexpr CalculationInitialization initialization = init;
};
using NewtonRaphsonLinearStart =
    SolverCase<NewtonRaphsonPFSolver<asymmetric_t>, CalculationInitialization::default_initialization>;
using NewtonRaphsonFlatStart = SolverCase<NewtonRaphsonPFSolver<asymmetric_t>, CalculationInitialization::flat>;
using NewtonRaphsonAverageSourceStart =
    SolverCase<NewtonRaphsonPFSolver<asymmetric_t>, CalculationInitialization::average_source>;
using IterativeCurrent =
    SolverCase<IterativeCurrentPFSolver<asymmetric_t>, CalculationInitialization::default_initialization>;
using Linear = SolverCase<LinearPFSolver<asymmetric_t>, CalculationInitialization::default_initialization>;
} // namespace
} // namespace power_grid_model::math_solver

TYPE_TO_STRING_AS("NewtonRaphsonPFSolver<asymmetric_t>, linear start",
                  power_grid_model::math_solver::NewtonRaphsonLinearStart);
TYPE_TO_STRING_AS("NewtonRaphsonPFSolver<asymmetric_t>, flat start",
                  power_grid_model::math_solver::NewtonRaphsonFlatStart);
TYPE_TO_STRING_AS("NewtonRaphsonPFSolver<asymmetric_t>, average source start",
                  power_grid_model::math_solver::NewtonRaphsonAverageSourceStart);
TYPE_TO_STRING_AS("IterativeCurrentPFSolver<asymmetric_t>", power_grid_model::math_solver::IterativeCurrent);
TYPE_TO_STRING_AS("LinearPFSolver<asymmetric_t>", power_grid_model::math_solver::Linear);

namespace power_grid_model::math_solver {
namespace {
using common::logging::NoLogger;
using Value = ComplexValue<asymmetric_t>;
using Tensor = ComplexTensor<asymmetric_t>;

constexpr double err_tol = 1e-13;
constexpr Idx max_iter = 50;
constexpr double result_tol = 1e-10;

constexpr double deg = deg_30 / 30.0;
constexpr DoubleComplex u_ref{1.05};
constexpr DoubleComplex y_ref{10.0, -50.0};

// symmetric source voltage, also the voltage of a phase that no branch loads
Value source_voltage() { return Value{u_ref}; }

Tensor single_phase_branch(DoubleComplex y) {
    Tensor result{0.0};
    result(0, 0) = y;
    return result;
}

Tensor two_phase_branch(DoubleComplex y_self, DoubleComplex y_mutual) {
    Tensor result{0.0};
    result(0, 0) = y_self;
    result(1, 1) = y_self;
    result(0, 1) = y_mutual;
    result(1, 0) = y_mutual;
    return result;
}

// source -- bus0 (abc) --branch-- bus1 (mask1)
// - the loads of each bus are split into const_pq, const_i and const_y for the iterative solvers, or are one const_y
//   load for the linear solver
// - the injections on the absent phases of bus1 are given by absent_injection; they must be ignored
struct ReducedPhaseCase {
    Tensor y_branch;
    PhaseMask mask1;
    Value u0;
    Value u1;
    SourceCalcParam source{.y1 = y_ref, .y0 = y_ref};
    std::optional<Tensor> y_shunt{};
    Value absent_injection{};
};

struct ReducedPhaseGrid {
    MathModelTopology topo;
    MathModelParam<asymmetric_t> param;
    PowerFlowInput<asymmetric_t> input;
    SolverOutput<asymmetric_t> expected;
};

// injection that scales with the voltage for a const_i (power 1) or const_y (power 2) load, only on present phases
Value scaled_injection(Value const& s, Value const& u, PhaseMask phases, int power) {
    Value result = s;
    for (Idx phase = 0; phase != 3; ++phase) {
        if (phases.has(phase)) {
            for (int k = 0; k != power; ++k) {
                result(phase) /= std::abs(u(phase));
            }
        }
    }
    return result;
}

ReducedPhaseGrid make_grid(ReducedPhaseCase const& rpc, bool const_y_only) {
    using enum LoadGenType;
    ReducedPhaseGrid grid{};
    Idx const n_load_per_bus = const_y_only ? 1 : 3;
    Idx const n_shunt = rpc.y_shunt.has_value() ? 1 : 0;

    // topology
    grid.topo.slack_bus = 0;
    grid.topo.phase_shift = {0.0, 0.0};
    grid.topo.branch_bus_idx = {{0, 1}};
    grid.topo.sources_per_bus = {from_sparse, {0, 1, 1}};
    grid.topo.shunts_per_bus = {from_sparse, {0, 0, n_shunt}};
    grid.topo.load_gens_per_bus = {from_sparse, {0, n_load_per_bus, 2 * n_load_per_bus}};
    grid.topo.load_gen_type = const_y_only
                                  ? std::vector<LoadGenType>{const_y, const_y}
                                  : std::vector<LoadGenType>{const_pq, const_i, const_y, const_pq, const_i, const_y};
    grid.topo.voltage_regulators_per_load_gen = {from_sparse,
                                                 IdxVector(static_cast<size_t>(2 * n_load_per_bus + 1), 0)};
    grid.topo.has_reduced_phase = true;
    grid.topo.bus_phase_mask = {phase_mask_abc, rpc.mask1};
    grid.topo.bus_slot_pattern = {SlotPattern::three_phase, SlotPattern::three_phase};

    // parameters
    grid.param.branch_param = {{rpc.y_branch, -rpc.y_branch, -rpc.y_branch, rpc.y_branch}};
    if (rpc.y_shunt.has_value()) {
        grid.param.shunt_param = {rpc.y_shunt.value()};
    }
    grid.param.source_param = {rpc.source};

    // Kirchhoff: currents into each bus from source, shunt and loads equal the currents out into the branch
    Value const i_f = dot(rpc.y_branch, Value{rpc.u0 - rpc.u1});
    Value const i_t = -i_f;
    Value const i_source = dot(rpc.source.y_ref<asymmetric_t>(), Value{source_voltage() - rpc.u0});
    Value i_shunt{}; // injection of the shunt into bus1, only on the present phases
    if (rpc.y_shunt.has_value()) {
        for (Idx row = 0; row != 3; ++row) {
            for (Idx col = 0; col != 3; ++col) {
                if (rpc.mask1.has(row) && rpc.mask1.has(col)) {
                    i_shunt(row) -= rpc.y_shunt.value()(row, col) * rpc.u1(col);
                }
            }
        }
    }
    Value const i_load0 = i_f - i_source;
    Value const i_load1 = i_t - i_shunt;
    Value const s_load0 = rpc.u0 * conj(i_load0);
    Value const s_load1 = rpc.u1 * conj(i_load1); // exactly zero on the absent phases, where u1 is zero

    // input
    grid.input.source = {u_ref};
    grid.input.initialization = CalculationInitialization::default_initialization;
    auto const add_bus_loads = [&grid, const_y_only](Value const& s, Value const& u, PhaseMask phases,
                                                     Value const& absent_injection) {
        auto const with_absent = [&phases, &absent_injection](Value value) {
            for (Idx phase = 0; phase != 3; ++phase) {
                if (!phases.has(phase)) {
                    value(phase) = absent_injection(phase);
                }
            }
            return value;
        };
        if (const_y_only) {
            grid.input.s_injection.push_back(with_absent(scaled_injection(s, u, phases, 2)));
            return;
        }
        Value const s_third = s / 3.0;
        for (int power = 0; power != 3; ++power) {
            grid.input.s_injection.push_back(with_absent(scaled_injection(s_third, u, phases, power)));
        }
    };
    add_bus_loads(s_load0, rpc.u0, phase_mask_abc, Value{});
    add_bus_loads(s_load1, rpc.u1, rpc.mask1, rpc.absent_injection);

    // expected output
    auto& expected = grid.expected;
    expected.u = {rpc.u0, rpc.u1};
    expected.branch = {{.s_f = rpc.u0 * conj(i_f), .s_t = rpc.u1 * conj(i_t), .i_f = i_f, .i_t = i_t}};
    expected.source = {{.s = rpc.u0 * conj(i_source), .i = i_source}};
    expected.bus_injection = {rpc.u0 * conj(i_f), rpc.u1 * conj(Value{i_t - i_shunt})};
    if (rpc.y_shunt.has_value()) {
        expected.shunt = {{.s = rpc.u1 * conj(i_shunt), .i = i_shunt}};
    }
    for (Idx load = 0; load != n_load_per_bus; ++load) {
        expected.load_gen.push_back(
            {.s = s_load0 / static_cast<double>(n_load_per_bus), .i = i_load0 / static_cast<double>(n_load_per_bus)});
    }
    for (Idx load = 0; load != n_load_per_bus; ++load) {
        expected.load_gen.push_back(
            {.s = s_load1 / static_cast<double>(n_load_per_bus), .i = i_load1 / static_cast<double>(n_load_per_bus)});
    }
    return grid;
}

template <typename SolverCaseType>
SolverOutput<asymmetric_t> run_case(MathModelTopology const& topo, MathModelParam<asymmetric_t> const& param,
                                    PowerFlowInput<asymmetric_t> input) {
    using Solver = SolverCaseType::Solver;
    input.initialization = SolverCaseType::initialization;
    YBus<asymmetric_t> const y_bus{topo, param};
    Solver solver{y_bus, topo};
    NoLogger log;
    return run_power_flow(solver, y_bus, input, err_tol, max_iter, log);
}

void check_absent_slots_zero(SolverOutput<asymmetric_t> const& output, MathModelTopology const& topo) {
    for (Idx bus = 0; bus != topo.n_bus(); ++bus) {
        PhaseMask const phases = topo.bus_phase(bus);
        for (Idx phase = 0; phase != 3; ++phase) {
            if (!phases.has(phase)) {
                CHECK(output.u[bus](phase) == 0.0);
                // a positive zero, so that the node output u_angle is arg(0) = 0 and not pi
                CHECK(std::arg(output.u[bus](phase)) == 0.0);
                CHECK(output.bus_injection[bus](phase) == 0.0);
            }
        }
    }
    for (auto const& [bus, load_gens] : enumerated_zip_sequence(topo.load_gens_per_bus)) {
        PhaseMask const phases = topo.bus_phase(bus);
        for (Idx const load_gen : load_gens) {
            for (Idx phase = 0; phase != 3; ++phase) {
                if (!phases.has(phase)) {
                    CHECK(output.load_gen[load_gen].s(phase) == 0.0);
                    CHECK(output.load_gen[load_gen].i(phase) == 0.0);
                }
            }
        }
    }
}

void check_reduced_phase_grid(ReducedPhaseGrid const& grid, SolverOutput<asymmetric_t> const& output) {
    REQUIRE(output.u.size() == grid.expected.u.size());
    REQUIRE(output.load_gen.size() == grid.expected.load_gen.size());
    REQUIRE(output.shunt.size() == grid.expected.shunt.size());
    assert_output(output, grid.expected, false, result_tol);
    check_absent_slots_zero(output, grid.topo);
    PhaseMask const mask1 = grid.topo.bus_phase(1);
    for (Idx phase = 0; phase != 3; ++phase) {
        if (!mask1.has(phase)) {
            CHECK(output.branch[0].i_t(phase) == 0.0);
            CHECK(output.branch[0].s_t(phase) == 0.0);
            for (auto const& shunt : output.shunt) {
                CHECK(shunt.i(phase) == 0.0);
                CHECK(shunt.s(phase) == 0.0);
            }
        }
    }
}

// U0a = 1.02 at -1 degree, U1a = 0.95 at -4 degree; phases b and c of bus0 are not loaded by any branch
ReducedPhaseCase single_phase_case() {
    return {.y_branch = single_phase_branch({2.0, -4.0}),
            .mask1 = phase_mask_a,
            .u0 = {1.02 * std::polar(1.0, -deg), source_voltage()(1), source_voltage()(2)},
            .u1 = {0.95 * std::polar(1.0, -4.0 * deg), 0.0, 0.0}};
}

// U1c = 0, phase c of bus0 is not loaded by any branch
ReducedPhaseCase two_phase_case() {
    return {.y_branch = two_phase_branch({2.0, -4.0}, {-0.4, 0.8}),
            .mask1 = phase_mask_ab,
            .u0 = {1.02 * std::polar(1.0, -deg), 1.03 * std::polar(1.0, -121.5 * deg), source_voltage()(2)},
            .u1 = {0.95 * std::polar(1.0, -4.0 * deg), 0.96 * std::polar(1.0, -125.0 * deg), 0.0}};
}

// meshed network with a single phase bus that is eliminated first, so that the LU factorization has a fill-in
// between two three-phase buses:
//   bus0 (a) --a-- bus1 (abc) --abc-- bus3 (abc, source)
//   bus0 (a) --a-- bus2 (abc) --abc-- bus3
// - loads: bus0 const_pq, const_i, const_y (with junk on the absent phases), bus1 const_pq, bus2 const_pq generator;
//   all const_y for the linear solver
// - with_regulator: the generator at bus2 regulates the voltage
MathModelTopology make_meshed_topo(bool const_y_only, bool with_regulator) {
    using enum LoadGenType;
    MathModelTopology topo;
    topo.slack_bus = 3;
    topo.is_radial = false;
    topo.phase_shift = {0.0, 0.0, 0.0, 0.0};
    topo.branch_bus_idx = {{0, 1}, {0, 2}, {1, 3}, {2, 3}};
    topo.fill_in = {{1, 2}};
    topo.sources_per_bus = {from_sparse, {0, 0, 0, 0, 1}};
    topo.shunts_per_bus = {from_sparse, {0, 0, 0, 0, 0}};
    topo.load_gens_per_bus = {from_sparse, {0, 3, 4, 5, 5}};
    topo.load_gen_type = const_y_only ? std::vector<LoadGenType>(5, const_y)
                                      : std::vector<LoadGenType>{const_pq, const_i, const_y, const_pq, const_pq};
    topo.voltage_regulators_per_load_gen = {from_sparse,
                                            with_regulator ? IdxVector{0, 0, 0, 0, 0, 1} : IdxVector{0, 0, 0, 0, 0, 0}};
    topo.has_reduced_phase = true;
    topo.bus_phase_mask = {phase_mask_a, phase_mask_abc, phase_mask_abc, phase_mask_abc};
    topo.bus_slot_pattern.assign(4, SlotPattern::three_phase);
    return topo;
}

MathModelParam<asymmetric_t> make_meshed_param() {
    MathModelParam<asymmetric_t> param;
    Tensor const y_a = single_phase_branch({3.0, -6.0});
    Tensor const y_abc{DoubleComplex{4.0, -9.0}, DoubleComplex{-0.5, 1.0}};
    param.branch_param = {
        {y_a, -y_a, -y_a, y_a}, {y_a, -y_a, -y_a, y_a}, {y_abc, -y_abc, -y_abc, y_abc}, {y_abc, -y_abc, -y_abc, y_abc}};
    param.source_param = {SourceCalcParam{.y1 = y_ref, .y0 = DoubleComplex{5.0, -30.0}}};
    return param;
}

PowerFlowInput<asymmetric_t> make_meshed_input(std::optional<double> q_max) {
    PowerFlowInput<asymmetric_t> input;
    input.source = {u_ref};
    // the loads of bus0 have junk on the absent phases b and c, including negative values
    input.s_injection = {Value{{-0.1, -0.05}, {-1.0, 0.5}, {0.7, -0.3}},
                         Value{{-0.05, -0.02}, {-0.3, 0.2}, {0.2, -0.1}},
                         Value{{-0.05, -0.01}, {-0.3, -0.2}, {-0.2, -0.1}},
                         Value{{-0.2, -0.1}, {-0.15, -0.05}, {-0.1, -0.02}}, Value{{0.1, 0.0}, {0.1, 0.0}, {0.1, 0.0}}};
    input.load_gen_status = {1, 1, 1, 1, 1};
    if (q_max.has_value()) {
        input.voltage_regulator = {
            {.status = 1, .u_ref = 1.06, .q_min = -q_max.value(), .q_max = q_max.value(), .generator_id = 1}};
    }
    return input;
}

// the power flow equation on every present phase: the specified (load_gen and source) injections equal the bus
// injection, which equals the sum of the flows into the branches
void check_meshed_power_balance(SolverOutput<asymmetric_t> const& output, MathModelTopology const& topo) {
    std::vector<Value> branch_flow(static_cast<size_t>(topo.n_bus()), Value{});
    for (auto const& [branch, buses] : std::views::zip(output.branch, topo.branch_bus_idx)) {
        branch_flow[static_cast<size_t>(buses[0])] += branch.s_f;
        branch_flow[static_cast<size_t>(buses[1])] += branch.s_t;
    }
    for (auto const& [bus, sources, load_gens] :
         enumerated_zip_sequence(topo.sources_per_bus, topo.load_gens_per_bus)) {
        CAPTURE(bus);
        Value specified{};
        for (Idx const source : sources) {
            specified += output.source[source].s;
        }
        for (Idx const load_gen : load_gens) {
            specified += output.load_gen[load_gen].s;
        }
        PhaseMask const phases = topo.bus_phase(bus);
        for (Idx phase = 0; phase != 3; ++phase) {
            if (phases.has(phase)) {
                CAPTURE(phase);
                CHECK(cabs(specified(phase) - output.bus_injection[bus](phase)) < result_tol);
                CHECK(cabs(branch_flow[static_cast<size_t>(bus)](phase) - output.bus_injection[bus](phase)) <
                      result_tol);
            }
        }
    }
}

void check_solver_output_exact(SolverOutput<asymmetric_t> const& output, SolverOutput<asymmetric_t> const& ref) {
    auto const equal = [](Value const& x, Value const& y) { return (x == y).all(); };
    REQUIRE(output.u.size() == ref.u.size());
    for (size_t i = 0; i != ref.u.size(); ++i) {
        CHECK(equal(output.u[i], ref.u[i]));
        CHECK(equal(output.bus_injection[i], ref.bus_injection[i]));
    }
    REQUIRE(output.branch.size() == ref.branch.size());
    for (size_t i = 0; i != ref.branch.size(); ++i) {
        CHECK(equal(output.branch[i].s_f, ref.branch[i].s_f));
        CHECK(equal(output.branch[i].s_t, ref.branch[i].s_t));
        CHECK(equal(output.branch[i].i_f, ref.branch[i].i_f));
        CHECK(equal(output.branch[i].i_t, ref.branch[i].i_t));
    }
    REQUIRE(output.source.size() == ref.source.size());
    for (size_t i = 0; i != ref.source.size(); ++i) {
        CHECK(equal(output.source[i].s, ref.source[i].s));
        CHECK(equal(output.source[i].i, ref.source[i].i));
    }
    REQUIRE(output.shunt.size() == ref.shunt.size());
    for (size_t i = 0; i != ref.shunt.size(); ++i) {
        CHECK(equal(output.shunt[i].s, ref.shunt[i].s));
        CHECK(equal(output.shunt[i].i, ref.shunt[i].i));
    }
    REQUIRE(output.load_gen.size() == ref.load_gen.size());
    for (size_t i = 0; i != ref.load_gen.size(); ++i) {
        CHECK(equal(output.load_gen[i].s, ref.load_gen[i].s));
        CHECK(equal(output.load_gen[i].i, ref.load_gen[i].i));
    }
}
} // namespace

TEST_CASE_TEMPLATE("Test math solver - PF with absent phases", SolverCaseType, NewtonRaphsonLinearStart,
                   NewtonRaphsonFlatStart, NewtonRaphsonAverageSourceStart, IterativeCurrent, Linear) {
    constexpr bool const_y_only = !SolverCaseType::Solver::is_iterative;

    SUBCASE("Single phase bus") {
        ReducedPhaseGrid const grid = make_grid(single_phase_case(), const_y_only);
        auto const output = run_case<SolverCaseType>(grid.topo, grid.param, grid.input);
        check_reduced_phase_grid(grid, output);
        CHECK(output.u[1](1) == 0.0);
        CHECK(output.u[1](2) == 0.0);
    }

    SUBCASE("Two phase bus") {
        ReducedPhaseGrid const grid = make_grid(two_phase_case(), const_y_only);
        auto const output = run_case<SolverCaseType>(grid.topo, grid.param, grid.input);
        check_reduced_phase_grid(grid, output);
        CHECK(output.u[1](2) == 0.0);
    }

    SUBCASE("Injections on absent phases are ignored") {
        ReducedPhaseCase rpc = single_phase_case();
        ReducedPhaseGrid const ref_grid = make_grid(rpc, const_y_only);
        rpc.absent_injection = {0.0, {0.3, 0.1}, {0.0, -0.2}};
        ReducedPhaseGrid const grid = make_grid(rpc, const_y_only);
        CHECK(grid.input.s_injection.back()(1) != 0.0);

        auto const output = run_case<SolverCaseType>(grid.topo, grid.param, grid.input);
        auto const ref_output = run_case<SolverCaseType>(ref_grid.topo, ref_grid.param, ref_grid.input);
        // the expected output is the same as without the injections on the absent phases
        check_reduced_phase_grid(grid, output);
        assert_output(output, ref_output, false, result_tol);
    }

    SUBCASE("Shunt on a single phase bus") {
        // sequence-based shunt with ys != ys0, so it has mutual terms on the absent phases
        ReducedPhaseCase rpc = single_phase_case();
        DoubleComplex const ys{0.1, 0.3};
        DoubleComplex const ys0{0.4, -0.2};
        rpc.y_shunt = Tensor{(2.0 * ys + ys0) / 3.0, (ys0 - ys) / 3.0};
        ReducedPhaseGrid const grid = make_grid(rpc, const_y_only);
        // shunt injection on phase a is -Y_aa U1a
        CHECK(cabs(grid.expected.shunt[0].i(0) + rpc.y_shunt.value()(0, 0) * rpc.u1(0)) < numerical_tolerance);

        auto const output = run_case<SolverCaseType>(grid.topo, grid.param, grid.input);
        check_reduced_phase_grid(grid, output);
    }

    SUBCASE("Equivalent to a full phase network") {
        // source with y1 != y0, so that phases b and c of bus0 depend on the current in phase a
        // in the full phase network, phases b and c of bus1 carry no load, so no current flows in the b and c
        // conductors of the branch, and the a phase results must be the same as in the reduced network
        DoubleComplex const y{2.0, -4.0};
        SourceCalcParam const source{.y1 = y_ref, .y0 = DoubleComplex{5.0, -30.0}};
        Value const s0{{-0.2, -0.1}, {-0.15, -0.05}, {-0.1, -0.02}};
        Value const s1{{-0.3, -0.1}, 0.0, 0.0};

        auto const make_topo = [](bool reduced) {
            using enum LoadGenType;
            MathModelTopology topo;
            topo.slack_bus = 0;
            topo.phase_shift = {0.0, 0.0};
            topo.branch_bus_idx = {{0, 1}};
            topo.sources_per_bus = {from_sparse, {0, 1, 1}};
            topo.shunts_per_bus = {from_sparse, {0, 0, 0}};
            topo.load_gens_per_bus = {from_sparse, {0, 1, 2}};
            topo.load_gen_type = {const_y_only ? const_y : const_pq, const_y_only ? const_y : const_pq};
            topo.voltage_regulators_per_load_gen = {from_sparse, {0, 0, 0}};
            if (reduced) {
                topo.has_reduced_phase = true;
                topo.bus_phase_mask = {phase_mask_abc, phase_mask_a};
                topo.bus_slot_pattern = {SlotPattern::three_phase, SlotPattern::three_phase};
            }
            return topo;
        };
        auto const make_param = [&source](Tensor const& y_branch) {
            MathModelParam<asymmetric_t> param;
            param.branch_param = {{y_branch, -y_branch, -y_branch, y_branch}};
            param.source_param = {source};
            return param;
        };
        PowerFlowInput<asymmetric_t> const input{.source = {u_ref},
                                                 .s_injection = {s0, s1},
                                                 .voltage_regulator = {},
                                                 .load_gen_status = {},
                                                 .initialization = CalculationInitialization::default_initialization};

        MathModelTopology const topo_reduced = make_topo(true);
        MathModelTopology const topo_full = make_topo(false);
        auto const reduced = run_case<SolverCaseType>(topo_reduced, make_param(single_phase_branch(y)), input);
        auto const full = run_case<SolverCaseType>(topo_full, make_param(Tensor{y}), input);

        constexpr double tol = 1e-12;
        check_close<asymmetric_t>(reduced.u[0], full.u[0], tol);
        check_close(reduced.u[1](0), full.u[1](0), tol);
        check_close(reduced.branch[0].i_f(0), full.branch[0].i_f(0), tol);
        check_close(reduced.branch[0].i_t(0), full.branch[0].i_t(0), tol);
        check_close(reduced.branch[0].s_f(0), full.branch[0].s_f(0), tol);
        check_close(reduced.branch[0].s_t(0), full.branch[0].s_t(0), tol);
        check_close<asymmetric_t>(reduced.source[0].i, full.source[0].i, tol);
        check_close<asymmetric_t>(reduced.source[0].s, full.source[0].s, tol);
        check_close<asymmetric_t>(reduced.load_gen[0].s, full.load_gen[0].s, tol);
        check_close<asymmetric_t>(reduced.load_gen[0].i, full.load_gen[0].i, tol);
        check_close(reduced.load_gen[1].s(0), full.load_gen[1].s(0), tol);
        check_close(reduced.load_gen[1].i(0), full.load_gen[1].i(0), tol);
        // the network is not trivially balanced: the b and c phases of bus0 are affected by the a phase
        CHECK(cabs(full.u[0](1) - source_voltage()(1)) > 1e-3);
        // in the full phase network phases b and c of bus1 follow bus0, in the reduced network they are absent
        check_close(full.u[1](1), full.u[0](1), tol);
        CHECK(reduced.u[1](1) == 0.0);
        CHECK(reduced.u[1](2) == 0.0);
    }
}

TEST_CASE_TEMPLATE("Test math solver - PF with all phases flagged as reduced is unchanged", SolverCaseType,
                   NewtonRaphsonLinearStart, NewtonRaphsonFlatStart, NewtonRaphsonAverageSourceStart, IterativeCurrent,
                   Linear) {
    // a math model with has_reduced_phase set but all buses with all phases runs through the gated code paths, which
    // must leave every value bit for bit identical to the unflagged run
    PFSolverTestGrid<asymmetric_t> const grid;
    MathModelTopology const topo_ref = grid.topo();
    MathModelTopology topo_flagged = grid.topo();
    topo_flagged.has_reduced_phase = true;
    topo_flagged.bus_phase_mask.assign(static_cast<size_t>(topo_flagged.n_bus()), phase_mask_abc);
    topo_flagged.bus_slot_pattern.assign(static_cast<size_t>(topo_flagged.n_bus()), SlotPattern::three_phase);

    YBus<asymmetric_t> const y_bus_ref{topo_ref, grid.param()};
    YBus<asymmetric_t> const y_bus_flagged{topo_flagged, grid.param()};
    REQUIRE(y_bus_ref.nnz() == y_bus_flagged.nnz());
    for (Idx k = 0; k != y_bus_ref.nnz(); ++k) {
        CHECK((y_bus_ref.admittance()[k] == y_bus_flagged.admittance()[k]).all());
    }

    for (auto const& input : {grid.pf_input(), grid.pf_input_z()}) {
        auto const output_ref = run_case<SolverCaseType>(topo_ref, grid.param(), input);
        auto const output_flagged = run_case<SolverCaseType>(topo_flagged, grid.param(), input);
        check_solver_output_exact(output_flagged, output_ref);
    }
}

TEST_CASE_TEMPLATE("Test math solver - PF with absent phases in a meshed network with fill-in", SolverCaseType,
                   NewtonRaphsonLinearStart, NewtonRaphsonFlatStart, NewtonRaphsonAverageSourceStart, IterativeCurrent,
                   Linear) {
    constexpr bool const_y_only = !SolverCaseType::Solver::is_iterative;
    MathModelTopology const topo = make_meshed_topo(const_y_only, false);
    MathModelParam<asymmetric_t> const param = make_meshed_param();
    PowerFlowInput<asymmetric_t> const input = make_meshed_input(std::nullopt);

    auto const output = run_case<SolverCaseType>(topo, param, input);
    check_absent_slots_zero(output, topo);
    check_meshed_power_balance(output, topo);
    for (Idx phase = 1; phase != 3; ++phase) {
        for (Idx branch = 0; branch != 2; ++branch) {
            CHECK(output.branch[branch].s_f(phase) == 0.0);
            CHECK(output.branch[branch].i_f(phase) == 0.0);
        }
    }

    // every solver and start agrees with the Newton-Raphson linear start; for the linear solver all loads are const_y,
    // for which the linear solver is exact
    auto const reference = run_case<NewtonRaphsonLinearStart>(topo, param, input);
    assert_output(output, reference, false, result_tol);
}

TEST_CASE("Test math solver - NR with absent phases and Q limits on a regulated bus") {
    // the regulated generator is at the three-phase bus2 of the meshed network
    MathModelTopology const topo = make_meshed_topo(false, true);
    MathModelParam<asymmetric_t> const param = make_meshed_param();

    SUBCASE("Limit not hit") {
        auto const output = run_case<NewtonRaphsonLinearStart>(topo, param, make_meshed_input(100.0));
        for (Idx phase = 0; phase != 3; ++phase) {
            CHECK(cabs(output.u[2](phase)) == doctest::Approx(1.06).epsilon(result_tol));
        }
        CHECK(output.bus[2].q_limit_violated == LimitViolation::none);
        check_absent_slots_zero(output, topo);
        check_meshed_power_balance(output, topo);
    }

    SUBCASE("Limit hit, partial rebuild") {
        auto const output = run_case<NewtonRaphsonLinearStart>(topo, param, make_meshed_input(0.01));
        CHECK(output.bus[2].q_limit_violated != LimitViolation::none);
        for (Idx phase = 0; phase != 3; ++phase) {
            CHECK(cabs(output.u[2](phase)) < 1.06);
        }
        check_absent_slots_zero(output, topo);
        check_meshed_power_balance(output, topo);
    }
}

TEST_CASE("Test math solver - PF with absent phases rejects voltage regulation") {
    using enum LoadGenType;
    // bus0 (abc, source) --a-- bus1 (a, regulated generator)
    MathModelTopology topo;
    topo.slack_bus = 0;
    topo.phase_shift = {0.0, 0.0};
    topo.branch_bus_idx = {{0, 1}};
    topo.sources_per_bus = {from_sparse, {0, 1, 1}};
    topo.shunts_per_bus = {from_sparse, {0, 0, 0}};
    topo.load_gens_per_bus = {from_sparse, {0, 0, 1}};
    topo.load_gen_type = {const_pq};
    topo.voltage_regulators_per_load_gen = {from_sparse, {0, 1}};
    topo.has_reduced_phase = true;
    topo.bus_phase_mask = {phase_mask_abc, phase_mask_a};
    topo.bus_slot_pattern = {SlotPattern::three_phase, SlotPattern::three_phase};

    MathModelParam<asymmetric_t> param;
    Tensor const y_branch = single_phase_branch({10.0, -20.0});
    param.branch_param = {{y_branch, -y_branch, -y_branch, y_branch}};
    param.source_param = {{.y1 = 1e6, .y0 = 1e6}};

    auto const input = [](IntS regulator_status) {
        return PowerFlowInput<asymmetric_t>{
            .source = {1.0},
            .s_injection = {Value{{0.1, 0.0}, 0.0, 0.0}},
            .voltage_regulator =
                {{.status = regulator_status, .u_ref = 1.0, .q_min = nan, .q_max = nan, .generator_id = 42}},
            .load_gen_status = {1}};
    };

    YBus<asymmetric_t> const y_bus{topo, param};
    NoLogger log;
    constexpr bool cache_run = false;

    SUBCASE("Active regulator at a single phase bus") {
        NewtonRaphsonPFSolver<asymmetric_t> solver{y_bus, topo};
        CHECK_THROWS_AS(solver.run_power_flow(y_bus, input(1), 1e-12, 20, cache_run, log), ReducedPhaseNotSupported);
    }
    SUBCASE("Inactive regulator at a single phase bus") {
        NewtonRaphsonPFSolver<asymmetric_t> solver{y_bus, topo};
        SolverOutput<asymmetric_t> output;
        REQUIRE_NOTHROW(output = solver.run_power_flow(y_bus, input(0), 1e-12, 20, cache_run, log));
        CHECK(output.u[1](1) == 0.0);
        CHECK(output.u[1](2) == 0.0);
    }
}

TEST_CASE("Test math solver - split-phase slot patterns of the start values") {
    // bus0 abc; bus1 ab behind a centre-tapped winding (split); bus2 a, a leg line off bus1 (split); bus3 ab behind a
    // phase-to-phase winding (three_phase)
    MathModelTopology topo;
    topo.phase_shift = {0.0, -0.4, -0.4, 0.7};
    topo.has_reduced_phase = true;
    topo.bus_phase_mask = {phase_mask_abc, phase_mask_ab, phase_mask_a, phase_mask_ab};
    topo.bus_slot_pattern = {SlotPattern::three_phase, SlotPattern::split_phase, SlotPattern::split_phase,
                             SlotPattern::three_phase};
    // a three-phase start: the reference angle plus the phase shift of each bus, 120 degrees between the slots
    DoubleComplex const u_ref_start = std::polar(1.02, 0.1);
    ComplexValueVector<asymmetric_t> start;
    for (double const shift : topo.phase_shift) {
        start.emplace_back(u_ref_start * std::exp(1.0i * shift));
    }
    auto const check_unchanged = [&start](ComplexValueVector<asymmetric_t> const& u, Idx bus) {
        for (Idx phase = 0; phase != 3; ++phase) {
            CHECK(u[bus](phase) == start[bus](phase));
        }
    };

    SUBCASE("Slot b of a split-phase bus is opposite to slot a; slot c is left to the mask") {
        ComplexValueVector<asymmetric_t> u = start;
        detail::apply_slot_patterns<asymmetric_t>(topo, u);
        check_unchanged(u, 0);
        check_unchanged(u, 3);
        for (Idx const bus : {1, 2}) {
            CAPTURE(bus);
            CHECK(u[bus](0) == start[bus](0));
            CHECK(u[bus](1) == -start[bus](0));
            CHECK(u[bus](2) == start[bus](2));
            CHECK(std::abs(std::arg(u[bus](1) / u[bus](0))) == doctest::Approx(pi));
        }
        detail::mask_absent_phase_values<asymmetric_t>(topo, u);
        CHECK(u[1](1) == -start[1](0));
        CHECK(u[1](2) == 0.0);
        CHECK(u[2](1) == 0.0);
        CHECK(u[2](2) == 0.0);
    }

    SUBCASE("No-op without reduced phases") {
        // the accessor ignores the pattern vector when the flag is off
        MathModelTopology topo_full = topo;
        topo_full.has_reduced_phase = false;
        ComplexValueVector<asymmetric_t> u = start;
        detail::apply_slot_patterns<asymmetric_t>(topo_full, u);
        for (Idx bus = 0; bus != topo.n_bus(); ++bus) {
            check_unchanged(u, bus);
        }
    }

    SUBCASE("No-op for a symmetric calculation") {
        ComplexValueVector<symmetric_t> const start_sym(4, u_ref_start);
        ComplexValueVector<symmetric_t> u_sym = start_sym;
        detail::apply_slot_patterns<symmetric_t>(topo, u_sym);
        CHECK(u_sym == start_sym);
    }
}

TEST_CASE("Test math solver - reduced phase guards") {
    // source -- bus0 (abc) --a-- bus1 (a)
    ReducedPhaseGrid const grid = make_grid(single_phase_case(), false);
    auto const topo_ptr = std::make_shared<MathModelTopology const>(grid.topo);
    NoLogger log;
    constexpr bool cache_run = false;

    SUBCASE("Asymmetric power flow is supported") {
        YBus<asymmetric_t> const y_bus{*topo_ptr, grid.param};
        MathSolver<asymmetric_t> solver{topo_ptr};
        for (auto const method : {CalculationMethod::newton_raphson, CalculationMethod::iterative_current}) {
            SolverOutput<asymmetric_t> output;
            REQUIRE_NOTHROW(output =
                                solver.run_power_flow(grid.input, err_tol, max_iter, cache_run, log, method, y_bus));
            check_reduced_phase_grid(grid, output);
        }
    }

    SUBCASE("State estimation is not supported") {
        YBus<asymmetric_t> const y_bus{*topo_ptr, grid.param};
        MathSolver<asymmetric_t> solver{topo_ptr};
        StateEstimationInput<asymmetric_t> const se_input{};
        for (auto const method : {CalculationMethod::default_method, CalculationMethod::iterative_linear,
                                  CalculationMethod::newton_raphson}) {
            CHECK_THROWS_AS(solver.run_state_estimation(se_input, err_tol, max_iter, log, method, y_bus),
                            ReducedPhaseNotSupported);
        }
    }

    SUBCASE("Symmetric calculations are not supported") {
        DoubleComplex const y_branch{2.0, -4.0};
        MathModelParam<symmetric_t> param_sym;
        param_sym.branch_param = {{y_branch, -y_branch, -y_branch, y_branch}};
        param_sym.source_param = {{.y1 = y_ref, .y0 = y_ref}};
        YBus<symmetric_t> const y_bus{*topo_ptr, param_sym};
        MathSolver<symmetric_t> solver{topo_ptr};

        PowerFlowInput<symmetric_t> const pf_input{.source = {u_ref},
                                                   .s_injection = ComplexValueVector<symmetric_t>(6, -0.1),
                                                   .voltage_regulator = {},
                                                   .load_gen_status = {},
                                                   .initialization = CalculationInitialization::default_initialization};
        CHECK_THROWS_AS(solver.run_power_flow(pf_input, err_tol, max_iter, cache_run, log,
                                              CalculationMethod::newton_raphson, y_bus),
                        ReducedPhaseNotSupported);

        StateEstimationInput<symmetric_t> const se_input{};
        CHECK_THROWS_AS(
            solver.run_state_estimation(se_input, err_tol, max_iter, log, CalculationMethod::iterative_linear, y_bus),
            ReducedPhaseNotSupported);

        ShortCircuitInput sc_input;
        sc_input.fault_buses = {from_sparse, {0, 0, 1}};
        sc_input.source = {u_ref};
        sc_input.faults = {{.y_fault = DoubleComplex{1.0, -1.0},
                            .fault_type = FaultType::three_phase,
                            .fault_phase = FaultPhase::abc}};
        CHECK_THROWS_AS(solver.run_short_circuit(sc_input, log, CalculationMethod::iec60909, y_bus),
                        ReducedPhaseNotSupported);
    }
}

} // namespace power_grid_model::math_solver
