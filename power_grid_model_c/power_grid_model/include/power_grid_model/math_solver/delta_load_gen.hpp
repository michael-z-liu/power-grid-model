// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Delta (phase-to-phase) connected load_gen for the asymmetric power-flow solvers
//
// The three entries of a delta load_gen are the branches ab, bc and ca (see delta_branches). Branch k = (p, q) has the
// voltage U_k = U_p - U_q and the voltage-dependent power
//     S_k(U) = S_k * (|U_k| / U_nom)^n,    n = 0 (const_pq), 1 (const_i), 2 (const_y),
// where U_nom is the nominal branch voltage in per unit of the phase-to-earth base: sqrt(3) at a three-phase bus and 2
// at a split-phase bus. The branch current I_k = conj(S_k(U) / U_k) flows into phase p and out of phase q, so the phase
// current injection is I = D^T I_delta with D = [[1, -1, 0], [0, 1, -1], [-1, 0, 1]].
//
// With alpha_p = U_p / U_k and alpha_q = -U_q / U_k (alpha_p + alpha_q = 1), the branch injects the phase powers
// g_p = S_k(U) * alpha_p and g_q = S_k(U) * alpha_q. With beta = S_k(U) * alpha_p * alpha_q and sigma = [[1, -1],
// [-1, 1]] over r, s in {p, q}, the derivatives are
//     dg_r/dtheta_s     = j * sigma_rs * beta - n * g_r * Im(alpha_s)
//     V_s * dg_r/dV_s   = sigma_rs * beta + n * g_r * Re(alpha_s)
// All of them stay within the bus, so they only add to the diagonal block of the bus.
//
// An entry is inactive if it is exactly zero or if one of its two phases is absent at the bus. Inactive entries are
// skipped before any division, so they contribute exactly zero and never produce NaN, even if U_k is zero. An active
// entry with U_k = 0 (two present phases at the same voltage, not physical in a power flow) gives a non-finite result.

#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/enum.hpp"
#include "../common/exception.hpp"
#include "../common/three_phase_tensor.hpp"

#include <array>
#include <cassert>
#include <complex>
#include <cstddef>

namespace power_grid_model::math_solver::delta_load_gen {

// nominal branch voltage in per unit of the phase-to-earth voltage base of the bus
inline double nominal_voltage(SlotPattern pattern) {
    using enum SlotPattern;
    switch (pattern) {
    case three_phase:
        return sqrt3; // the phase-to-phase voltage u_rated
    case split_phase:
        return 2.0; // the voltage across both legs of the centre-tapped winding
    default:
        throw MissingCaseForEnumError{"Delta load_gen nominal voltage", pattern};
    }
}

// exponent n of the voltage dependence S(U) = S * (|U| / U_nom)^n
inline int voltage_exponent(LoadGenType type) {
    using enum LoadGenType;
    switch (type) {
    case const_pq:
        return 0;
    case const_i:
        return 1;
    case const_y:
        return 2;
    default:
        throw MissingCaseForEnumError{"Delta load_gen voltage exponent", type};
    }
}

// a delta entry contributes iff it is not exactly zero (NaN counts as active) and both of its phases are present
inline bool is_active(DoubleComplex const& s, PhaseMask bus_phases, DeltaBranch const& branch) {
    return s != DoubleComplex{} && bus_phases.contains(branch.phases);
}

// power of a branch at the branch voltage u_k: s * (|u_k| / u_nom)^n
inline DoubleComplex branch_power(DoubleComplex const& s, int n, DoubleComplex const& u_k, double u_nom) {
    assert(n >= 0 && n <= 2);
    switch (n) {
    case 0:
        return s;
    case 1:
        return s * (std::abs(u_k) / u_nom);
    default:
        return s * (std::norm(u_k) / (u_nom * u_nom));
    }
}

// phase current injection I = D^T I_delta of a delta load_gen with the entries s_delta (ab, bc, ca)
inline ComplexValue<asymmetric_t> phase_current(ComplexValue<asymmetric_t> const& s_delta, LoadGenType type,
                                                ComplexValue<asymmetric_t> const& u, PhaseMask bus_phases,
                                                double u_nom) {
    int const n = voltage_exponent(type);
    ComplexValue<asymmetric_t> i{}; // all slots exactly zero
    for (std::size_t k = 0; k != delta_branches.size(); ++k) {
        DeltaBranch const& branch = delta_branches[k];
        DoubleComplex const s = s_delta(static_cast<Idx>(k));
        if (!is_active(s, bus_phases, branch)) {
            continue;
        }
        DoubleComplex const u_k = u(branch.from) - u(branch.to);
        DoubleComplex const i_k = std::conj(branch_power(s, n, u_k, u_nom) / u_k);
        i(branch.from) += i_k;
        i(branch.to) -= i_k;
    }
    return i;
}

// phase power injection s = U o conj(I) of a delta load_gen; the slots without an active entry are exactly zero
inline ComplexValue<asymmetric_t> phase_power(ComplexValue<asymmetric_t> const& s_delta, LoadGenType type,
                                              ComplexValue<asymmetric_t> const& u, PhaseMask bus_phases, double u_nom) {
    ComplexValue<asymmetric_t> const i = phase_current(s_delta, type, u, bus_phases, u_nom);
    PhaseMask loaded{};
    for (std::size_t k = 0; k != delta_branches.size(); ++k) {
        if (is_active(s_delta(static_cast<Idx>(k)), bus_phases, delta_branches[k])) {
            loaded = loaded | delta_branches[k].phases;
        }
    }
    ComplexValue<asymmetric_t> s{}; // all slots exactly zero
    for (Idx phase = 0; phase != 3; ++phase) {
        if (loaded.has(phase)) {
            s(phase) = u(phase) * std::conj(i(phase));
        }
    }
    return s;
}

// phase power injection and its derivatives with respect to the voltage angles (d/dtheta) and the scaled voltage
// magnitudes (V * d/dV); row r is the injection at phase r, column s the voltage of phase s
struct PowerAndJacobian {
    ComplexValue<asymmetric_t> s{};
    RealTensor<asymmetric_t> dp_dtheta{};
    RealTensor<asymmetric_t> dp_dv{};
    RealTensor<asymmetric_t> dq_dtheta{};
    RealTensor<asymmetric_t> dq_dv{};
};

inline PowerAndJacobian power_and_jacobian(ComplexValue<asymmetric_t> const& s_delta, LoadGenType type,
                                           ComplexValue<asymmetric_t> const& u, PhaseMask bus_phases, double u_nom) {
    int const n = voltage_exponent(type);
    auto const n_real = static_cast<double>(n);
    PowerAndJacobian result{}; // all entries exactly zero
    for (std::size_t k = 0; k != delta_branches.size(); ++k) {
        DeltaBranch const& branch = delta_branches[k];
        DoubleComplex const s = s_delta(static_cast<Idx>(k));
        if (!is_active(s, bus_phases, branch)) {
            continue;
        }
        DoubleComplex const u_k = u(branch.from) - u(branch.to);
        DoubleComplex const s_u = branch_power(s, n, u_k, u_nom);
        std::array<Idx, 2> const idx{branch.from, branch.to};
        std::array<DoubleComplex, 2> const alpha{u(branch.from) / u_k, -u(branch.to) / u_k};
        std::array<DoubleComplex, 2> const g{s_u * alpha[0], s_u * alpha[1]};
        DoubleComplex const beta = s_u * alpha[0] * alpha[1];
        for (std::size_t r = 0; r != 2; ++r) {
            result.s(idx[r]) += g[r];
            for (std::size_t c = 0; c != 2; ++c) {
                double const sigma = r == c ? 1.0 : -1.0;
                result.dp_dtheta(idx[r], idx[c]) += -sigma * beta.imag() - n_real * g[r].real() * alpha[c].imag();
                result.dq_dtheta(idx[r], idx[c]) += sigma * beta.real() - n_real * g[r].imag() * alpha[c].imag();
                result.dp_dv(idx[r], idx[c]) += sigma * beta.real() + n_real * g[r].real() * alpha[c].real();
                result.dq_dv(idx[r], idx[c]) += sigma * beta.imag() + n_real * g[r].imag() * alpha[c].real();
            }
        }
    }
    return result;
}

// constant admittance of a delta load_gen at the nominal branch voltage, Y = -D^T diag(conj(s_delta) / U_nom^2) D, so
// that Y U = -phase_current(s_delta, const_y, U); the injected current is -Y U
inline ComplexTensor<asymmetric_t> linear_admittance(ComplexValue<asymmetric_t> const& s_delta, PhaseMask bus_phases,
                                                     double u_nom) {
    ComplexTensor<asymmetric_t> y{}; // all entries exactly zero
    for (std::size_t k = 0; k != delta_branches.size(); ++k) {
        DeltaBranch const& branch = delta_branches[k];
        DoubleComplex const s = s_delta(static_cast<Idx>(k));
        if (!is_active(s, bus_phases, branch)) {
            continue;
        }
        DoubleComplex const y_k = -std::conj(s) / (u_nom * u_nom);
        y(branch.from, branch.from) += y_k;
        y(branch.to, branch.to) += y_k;
        y(branch.from, branch.to) -= y_k;
        y(branch.to, branch.from) -= y_k;
    }
    return y;
}

} // namespace power_grid_model::math_solver::delta_load_gen
