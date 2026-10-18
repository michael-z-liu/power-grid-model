// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

// Check if the name means anything
#pragma once

/*
Iterative Power Flow

Description:
    The algorithm similar to jacobi way of solving linear equations.
    Only I_inj is calculated fresh on each iteration based on latest values of U.
    Linear equation here: I_inj = YU

Prefactorization:
    If the Y bus matrix does not change, then there is no need for factorizing it again to solve linear equations.
    Hence it is done only once in the first iteration and same result is used in subsequent iterations.
    Same factorization is also used in subsequent batches

Steps:
    Initialize U with averaged u_ref, ie source voltage and phase shifts accounted
    Initialize solver
    while maximum deviation > error tolerance
        Calculate I_inj with U of previous iteration as per load/gen types.
        Solve YU = I_inj using prefactorization.
        Find maximum deviation in voltage buses U
        Update U
    Calculate output values from U result

    Initialize solver:
        Source admittance is not included in Y bus matrix here. Include that to complete the Y bus matrix.
        Invalidate prefactorization if parameters change, ie y bus values changes

    Calculating Injected current:
        Initialize I_inj = 0
        For each bus i
            For source on bus i, I_inj_i = y_ref * u_uref
            For Loads on bus i:
                If type is constant PQ: I_inj_i += conj(S_inj_j/U_i)
                If type is constant impedance: I_inj_i += conj((S_inj_j * abs(U_i)^2) / U_i) = conj((S_inj_j) * U_i
                If type is constant current: I_inj_i += conj(S_inj_j*abs(U_i)/U_i)
                If the load is delta (phase-to-phase) connected: I_inj_i += D^T I_delta, with the branch currents
                    I_delta_k = conj(S_k(U) / U_k) of the type-dependent branch powers (see delta_load_gen.hpp)

    Constant impedance loads in a network with an absent phase (a bus with a one- or two-phase slot pattern):
        The current of a constant impedance load is linear in U, I_inj = -Y_load U, with Y_load = -conj(S_inj) on the
        diagonal (wye) or the branch admittances of a delta load (see add_linear_loads). In such a network these loads
        are moved from I_inj into the matrix, which is then Y bus + y_ref + Y_load; the solution is the same.
        Explicitly, a constant impedance load to earth is the only earth reference of a floating common mode (e.g. the
        ab secondary of an ab-ab single phase transformer, held to earth only by a low admittance), and the iteration
        diverges there with a gain of about Y_load / (low admittance) per iteration; in the matrix, the load anchors it.
        - The matrix is refactorized whenever the constant impedance loads of the input change, so that it always
          holds those of the current input (deterministic in a batch, whatever the order of the scenarios).
        - A floating common mode with only constant power or constant current loads still diverges.
        - Linear current (one iteration) then treats constant impedance loads exactly in such a network.
        - Networks without an absent phase, including a network flagged reduced-phase with all buses abc, keep the
          explicit treatment, bit for bit.

Nomenclature:
    I_inj : Injected current
    S_inj : Injected power
    Y : Y bus matrix
    U : Bus Voltage
    u_ref : reference voltage for source
    y_ref : Source admittance


*/

#include "common_solver_functions.hpp"
#include "delta_load_gen.hpp"
#include "iterative_pf_solver.hpp"
#include "sparse_lu_solver.hpp"
#include "y_bus.hpp"

#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/counting_iterator.hpp"
#include "../common/enum.hpp"
#include "../common/exception.hpp"
#include "../common/three_phase_tensor.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

namespace power_grid_model::math_solver {

// hide implementation in inside namespace
namespace iterative_current_pf {

// solver
template <symmetry_tag sym_type>
class IterativeCurrentPFSolver : public IterativePFSolver<sym_type, IterativeCurrentPFSolver<sym_type>> {
  public:
    using sym = sym_type;

    using SparseSolverType = SparseLUSolver<ComplexTensor<sym>, ComplexValue<sym>, ComplexValue<sym>>;
    using BlockPermArray = SparseSolverType::BlockPermArray;

    static constexpr auto is_iterative = true;

    IterativeCurrentPFSolver(YBus<sym> const& y_bus, MathModelTopology const& topo)
        : IterativePFSolver<sym, IterativeCurrentPFSolver>{y_bus, topo},
          rhs_u_(y_bus.size()),
          sparse_solver_{y_bus.row_indptr_lu(), y_bus.col_indices_lu(), y_bus.lu_diag()},
          absent_phase_{has_absent_phase(topo)} {}

    // Add source admittance to Y bus and set variable for prepared y bus to true
    void initialize_derived_solver(YBus<sym> const& y_bus, PowerFlowInput<sym> const& input,
                                   SolverOutput<sym>& output) {
        this->make_average_source_start(input, output.u);
        detail::apply_slot_patterns<sym>(y_bus.math_topology(), output.u);
        detail::mask_absent_phase_values<sym>(y_bus.math_topology(), output.u);

        auto const& sources_per_bus = this->sources_per_bus_.get();
        IdxVector const& bus_entry = y_bus.lu_diag();
        // at an absent-phase network, the matrix holds the const_y load_gens of this input, so it is refactorized
        // whenever they change
        if (absent_phase_) {
            update_const_y_admittance(y_bus.math_topology(), input);
        }
        // if Y bus is not up to date
        // re-build matrix and prefactorize Build y bus data with source admittance
        if (parameters_changed_) {
            ComplexTensorVector<sym> mat_data(y_bus.nnz_lu());
            detail::copy_y_bus<sym>(y_bus, mat_data);

            for (auto const& [bus_number, sources] : enumerated_zip_sequence(sources_per_bus)) {
                Idx const data_sequence = bus_entry[bus_number];
                for (auto source_number : sources) {
                    // YBus_diag += Y_source // NOSONAR
                    mat_data[data_sequence] +=
                        y_bus.math_model_param().source_param[source_number].template y_ref<sym>();
                }
            }
            if (absent_phase_) {
                // YBus_diag += Y_load of the const_y load_gens // NOSONAR(S125)
                for (auto const& [bus_number, load_gens] : enumerated_zip_sequence(this->load_gens_per_bus_.get())) {
                    if (has_const_y_load(load_gens)) {
                        mat_data[bus_entry[bus_number]] += const_y_admittance_[bus_number];
                    }
                }
            }
            // prefactorize
            BlockPermArray perm(this->n_bus_);
            sparse_solver_.prefactorize(mat_data, perm);
            // move pre-factorized version into shared ptr
            mat_data_ = std::make_shared<ComplexTensorVector<sym> const>(std::move(mat_data));
            perm_ = std::make_shared<BlockPermArray const>(std::move(perm));
        }
        parameters_changed_ = false;
    }

    // Prepare matrix calculates injected current, i.e., RHS of solver for each iteration.
    void prepare_matrix_and_rhs(YBus<sym> const& y_bus, PowerFlowInput<sym> const& input,
                                ComplexValueVector<sym> const& u) {
        std::vector<LoadGenType> const& load_gen_type = this->load_gen_type_.get();

        // set rhs to zero for iteration start
        std::ranges::fill(rhs_u_, ComplexValue<sym>{0.0});

        auto const& sources_per_bus = this->sources_per_bus_.get();
        auto const& load_gens_per_bus = this->load_gens_per_bus_.get();
        // loop buses: i
        for (auto const& [bus_number, load_gens, sources] :
             enumerated_zip_sequence(load_gens_per_bus, sources_per_bus)) {
            add_loads(load_gens, bus_number, input, load_gen_type, u, y_bus.math_topology());
            add_sources(sources, bus_number, y_bus, input);
        }
        // a load on an absent phase is ignored; its current conj(s / u) is not finite at u = 0
        detail::mask_absent_phase_values<sym>(y_bus.math_topology(), rhs_u_);
    }

    // Solve the linear equations I_inj = YU
    // inplace
    void solve_matrix() { sparse_solver_.solve_with_prefactorized_matrix(*mat_data_, *perm_, rhs_u_, rhs_u_); }

    // Find maximum deviation in voltage among all buses
    double iterate_unknown(ComplexValueVector<sym>& u, double /*err_tol*/, bool /*cache_run*/) {
        double max_dev = 0.0;
        // loop all buses
        for (Idx bus_number = 0; bus_number != this->n_bus_; ++bus_number) {
            // Keep maximum deviation of all buses; a NaN deviation is kept
            max_dev = detail::fold_max_deviation<sym>(max_dev, cabs(rhs_u_[bus_number] - u[bus_number]));
            // assign updated values
            u[bus_number] = rhs_u_[bus_number];
        }
        return max_dev;
    }

    void parameters_changed(bool changed) { parameters_changed_ = parameters_changed_ || changed; }

  private:
    ComplexValueVector<sym> rhs_u_;
    std::shared_ptr<ComplexTensorVector<sym> const> mat_data_;
    // sparse solver
    SparseSolverType sparse_solver_;
    std::shared_ptr<BlockPermArray const> perm_;
    bool parameters_changed_ = true;
    // a bus has an absent phase (asymmetric only); then the const_y load_gens are part of the matrix instead of the
    // right-hand side
    bool absent_phase_;
    // per bus, the admittance of the const_y load_gens in the current matrix; empty unless absent_phase_
    ComplexTensorVector<sym> const_y_admittance_;

    // the gate of the const_y handling: not has_reduced_phase alone, which a math model with all buses abc may have
    static bool has_absent_phase(MathModelTopology const& topo) {
        if constexpr (is_symmetric_v<sym>) {
            return false;
        } else {
            return topo.has_reduced_phase &&
                   std::ranges::any_of(topo.bus_phase_mask, [](PhaseMask phases) { return !phases.all(); });
        }
    }

    // the const_y admittance per bus of this input; if any entry differs (exactly) from the one in the matrix, store it
    // and mark the matrix for refactorization, so that the matrix always holds the const_y load_gens of the current
    // input, whatever the order of the scenarios in a batch
    void update_const_y_admittance(MathModelTopology const& topo, PowerFlowInput<sym> const& input) {
        std::vector<LoadGenType> const& load_gen_type = this->load_gen_type_.get();
        auto const& load_gens_per_bus = this->load_gens_per_bus_.get();
        ComplexTensorVector<sym> y_load(static_cast<std::size_t>(this->n_bus_), ComplexTensor<sym>{});
        for (auto const& [bus_number, load_gens] : enumerated_zip_sequence(load_gens_per_bus)) {
            detail::add_const_y_loads(load_gens, bus_number, input, load_gen_type, y_load[bus_number], topo);
        }
        auto const exactly_equal = [](ComplexTensor<sym> const& x, ComplexTensor<sym> const& y) {
            if constexpr (is_symmetric_v<sym>) {
                return x == y;
            } else {
                return (x == y).all();
            }
        };
        if (!std::ranges::equal(y_load, const_y_admittance_, exactly_equal)) {
            const_y_admittance_ = std::move(y_load);
            parameters_changed_ = true;
        }
    }

    // a bus without const_y load_gens keeps its matrix entry untouched
    bool has_const_y_load(IdxRange const& load_gens) const {
        std::vector<LoadGenType> const& load_gen_type = this->load_gen_type_.get();
        return std::ranges::any_of(load_gens,
                                   [&load_gen_type](Idx load) { return load_gen_type[load] == LoadGenType::const_y; });
    }

    void add_loads(IdxRange const& load_gens, Idx bus_number, PowerFlowInput<sym> const& input,
                   std::vector<LoadGenType> const& load_gen_type, ComplexValueVector<sym> const& u,
                   [[maybe_unused]] MathModelTopology const& topo) {
        for (Idx const load_number : load_gens) {
            // load type
            LoadGenType const type = load_gen_type[load_number];
            // at an absent-phase network the const_y load_gens are in the matrix
            if (absent_phase_ && type == LoadGenType::const_y) {
                continue;
            }
            if constexpr (!is_symmetric_v<sym>) {
                if (topo.is_delta_load_gen(load_number)) {
                    // I_inj_i = D^T conj(S_k(U) / U_k) over the branches k of the delta load
                    rhs_u_[bus_number] += delta_load_gen::phase_current(
                        input.s_injection[load_number], type, u[bus_number], topo.bus_phase(bus_number),
                        delta_load_gen::nominal_voltage(topo.slot_pattern(bus_number)));
                    continue;
                }
            }
            switch (type) {
                using enum LoadGenType;

            case const_pq:
                // I_inj_i = conj(S_inj_j/U_i) for constant PQ type
                rhs_u_[bus_number] += conj(input.s_injection[load_number] / u[bus_number]);
                break;
            case const_y:
                // I_inj_i = conj((S_inj_j * abs(U_i)^2) / U_i) = conj((S_inj_j) * U_i for const impedance type
                rhs_u_[bus_number] += conj(input.s_injection[load_number]) * u[bus_number];
                break;
            case const_i:
                // I_inj_i = conj(S_inj_j*abs(U_i)/U_i) for const current type
                rhs_u_[bus_number] += conj(input.s_injection[load_number] * cabs(u[bus_number]) / u[bus_number]);
                break;
            default:
                throw MissingCaseForEnumError("Injection current calculation", type);
            }
        }
    }

    void add_sources(IdxRange const& sources, Idx bus_number, YBus<sym> const& y_bus,
                     PowerFlowInput<sym> const& input) {
        for (Idx const source_number : sources) {
            // I_inj_i += Y_source_j * U_ref_j // NOSONAR(S125)
            rhs_u_[bus_number] += dot(y_bus.math_model_param().source_param[source_number].template y_ref<sym>(),
                                      ComplexValue<sym>{input.source[source_number]});
        }
    }
};

} // namespace iterative_current_pf

using iterative_current_pf::IterativeCurrentPFSolver;

} // namespace power_grid_model::math_solver
