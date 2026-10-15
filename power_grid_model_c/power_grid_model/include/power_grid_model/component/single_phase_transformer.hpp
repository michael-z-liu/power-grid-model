// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "branch.hpp"
#include "edge.hpp"
#include "line_utils.hpp"
#include "transformer_utils.hpp"

#include "../auxiliary/input.hpp"
#include "../auxiliary/output.hpp"
#include "../auxiliary/update.hpp"
#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/enum.hpp"
#include "../common/exception.hpp"
#include "../common/three_phase_tensor.hpp"
#include "component.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <format>
#include <string_view>
#include <utility>

namespace power_grid_model {

// single-phase two-winding transformer: each winding is connected between two phase slots of its node (ab, ac, bc)
// or between one phase slot and earth (a, b, c), the latter directly or through an earthing impedance
// the two coupled windings are stamped onto the slots of both nodes (C^T Y_prim C, in physical units) and scaled to
// per unit with the phase-to-earth base voltage u_rated / sqrt3 of each node
class SinglePhaseTransformer : public Branch {
  private:
    using Matrix2 = Eigen::Matrix<DoubleComplex, 2, 2>;
    // terminals 0, 1, 2: from slots a, b, c; 3, 4, 5: to slots a, b, c; 6, 7: earthed terminal of the from and to
    // winding behind an earthing impedance
    using Matrix8 = Eigen::Matrix<DoubleComplex, 8, 8>;
    static constexpr Idx earth{-1};
    static constexpr Idx to_offset{3};
    static constexpr Idx electrode_from{6};
    static constexpr Idx electrode_to{7};

    // winding of one side: the phase slots of the terminals and the angle of the winding voltage of a balanced
    // positive-sequence set (U_a = 1, U_b = a^2, U_c = a)
    struct Winding {
        PhaseMask slots;
        Idx positive{};      // slot 0, 1, 2
        Idx negative{earth}; // slot 0, 1, 2, or earth for a phase-to-earth winding
        double angle{};

        constexpr bool phase_to_phase() const { return negative != earth; }
    };

    // terminals of a winding in the index space of Matrix8; earth means the true earth (eliminated)
    struct Terminals {
        Idx positive{};
        Idx negative{};
    };

  public:
    using InputType = SinglePhaseTransformerInput;
    using UpdateType = TransformerUpdate;
    static constexpr char const* name = "single_phase_transformer";

    explicit SinglePhaseTransformer(SinglePhaseTransformerInput const& transformer_input, double u1_rated,
                                    double u2_rated)
        : Branch{transformer_input},
          u1_{transformer_input.u1},
          u2_{transformer_input.u2},
          sn_{transformer_input.sn},
          tap_size_{is_nan(transformer_input.tap_size) ? 0.0 : transformer_input.tap_size},
          uk_{transformer_input.uk},
          pk_{transformer_input.pk},
          i0_{transformer_input.i0},
          p0_{transformer_input.p0},
          from_{winding_of(transformer_input.from_connection)},
          to_{winding_of(transformer_input.to_connection)},
          clock_{checked_clock(id(), transformer_input.clock)},
          tap_side_{is_nan(transformer_input.tap_side) ? BranchSide::from : transformer_input.tap_side},
          tap_min_{transformer_input.tap_min},
          tap_max_{transformer_input.tap_max},
          tap_nom_{zero_if_nan(transformer_input.tap_nom)},
          // the tap dependent short circuit parameters at a limit only apply when that limit is provided
          uk_min_{is_nan(transformer_input.uk_min) || tap_min_ == na_IntS ? uk_ : transformer_input.uk_min},
          uk_max_{is_nan(transformer_input.uk_max) || tap_max_ == na_IntS ? uk_ : transformer_input.uk_max},
          pk_min_{is_nan(transformer_input.pk_min) || tap_min_ == na_IntS ? pk_ : transformer_input.pk_min},
          pk_max_{is_nan(transformer_input.pk_max) || tap_max_ == na_IntS ? pk_ : transformer_input.pk_max},
          base_u_from_{u1_rated / sqrt3},
          base_u_to_{u2_rated / sqrt3},
          base_i_from_{base_power_3p / u1_rated / sqrt3},
          base_i_to_{base_power_3p / u2_rated / sqrt3},
          z_grounding_from_{checked_grounding(id(), from_, transformer_input.r_grounding_from,
                                              transformer_input.x_grounding_from, "from")},
          z_grounding_to_{
              checked_grounding(id(), to_, transformer_input.r_grounding_to, transformer_input.x_grounding_to, "to")} {
        // not provided means the nominal tap
        tap_pos_ = tap_limit(transformer_input.tap_pos == na_IntS ? tap_nom_ : transformer_input.tap_pos);
    }

    // override getter
    double base_i_from() const final { return base_i_from_; }
    double base_i_to() const final { return base_i_to_; }
    // fallback only: get_output below computes the loading from the complex power of each side
    double loading(double max_s, double /* max_i */) const final { return max_s / sn_; }
    // phase shift is theta_from - theta_to: the angle of the to winding voltage minus that of the from winding
    // voltage, plus pi for reversed polarity (clock 6)
    double phase_shift() const final { return to_.angle - from_.angle + clock_ * deg_30; }
    BranchPhaseIncidence phase_incidence() const final {
        return {.from = from_.slots, .to = to_.slots, .coupling = PhaseCoupling::winding};
    }

    // hides Edge::get_output: the loading is the larger apparent power of the two sides, each from the complex sum of
    // its per-phase powers; the sum of the per-phase apparent powers would overstate the loading of a phase-to-phase
    // winding by 2 / sqrt3
    template <symmetry_tag sym> BranchOutput<sym> get_output(BranchSolverOutput<sym> const& solver_output) const {
        BranchOutput<sym> output = Edge::template get_output<sym>(solver_output);
        output.loading = std::max(side_power(solver_output.s_f), side_power(solver_output.s_t)) * base_power<sym> / sn_;
        return output;
    }

    // getters
    constexpr IntS tap_pos() const { return tap_pos_; }
    constexpr BranchSide tap_side() const { return tap_side_; }
    // a tap limit that is not provided follows the tap position, so it never limits it
    constexpr IntS tap_min() const { return tap_min_ == na_IntS ? std::min(tap_pos_, tap_nom_) : tap_min_; }
    constexpr IntS tap_max() const { return tap_max_ == na_IntS ? std::max(tap_pos_, tap_nom_) : tap_max_; }
    constexpr IntS tap_nom() const { return tap_nom_; }
    constexpr IntS clock() const { return clock_; }

    // setter
    constexpr bool set_tap(IntS new_tap) {
        if (new_tap == na_IntS || new_tap == tap_pos_) {
            return false;
        }
        tap_pos_ = tap_limit(new_tap);
        return true;
    }

    // update for transformer, hide default update for branch
    UpdateChange update(TransformerUpdate const& update_data) {
        assert(update_data.id == this->id() || is_nan(update_data.id));
        bool const topo_changed = set_status(update_data.from_status, update_data.to_status);
        bool const param_changed = set_tap(update_data.tap_pos) || topo_changed;
        return {.topo = topo_changed, .param = param_changed};
    }

    TransformerUpdate inverse(TransformerUpdate update_data) const {
        assert(update_data.id == this->id() || is_nan(update_data.id));

        update_data = Branch::inverse(update_data);
        set_if_not_nan(update_data.tap_pos, tap_pos_);

        return update_data;
    }

  private:
    // transformer parameter
    double u1_;
    double u2_;
    double sn_;
    double tap_size_;
    double uk_;
    double pk_;
    double i0_;
    double p0_;
    Winding from_;
    Winding to_;
    IntS clock_;
    BranchSide tap_side_;
    IntS tap_pos_{};
    IntS tap_min_; // na_IntS: not provided, see tap_min()
    IntS tap_max_; // na_IntS: not provided, see tap_max()
    IntS tap_nom_;
    double uk_min_;
    double uk_max_;
    double pk_min_;
    double pk_max_;

    // calculation parameter
    double base_u_from_;
    double base_u_to_;
    double base_i_from_;
    double base_i_to_;
    DoubleComplex z_grounding_from_; // ohm
    DoubleComplex z_grounding_to_;   // ohm

    static constexpr IntS zero_if_nan(IntS value) { return value == na_IntS ? IntS{0} : value; }

    static double side_power(ComplexValue<symmetric_t> const& s) { return cabs(s); }
    static double side_power(ComplexValue<asymmetric_t> const& s) { return cabs(s.sum()); }

    // the first letter is the positive terminal; angles of U_a - U_b etc. for U_a = 1, U_b = a^2, U_c = a
    static Winding winding_of(WindingConnection1P connection) {
        using enum WindingConnection1P;
        switch (connection) {
        case a:
            return {.slots = phase_mask_a, .positive = 0, .negative = earth, .angle = 0.0};
        case b:
            return {.slots = phase_mask_b, .positive = 1, .negative = earth, .angle = -deg_120};
        case c:
            return {.slots = phase_mask_c, .positive = 2, .negative = earth, .angle = deg_120};
        case ab:
            return {.slots = phase_mask_ab, .positive = 0, .negative = 1, .angle = deg_30};
        case ac:
            return {.slots = phase_mask_ac, .positive = 0, .negative = 2, .angle = -deg_30};
        case bc:
            return {.slots = phase_mask_bc, .positive = 1, .negative = 2, .angle = -3.0 * deg_30};
        default:
            throw MissingCaseForEnumError{"single_phase_transformer connection", connection};
        }
    }

    // not provided means 0; 0 and 6 (modulo 12) are the only valid values, 6 reverses the polarity of the to winding
    static IntS checked_clock(ID transformer_id, IntS clock) {
        if (clock == na_IntS) {
            return IntS{0};
        }
        auto const mapped = map_to_cyclic_range(clock, IntS{12});
        if (mapped != IntS{0} && mapped != IntS{6}) {
            throw InvalidTransformerClock{transformer_id, clock};
        }
        return mapped;
    }

    // earthing impedance in ohm, not provided means solidly earthed; only a phase-to-earth winding has an earthed
    // terminal, so the earthing impedance of a phase-to-phase winding must be left empty
    static DoubleComplex checked_grounding(ID transformer_id, Winding const& winding, double r, double x,
                                           std::string_view side) {
        if (winding.phase_to_phase() && (!is_nan(r) || !is_nan(x))) {
            throw InvalidPhaseConnection{
                name, transformer_id, phase_letters(winding.slots),
                std::format("{}_grounding_{} refers to a phase-to-phase winding and must be left empty",
                            is_nan(r) ? "x" : "r", side)};
        }
        return {is_nan(r) ? 0.0 : r, is_nan(x) ? 0.0 : x};
    }

    // limit the tap position as transformer does; a limit that is not provided takes the value of tap_min() or
    // tap_max() for the new tap position, so only the provided limits can limit it
    constexpr IntS tap_limit(IntS new_tap) const {
        IntS const lower = tap_min_ == na_IntS ? std::min(new_tap, tap_nom_) : tap_min_;
        IntS const upper = tap_max_ == na_IntS ? std::max(new_tap, tap_nom_) : tap_max_;
        new_tap = std::min(new_tap, std::max(upper, lower));
        new_tap = std::max(new_tap, std::min(upper, lower));
        return new_tap;
    }

    // as transformer: a higher tap position raises the voltage of the tapped winding when tap_max > tap_min
    constexpr IntS tap_direction() const { return tap_max() > tap_min() ? IntS{1} : IntS{-1}; }

    // terminals of a winding at the side with the given slot offset and earthing electrode
    static Terminals terminals_of(Winding const& winding, Idx offset, Idx electrode, DoubleComplex const& z_grounding) {
        if (winding.phase_to_phase()) {
            return {.positive = offset + winding.positive, .negative = offset + winding.negative};
        }
        return {.positive = offset + winding.positive, .negative = z_grounding == DoubleComplex{} ? earth : electrode};
    }

    static void add_entry(Matrix8& y, Idx row, Idx col, DoubleComplex const& value) {
        if (row != earth && col != earth) {
            y(row, col) += value;
        }
    }

    // MNA stamp of the admittance value between winding row and winding col (row of C^T Y_prim C)
    static void stamp(Matrix8& y, Terminals const& row, Terminals const& col, DoubleComplex const& value) {
        add_entry(y, row.positive, col.positive, value);
        add_entry(y, row.positive, col.negative, -value);
        add_entry(y, row.negative, col.positive, -value);
        add_entry(y, row.negative, col.negative, value);
    }

    // primitive admittance of the two windings in siemens (winding currents from winding voltages), including the
    // off-nominal tap; with one side open, the open winding carries no current and the connected winding sees the
    // open-end shunt 1/2 y_m + (1/y_s + 2/y_m)^-1
    Matrix2 primitive_admittance() const {
        double u1 = u1_;
        double u2 = u2_;
        double const tap_offset = tap_direction() * (tap_pos_ - tap_nom_) * tap_size_;
        if (tap_side_ == BranchSide::from) {
            u1 += tap_offset;
        } else {
            u2 += tap_offset;
        }
        double const uk = tap_adjust_impedance(tap_pos_, tap_min(), tap_max(), tap_nom_, uk_, uk_min_, uk_max_);
        double const pk = tap_adjust_impedance(tap_pos_, tap_min(), tap_max(), tap_nom_, pk_, pk_min_, pk_max_);

        // series impedance referred to side 2, as for transformer
        double const z_series_abs = std::abs(uk) * u2 * u2 / sn_;
        double const r_series = pk * u2 * u2 / sn_ / sn_;
        double const x_series_squared = z_series_abs * z_series_abs - r_series * r_series;
        double const x_series = (uk >= 0 ? 1.0 : -1.0) * (x_series_squared > 0.0 ? std::sqrt(x_series_squared) : 0.0);
        DoubleComplex const y_series = 1.0 / DoubleComplex{r_series, x_series};

        // magnetising admittance referred to side 2
        double const y_shunt_abs = i0_ * sn_ / u2 / u2;
        double const g_shunt = p0_ / u2 / u2;
        double const b_shunt_squared = y_shunt_abs * y_shunt_abs - g_shunt * g_shunt;
        DoubleComplex const y_shunt{g_shunt, b_shunt_squared > 0.0 ? -std::sqrt(b_shunt_squared) : 0.0};

        double const ratio = u1 / u2;
        Matrix2 y_prim = Matrix2::Zero();
        if (edge_status()) {
            y_prim(0, 0) = (y_series + 0.5 * y_shunt) / ratio / ratio;
            y_prim(0, 1) = -y_series / ratio;
            y_prim(1, 0) = -y_series / ratio;
            y_prim(1, 1) = y_series + 0.5 * y_shunt;
            return y_prim;
        }
        // the tolerance applies to the magnetising admittance in per unit of the to side, as in Edge::calc_param_y_sym
        if (cabs(y_shunt) * base_u_to_ * base_u_to_ / base_power_1p < numerical_tolerance) {
            return y_prim;
        }
        DoubleComplex const open_end = 0.5 * y_shunt + 1.0 / (1.0 / y_series + 2.0 / y_shunt);
        if (from_status()) {
            y_prim(0, 0) = open_end / ratio / ratio;
        }
        if (to_status()) {
            y_prim(1, 1) = open_end;
        }
        return y_prim;
    }

    // nodal admittance in siemens of the from and to slots: the windings stamped onto their terminals, the earthing
    // impedances added and their terminals eliminated
    Matrix8 nodal_admittance() const {
        Matrix2 const y_prim = primitive_admittance();
        Terminals const from_terminals = terminals_of(from_, 0, electrode_from, z_grounding_from_);
        Terminals to_terminals = terminals_of(to_, to_offset, electrode_to, z_grounding_to_);
        if (clock_ == IntS{6}) {
            std::swap(to_terminals.positive, to_terminals.negative);
        }

        Matrix8 y = Matrix8::Zero();
        stamp(y, from_terminals, from_terminals, y_prim(0, 0));
        stamp(y, from_terminals, to_terminals, y_prim(0, 1));
        stamp(y, to_terminals, from_terminals, y_prim(1, 0));
        stamp(y, to_terminals, to_terminals, y_prim(1, 1));

        // a non-zero earthing impedance is only possible on a phase-to-earth winding, whose earthed terminal is then
        // the electrode
        if (z_grounding_from_ != DoubleComplex{}) {
            y(electrode_from, electrode_from) += 1.0 / z_grounding_from_;
            kron_eliminate(y, electrode_from);
        }
        if (z_grounding_to_ != DoubleComplex{}) {
            y(electrode_to, electrode_to) += 1.0 / z_grounding_to_;
            kron_eliminate(y, electrode_to);
        }
        return y;
    }

    // a 3x3 block of the nodal admittance in per unit: Y_pu(s, t) = Y(s, t) U_b,s U_b,t / S_b,1p
    static ComplexTensor<asymmetric_t> per_unit_block(Matrix8 const& y, Idx row_offset, Idx col_offset,
                                                      double base_u_row, double base_u_col) {
        double const scale = base_u_row * base_u_col / base_power_1p;
        return ComplexTensor<asymmetric_t>{(y.block<3, 3>(row_offset, col_offset) * scale).array()};
    }

    // the common mode of a phase-to-phase winding has no earth reference: add a small admittance to earth of
    // y_low / 2 on each entry of its two slots, the analogue of the zero-sequence low susceptance of transformer
    void add_low_susceptance(ComplexTensor<asymmetric_t>& y_block, Winding const& winding) const {
        if (!winding.phase_to_phase()) {
            return;
        }
        DoubleComplex const half_low_susceptance{0.0,
                                                 -0.5 * transformer_low_susceptance_ratio * sn_ / base_power_1p / uk_};
        y_block(winding.positive, winding.positive) += half_low_susceptance;
        y_block(winding.positive, winding.negative) += half_low_susceptance;
        y_block(winding.negative, winding.positive) += half_low_susceptance;
        y_block(winding.negative, winding.negative) += half_low_susceptance;
    }

    BranchCalcParam<symmetric_t> sym_calc_param() const final {
        throw ReducedPhaseNotSupported{"Symmetric calculation"};
    }

    BranchCalcParam<asymmetric_t> asym_calc_param() const final {
        Matrix8 const y = nodal_admittance();
        BranchCalcParam<asymmetric_t> param{};
        param.yff() = per_unit_block(y, 0, 0, base_u_from_, base_u_from_);
        param.yft() = per_unit_block(y, 0, to_offset, base_u_from_, base_u_to_);
        param.ytf() = per_unit_block(y, to_offset, 0, base_u_to_, base_u_from_);
        param.ytt() = per_unit_block(y, to_offset, to_offset, base_u_to_, base_u_to_);
        if (from_status()) {
            add_low_susceptance(param.yff(), from_);
        }
        if (to_status()) {
            add_low_susceptance(param.ytt(), to_);
        }
        // every slot outside the windings is exactly zero, by assignment
        mask_absent_phases(param.yff(), from_.slots, from_.slots);
        mask_absent_phases(param.yft(), from_.slots, to_.slots);
        mask_absent_phases(param.ytf(), to_.slots, from_.slots);
        mask_absent_phases(param.ytt(), to_.slots, to_.slots);
        return param;
    }
};

static_assert(transformer_c<SinglePhaseTransformer>);

} // namespace power_grid_model
