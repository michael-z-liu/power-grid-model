// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "branch.hpp"
#include "line_utils.hpp"

#include "../auxiliary/input.hpp"
#include "../auxiliary/update.hpp"
#include "../calculation_parameters.hpp"
#include "../common/common.hpp"
#include "../common/enum.hpp"
#include "../common/exception.hpp"
#include "../common/matrix_utils.hpp"
#include "../common/three_phase_tensor.hpp"

#include <array>
#include <complex>
#include <format>
#include <string_view>

namespace power_grid_model {

class AsymLine : public Branch {
  public:
    using InputType = AsymLineInput;
    using UpdateType = BranchUpdate;
    static constexpr char const* name = "asym_line";

    explicit AsymLine(AsymLineInput const& asym_line_input, double system_frequency, double u1, double u2)
        : Branch{asym_line_input},
          i_n_{asym_line_input.i_n},
          base_i_{base_power_3p / u1 / sqrt3},
          phases_{phase_connection_mask(asym_line_input.phases)} {
        if (cabs(u1 - u2) > numerical_tolerance) {
            throw ConflictVoltage{id(), from_node(), to_node(), u1, u2};
        }

        if (phases_.all()) {
            ComplexTensor<asymmetric_t> const c_matrix = compute_c_matrix_from_input(asym_line_input);
            ComplexTensor<asymmetric_t> const z_series = compute_z_series_from_input(asym_line_input);

            double const base_y = base_i_ / (u1 / sqrt3);

            y_series_abc_ = 1 / base_y * inv(z_series);
            y_shunt_abc_ = 1 / base_y * (2.0i * pi * system_frequency * c_matrix);
            return;
        }

        // one- or two-phase line: only the present phases (and the neutral, if any) have values, the absent rows and
        // columns of every matrix are exactly zero
        check_reduced_phase_input(asym_line_input, id(), phases_);
        ComplexTensor<asymmetric_t> c_matrix = compute_c_matrix_from_input(asym_line_input);
        mask_absent_phases(c_matrix, phases_, phases_);
        ComplexTensor<asymmetric_t> const z_series = compute_reduced_z_series_from_input(asym_line_input, phases_);

        double const base_y = base_i_ / (u1 / sqrt3);

        y_series_abc_ = 1 / base_y * masked_inv(z_series, phases_);
        y_shunt_abc_ = 1 / base_y * (2.0i * pi * system_frequency * c_matrix);
    }

    // override getter
    constexpr double base_i_from() const override { return base_i_; }
    constexpr double base_i_to() const override { return base_i_; }
    constexpr double loading(double /* max_s */, double max_i) const override { return max_i / i_n_; };
    constexpr double phase_shift() const override { return 0.0; }
    // for abc this equals the default abc -> abc, per phase
    BranchPhaseIncidence phase_incidence() const final {
        return {.from = phases_, .to = phases_, .coupling = PhaseCoupling::per_phase};
    }

  private:
    double i_n_{};
    double base_i_{};
    PhaseMask phases_{phase_mask_abc};
    ComplexTensor<asymmetric_t> y_series_abc_{};
    ComplexTensor<asymmetric_t> y_shunt_abc_{};

    // phases with a conductor; not provided means abc
    static PhaseMask phase_connection_mask(PhaseConnection phases) {
        if (is_nan(phases)) {
            return phase_mask_abc;
        }
        switch (phases) {
            using enum PhaseConnection;
        case abc:
            return phase_mask_abc;
        case a:
            return phase_mask_a;
        case b:
            return phase_mask_b;
        case c:
            return phase_mask_c;
        case ab:
            return phase_mask_ab;
        case ac:
            return phase_mask_ac;
        case bc:
            return phase_mask_bc;
        default:
            throw MissingCaseForEnumError{"asym_line phases", phases};
        }
    }

    // the neutral of a one- or two-phase line is given by its self impedance
    static bool reduced_phase_has_neutral(AsymLineInput const& asym_line_input) {
        return !is_nan(asym_line_input.r_nn) || !is_nan(asym_line_input.x_nn);
    }

    // an entry of the r, x or c matrix between two conductors: 0, 1, 2 for a, b, c and 3 for the neutral
    struct ConductorPairEntry {
        std::string_view name;
        double value{};
        Idx row{};
        Idx col{};
    };

    // an entry is required if both of its conductors are present and must be empty (nan) otherwise
    static void check_reduced_phase_input(AsymLineInput const& in, ID line_id, PhaseMask phases) {
        if (!is_nan(in.c0) || !is_nan(in.c1)) {
            throw InvalidPhaseConnection{name, line_id, phase_letters(phases),
                                         "c0 and c1 are only supported with all three phases; use the c matrix"};
        }
        bool const has_neutral = reduced_phase_has_neutral(in);
        auto const is_present = [phases, has_neutral](Idx conductor) {
            return conductor == 3 ? has_neutral : phases.has(conductor);
        };
        std::array const entries{
            ConductorPairEntry{.name = "r_aa", .value = in.r_aa, .row = 0, .col = 0},
            ConductorPairEntry{.name = "r_ba", .value = in.r_ba, .row = 1, .col = 0},
            ConductorPairEntry{.name = "r_bb", .value = in.r_bb, .row = 1, .col = 1},
            ConductorPairEntry{.name = "r_ca", .value = in.r_ca, .row = 2, .col = 0},
            ConductorPairEntry{.name = "r_cb", .value = in.r_cb, .row = 2, .col = 1},
            ConductorPairEntry{.name = "r_cc", .value = in.r_cc, .row = 2, .col = 2},
            ConductorPairEntry{.name = "r_na", .value = in.r_na, .row = 3, .col = 0},
            ConductorPairEntry{.name = "r_nb", .value = in.r_nb, .row = 3, .col = 1},
            ConductorPairEntry{.name = "r_nc", .value = in.r_nc, .row = 3, .col = 2},
            ConductorPairEntry{.name = "r_nn", .value = in.r_nn, .row = 3, .col = 3},
            ConductorPairEntry{.name = "x_aa", .value = in.x_aa, .row = 0, .col = 0},
            ConductorPairEntry{.name = "x_ba", .value = in.x_ba, .row = 1, .col = 0},
            ConductorPairEntry{.name = "x_bb", .value = in.x_bb, .row = 1, .col = 1},
            ConductorPairEntry{.name = "x_ca", .value = in.x_ca, .row = 2, .col = 0},
            ConductorPairEntry{.name = "x_cb", .value = in.x_cb, .row = 2, .col = 1},
            ConductorPairEntry{.name = "x_cc", .value = in.x_cc, .row = 2, .col = 2},
            ConductorPairEntry{.name = "x_na", .value = in.x_na, .row = 3, .col = 0},
            ConductorPairEntry{.name = "x_nb", .value = in.x_nb, .row = 3, .col = 1},
            ConductorPairEntry{.name = "x_nc", .value = in.x_nc, .row = 3, .col = 2},
            ConductorPairEntry{.name = "x_nn", .value = in.x_nn, .row = 3, .col = 3},
            ConductorPairEntry{.name = "c_aa", .value = in.c_aa, .row = 0, .col = 0},
            ConductorPairEntry{.name = "c_ba", .value = in.c_ba, .row = 1, .col = 0},
            ConductorPairEntry{.name = "c_bb", .value = in.c_bb, .row = 1, .col = 1},
            ConductorPairEntry{.name = "c_ca", .value = in.c_ca, .row = 2, .col = 0},
            ConductorPairEntry{.name = "c_cb", .value = in.c_cb, .row = 2, .col = 1},
            ConductorPairEntry{.name = "c_cc", .value = in.c_cc, .row = 2, .col = 2},
        };
        for (auto const& entry : entries) {
            bool const required = is_present(entry.row) && is_present(entry.col);
            if (required && is_nan(entry.value)) {
                throw InvalidPhaseConnection{name, line_id, phase_letters(phases),
                                             std::format("{} is required", entry.name)};
            }
            if (!required && !is_nan(entry.value)) {
                throw InvalidPhaseConnection{
                    name, line_id, phase_letters(phases),
                    std::format("{} refers to an absent phase or a missing neutral and must be left empty",
                                entry.name)};
            }
        }
    }

    // the absent entries are nan; they are discarded by assignment, here or in masked_inv
    static ComplexTensor<asymmetric_t> compute_reduced_z_series_from_input(AsymLineInput const& asym_line_input,
                                                                           PhaseMask phases) {
        if (!reduced_phase_has_neutral(asym_line_input)) {
            ComplexTensor<asymmetric_t> const r_matrix{asym_line_input.r_aa, asym_line_input.r_bb,
                                                       asym_line_input.r_cc, asym_line_input.r_ba,
                                                       asym_line_input.r_ca, asym_line_input.r_cb};
            ComplexTensor<asymmetric_t> const x_matrix{asym_line_input.x_aa, asym_line_input.x_bb,
                                                       asym_line_input.x_cc, asym_line_input.x_ba,
                                                       asym_line_input.x_ca, asym_line_input.x_cb};
            return r_matrix + 1.0i * x_matrix;
        }
        ComplexTensor4 const r_matrix{asym_line_input.r_aa, asym_line_input.r_bb, asym_line_input.r_cc,
                                      asym_line_input.r_nn, asym_line_input.r_ba, asym_line_input.r_ca,
                                      asym_line_input.r_na, asym_line_input.r_cb, asym_line_input.r_nb,
                                      asym_line_input.r_nc};
        ComplexTensor4 const x_matrix{asym_line_input.x_aa, asym_line_input.x_bb, asym_line_input.x_cc,
                                      asym_line_input.x_nn, asym_line_input.x_ba, asym_line_input.x_ca,
                                      asym_line_input.x_na, asym_line_input.x_cb, asym_line_input.x_nb,
                                      asym_line_input.x_nc};
        ComplexTensor4 const z = r_matrix + 1.0i * x_matrix;
        return kron_reduction(mask_absent_phase_conductors(z, phases));
    }

    // every entry between two present phases is non-zero
    static bool present_entries_nonzero(ComplexTensor<asymmetric_t> const& matrix, PhaseMask phases) {
        for (Idx row = 0; row != 3; ++row) {
            for (Idx col = 0; col != 3; ++col) {
                if (phases.has(row) && phases.has(col) && cabs(matrix(row, col)) < numerical_tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    static ComplexTensor<asymmetric_t>
    compute_z_series_from_input(const power_grid_model::AsymLineInput& asym_line_input) {
        ComplexTensor<asymmetric_t> z_series_abc;
        if (is_nan(asym_line_input.r_na) && is_nan(asym_line_input.x_na)) {
            ComplexTensor<asymmetric_t> const r_matrix{asym_line_input.r_aa, asym_line_input.r_bb,
                                                       asym_line_input.r_cc, asym_line_input.r_ba,
                                                       asym_line_input.r_ca, asym_line_input.r_cb};
            ComplexTensor<asymmetric_t> const x_matrix{asym_line_input.x_aa, asym_line_input.x_bb,
                                                       asym_line_input.x_cc, asym_line_input.x_ba,
                                                       asym_line_input.x_ca, asym_line_input.x_cb};
            z_series_abc = r_matrix + 1.0i * x_matrix;
        } else {
            ComplexTensor4 const r_matrix{asym_line_input.r_aa, asym_line_input.r_bb, asym_line_input.r_cc,
                                          asym_line_input.r_nn, asym_line_input.r_ba, asym_line_input.r_ca,
                                          asym_line_input.r_na, asym_line_input.r_cb, asym_line_input.r_nb,
                                          asym_line_input.r_nc};
            ComplexTensor4 const x_matrix{asym_line_input.x_aa, asym_line_input.x_bb, asym_line_input.x_cc,
                                          asym_line_input.x_nn, asym_line_input.x_ba, asym_line_input.x_ca,
                                          asym_line_input.x_na, asym_line_input.x_cb, asym_line_input.x_nb,
                                          asym_line_input.x_nc};
            ComplexTensor4 const z = r_matrix + 1.0i * x_matrix;
            z_series_abc = kron_reduction(z);
        }
        return z_series_abc;
    }

    static ComplexTensor<asymmetric_t>
    compute_c_matrix_from_input(const power_grid_model::AsymLineInput& asym_line_input) {
        ComplexTensor<asymmetric_t> c_matrix;
        if (!is_nan(asym_line_input.c0) && !is_nan(asym_line_input.c1)) {
            c_matrix = ComplexTensor<asymmetric_t>{(2.0 * asym_line_input.c1 + asym_line_input.c0) / 3.0,
                                                   (asym_line_input.c0 - asym_line_input.c1) / 3.0};
        } else {
            c_matrix = ComplexTensor<asymmetric_t>(asym_line_input.c_aa, asym_line_input.c_bb, asym_line_input.c_cc,
                                                   asym_line_input.c_ba, asym_line_input.c_ca, asym_line_input.c_cb);
        }
        return c_matrix;
    }

    BranchCalcParam<symmetric_t> sym_calc_param() const final {
        if (!phases_.all()) {
            throw ReducedPhaseNotSupported{"Symmetric calculation"};
        }
        DoubleComplex const y1_series =
            average_of_diagonal_of_matrix(y_series_abc_) - average_of_off_diagonal_of_matrix(y_series_abc_);
        DoubleComplex const y1_shunt =
            average_of_diagonal_of_matrix(y_shunt_abc_) - average_of_off_diagonal_of_matrix(y_shunt_abc_);
        return calc_param_y_sym(y1_series, y1_shunt, 1.0);
    }

    BranchCalcParam<asymmetric_t> asym_calc_param() const final {
        BranchCalcParam<asymmetric_t> param{};
        // not both connected
        if (!edge_status()) {
            // single connected
            if (from_status() || to_status()) {
                // branch_shunt = 0.5 * y_shunt + 1.0 / (1.0 / y_series + 2.0 / y_shunt); // NOSONAR(S125)
                auto branch_shunt = ComplexTensor<asymmetric_t>();
                if (phases_.all()) {
                    if ((cabs(y_shunt_abc_) >= numerical_tolerance).all()) {
                        branch_shunt = 0.5 * y_shunt_abc_ + inv(inv(y_series_abc_) + 2.0 * inv(y_shunt_abc_));
                    }
                } else if (present_entries_nonzero(y_shunt_abc_, phases_)) {
                    // same as above on the present phases; the absent rows and columns stay exactly zero
                    branch_shunt = 0.5 * y_shunt_abc_ +
                                   masked_inv(ComplexTensor<asymmetric_t>{masked_inv(y_series_abc_, phases_) +
                                                                          2.0 * masked_inv(y_shunt_abc_, phases_)},
                                              phases_);
                }
                // from or to connected
                param.yff() = from_status() ? branch_shunt : ComplexTensor<asymmetric_t>();
                param.ytt() = to_status() ? branch_shunt : ComplexTensor<asymmetric_t>();
            }
        }
        // both connected
        else {
            param.ytt() = y_series_abc_ + 0.5 * y_shunt_abc_;
            param.yff() = param.ytt();
            param.yft() = -y_series_abc_;
            param.ytf() = -y_series_abc_;
        }
        return param;
    }
};
} // namespace power_grid_model
