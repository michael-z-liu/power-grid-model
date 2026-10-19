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
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <format>
#include <string_view>
#include <utility>

namespace power_grid_model {

// single-phase two-winding transformer: each winding is connected between two phase slots of its node (ab, ac, bc)
// or between one phase slot and earth (a, b, c), the latter directly or through an earthing impedance
// the two coupled windings are stamped onto the slots of both nodes (C^T Y_prim C, in physical units) and scaled to
// per unit with the phase-to-earth base voltage u_rated / sqrt3 of each node
// the to winding can instead be centre-tapped (ab_center_tapped, split phase): two halves, a to the centre tap and the
// centre tap to b, with the centre tap earthed directly or through an earthing impedance; the primary and the two
// halves are then three coupled windings
class SinglePhaseTransformer : public Branch {
  private:
    using Matrix2 = Eigen::Matrix<DoubleComplex, 2, 2>;
    using Matrix3 = Eigen::Matrix<DoubleComplex, 3, 3>;
    // terminals 0, 1, 2: from slots a, b, c; 3, 4, 5: to slots a, b, c; 6, 7: earthed terminal of the from and to
    // winding behind an earthing impedance
    using Matrix8 = Eigen::Matrix<DoubleComplex, 8, 8>;
    static constexpr Idx earth{-1};
    static constexpr Idx to_offset{3};
    static constexpr Idx electrode_from{6};
    static constexpr Idx electrode_to{7};

    enum class WindingKind : std::uint8_t { phase_to_earth, phase_to_phase, center_tapped };

    // winding of one side: the phase slots of the terminals and the angle of the winding voltage of a balanced
    // positive-sequence set (U_a = 1, U_b = a^2, U_c = a); for a centre-tapped winding, the angle of the voltage of
    // half 1 (slot a to the centre tap)
    struct Winding {
        PhaseMask slots;
        Idx positive{};      // slot 0, 1, 2
        Idx negative{earth}; // slot 0, 1, 2, or earth for a phase-to-earth winding; slot b for a centre-tapped winding
        double angle{};
        WindingKind kind{WindingKind::phase_to_earth};

        constexpr bool phase_to_phase() const { return kind == WindingKind::phase_to_phase; }
        constexpr bool center_tapped() const { return kind == WindingKind::center_tapped; }
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
          from_{checked_from_winding(id(), transformer_input.from_connection)},
          to_{winding_of(transformer_input.to_connection)},
          uk_half_{checked_half(id(), to_, transformer_input.uk_half, "uk_half")},
          pk_half_{checked_half(id(), to_, transformer_input.pk_half, "pk_half")},
          clock_{checked_clock(id(), transformer_input.clock)},
          tap_side_{is_nan(transformer_input.tap_side) ? BranchSide::from : transformer_input.tap_side},
          tap_min_{transformer_input.tap_min},
          tap_max_{transformer_input.tap_max},
          tap_nom_{zero_if_nan(transformer_input.tap_nom)},
          tap_direction_{tap_direction_of(tap_min_, tap_max_, tap_nom_)},
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
    // a centre-tapped to winding gives its node the split-phase slot pattern (legs a and b 180 degrees apart)
    BranchPhaseIncidence phase_incidence() const final {
        return {.from = from_.slots,
                .to = to_.slots,
                .coupling = PhaseCoupling::winding,
                .to_pattern = to_.center_tapped() ? SlotPattern::split_phase : SlotPattern::three_phase};
    }

    // hides Edge::get_output: the loading is the larger apparent power of the two sides, each from the complex sum of
    // its per-phase powers; the sum of the per-phase apparent powers would overstate the loading of a phase-to-phase
    // winding by 2 / sqrt3; each half of a centre-tapped winding is rated sn / 2, so that side counts twice the
    // apparent power of its larger leg
    template <symmetry_tag sym> BranchOutput<sym> get_output(BranchSolverOutput<sym> const& solver_output) const {
        BranchOutput<sym> output = Edge::template get_output<sym>(solver_output);
        output.loading =
            std::max(side_power(solver_output.s_f), to_side_power(solver_output.s_t, to_.center_tapped())) *
            base_power<sym> / sn_;
        return output;
    }

    // getters
    constexpr IntS tap_pos() const { return tap_pos_; }
    constexpr BranchSide tap_side() const { return tap_side_; }
    // a tap limit that is not provided follows the tap position on its own side of tap_nom, so it never limits it:
    // below tap_nom for tap_min of a normal tap changer, above tap_nom for tap_min of a reversed one, and conversely
    // for tap_max
    constexpr IntS tap_min() const {
        if (tap_min_ != na_IntS) {
            return tap_min_;
        }
        return tap_direction_ > 0 ? std::min(tap_pos_, tap_nom_) : std::max(tap_pos_, tap_nom_);
    }
    constexpr IntS tap_max() const {
        if (tap_max_ != na_IntS) {
            return tap_max_;
        }
        return tap_direction_ > 0 ? std::max(tap_pos_, tap_nom_) : std::min(tap_pos_, tap_nom_);
    }
    constexpr IntS tap_nom() const { return tap_nom_; }
    constexpr IntS clock() const { return clock_; }
    // both tap limits are provided; an automatic tap regulator needs them, because a limit that is not provided
    // follows the tap position (see tap_min() and tap_max())
    constexpr bool tap_limits_provided() const { return tap_min_ != na_IntS && tap_max_ != na_IntS; }

    // the voltage that a tap regulator controls at the winding of the given side, from the voltages v of the slots of
    // its node (the compensated voltages U + Z_c I for line drop compensation), in per unit of the phase-to-earth
    // base of that node: |V_p| for a phase-to-earth winding, |V_p - V_q| for a phase-to-phase winding and the mean
    // (|V_a| + |V_b|) / 2 of the two legs for a centre-tapped winding; only magnitudes are used, so the clock does not
    // matter
    double controlled_voltage(BranchSide side, ComplexValue<asymmetric_t> const& v) const {
        Winding const& winding = [this, side]() -> Winding const& {
            switch (side) {
            case BranchSide::from:
                return from_;
            case BranchSide::to:
                return to_;
            default:
                throw MissingCaseForEnumError{"single_phase_transformer controlled_voltage", side};
            }
        }();
        switch (winding.kind) {
        case WindingKind::phase_to_earth:
            return cabs(v(winding.positive));
        case WindingKind::phase_to_phase:
            return cabs(v(winding.positive) - v(winding.negative));
        case WindingKind::center_tapped:
            return 0.5 * (cabs(v(winding.positive)) + cabs(v(winding.negative)));
        default:
            throw MissingCaseForEnumError{"single_phase_transformer controlled_voltage", winding.kind};
        }
    }
    // a single-phase transformer only takes part in asymmetric calculations; symmetric power flow is rejected before
    // any tap regulation
    static double controlled_voltage(BranchSide /* side */, ComplexValue<symmetric_t> const& /* v */) {
        throw ReducedPhaseNotSupported{"Symmetric power flow"};
    }

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
    double uk_half_; // nan: not provided, the default split of uk and pk
    double pk_half_; // nan: not provided, the default split of uk and pk
    IntS clock_;
    BranchSide tap_side_;
    IntS tap_pos_{};
    IntS tap_min_; // na_IntS: not provided, see tap_min()
    IntS tap_max_; // na_IntS: not provided, see tap_max()
    IntS tap_nom_;
    IntS tap_direction_; // 1: a higher tap position raises the voltage of the tapped winding, -1: reversed
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
    static double to_side_power(ComplexValue<symmetric_t> const& s, bool /* center_tapped */) { return cabs(s); }
    static double to_side_power(ComplexValue<asymmetric_t> const& s, bool center_tapped) {
        return center_tapped ? 2.0 * std::max(cabs(s(0)), cabs(s(1))) : side_power(s);
    }

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
            return {.slots = phase_mask_ab,
                    .positive = 0,
                    .negative = 1,
                    .angle = deg_30,
                    .kind = WindingKind::phase_to_phase};
        case ac:
            return {.slots = phase_mask_ac,
                    .positive = 0,
                    .negative = 2,
                    .angle = -deg_30,
                    .kind = WindingKind::phase_to_phase};
        case bc:
            return {.slots = phase_mask_bc,
                    .positive = 1,
                    .negative = 2,
                    .angle = -3.0 * deg_30,
                    .kind = WindingKind::phase_to_phase};
        case ab_center_tapped:
            return {
                .slots = phase_mask_ab, .positive = 0, .negative = 1, .angle = 0.0, .kind = WindingKind::center_tapped};
        default:
            throw MissingCaseForEnumError{"single_phase_transformer connection", connection};
        }
    }

    // a centre-tapped winding is only allowed on the to side
    static Winding checked_from_winding(ID transformer_id, WindingConnection1P connection) {
        Winding const winding = winding_of(connection);
        if (winding.center_tapped()) {
            throw InvalidPhaseConnection{name, transformer_id, phase_letters(winding.slots),
                                         "a centre-tapped winding is only allowed on the to side"};
        }
        return winding;
    }

    // uk_half and pk_half only apply to a centre-tapped to winding; on any other they must be left empty
    static double checked_half(ID transformer_id, Winding const& to, double value, std::string_view attribute) {
        if (!to.center_tapped() && !is_nan(value)) {
            throw InvalidPhaseConnection{
                name, transformer_id, phase_letters(to.slots),
                std::format("{} refers to a centre-tapped winding and must be left empty", attribute)};
        }
        return value;
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

    // earthing impedance in ohm, not provided means solidly earthed; only a phase-to-earth winding (its earthed end)
    // and a centre-tapped winding (its centre tap) have an earthed terminal, so the earthing impedance of a
    // phase-to-phase winding must be left empty
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

    // the tap direction, decided once from the provided limits: reversed (tap_min at the high end of the numbering)
    // when both limits are provided and tap_max < tap_min, when only tap_min is provided and tap_min > tap_nom, or
    // when only tap_max is provided and tap_max < tap_nom; otherwise, including without limits, normal
    static constexpr IntS tap_direction_of(IntS tap_min, IntS tap_max, IntS tap_nom) {
        bool const has_min = tap_min != na_IntS;
        bool const has_max = tap_max != na_IntS;
        bool reversed = false;
        if (has_min && has_max) {
            reversed = tap_max < tap_min;
        } else if (has_min) {
            reversed = tap_min > tap_nom;
        } else if (has_max) {
            reversed = tap_max < tap_nom;
        }
        return reversed ? IntS{-1} : IntS{1};
    }

    // limit the tap position to the provided limits, in the tap direction; a limit that is not provided leaves its
    // side unbounded
    constexpr IntS tap_limit(IntS new_tap) const {
        // the provided limits at the low and high end of the numbering
        IntS const low = tap_direction_ > 0 ? tap_min_ : tap_max_;
        IntS const high = tap_direction_ > 0 ? tap_max_ : tap_min_;
        if (high != na_IntS) {
            new_tap = std::min(new_tap, high);
        }
        if (low != na_IntS) {
            new_tap = std::max(new_tap, low);
        }
        return new_tap;
    }

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

    // the winding voltages u1' and u2' at the tap position
    std::pair<double, double> tap_adjusted_voltages() const {
        double u1 = u1_;
        double u2 = u2_;
        double const tap_offset = tap_direction_ * (tap_pos_ - tap_nom_) * tap_size_;
        if (tap_side_ == BranchSide::from) {
            u1 += tap_offset;
        } else {
            u2 += tap_offset;
        }
        return {u1, u2};
    }

    // the short circuit parameters uk' and pk' at the tap position
    std::pair<double, double> tap_adjusted_short_circuit() const {
        return {tap_adjust_impedance(tap_pos_, tap_min(), tap_max(), tap_nom_, uk_, uk_min_, uk_max_),
                tap_adjust_impedance(tap_pos_, tap_min(), tap_max(), tap_nom_, pk_, pk_min_, pk_max_)};
    }

    // primitive admittance of the two windings in siemens (winding currents from winding voltages), including the
    // off-nominal tap; with one side open, the open winding carries no current and the connected winding sees the
    // open-end shunt 1/2 y_m + (1/y_s + 2/y_m)^-1
    Matrix2 primitive_admittance() const {
        auto const [u1, u2] = tap_adjusted_voltages();
        auto const [uk, pk] = tap_adjusted_short_circuit();

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

    // series impedance in per unit on sn from the short circuit voltage and the load losses, as for transformer
    static DoubleComplex short_circuit_impedance(double uk, double pk, double sn) {
        double const r = pk / sn;
        double const x_squared = uk * uk - r * r;
        return {r, (uk >= 0 ? 1.0 : -1.0) * (x_squared > 0.0 ? std::sqrt(x_squared) : 0.0)};
    }

    // short circuit impedance between the primary and one half of a centre-tapped winding (the other half open), in
    // per unit on sn and the half voltage, at the nominal tap; uk_half and pk_half are on the pair rating sn / 2, so
    // r = 4 pk_half / sn and |z| = 2 uk_half; not provided, each follows the split of Kersting (primary 0.5 r + j0.8 x,
    // each half r + j0.4 x), which gives 1.5 r + j1.2 x
    DoubleComplex half_winding_impedance() const {
        DoubleComplex const z_nominal = short_circuit_impedance(uk_, pk_, sn_);
        double const uk_half =
            is_nan(uk_half_)
                ? (uk_ >= 0 ? 0.5 : -0.5) * cabs(DoubleComplex{1.5 * z_nominal.real(), 1.2 * z_nominal.imag()})
                : uk_half_;
        double const pk_half = is_nan(pk_half_) ? 0.375 * pk_ : pk_half_;
        return short_circuit_impedance(2.0 * uk_half, 4.0 * pk_half, sn_);
    }

    // primitive admittance in siemens of the primary (0) and the two halves (1, 2) of a centre-tapped winding, in the
    // N-winding form A^T Z_B^-1 A of OpenDSS, with Z_B = z_0 1 + z_b I in per unit on sn; z_b is the branch of each
    // half, fixed at the nominal tap, and z_0 = z_f - z_b / 2 the branch of the primary, which takes the tap
    // dependence of the full-winding impedance z_f; the magnetising admittance is diag(y_m / 2, y_m / 4, y_m / 4)
    // with the from side open, the primary is eliminated and the halves stay coupled (a balancer); with the to side
    // open, both halves are eliminated
    Matrix3 primitive_admittance_center_tapped() const {
        auto const [u1, u2] = tap_adjusted_voltages();
        auto const [uk, pk] = tap_adjusted_short_circuit();
        DoubleComplex const z_full = short_circuit_impedance(uk, pk, sn_);
        DoubleComplex const z_half_branch = 2.0 * (half_winding_impedance() - short_circuit_impedance(uk_, pk_, sn_));
        DoubleComplex const z_primary_branch = z_full - 0.5 * z_half_branch;
        double const g_shunt = p0_ / sn_;
        double const b_shunt_squared = i0_ * i0_ - g_shunt * g_shunt;
        DoubleComplex const y_shunt{g_shunt, b_shunt_squared > 0.0 ? -std::sqrt(b_shunt_squared) : 0.0};

        // per unit on sn and the winding voltages u1' and u2' / 2, scaled to siemens with sn / (V_i V_j)
        double const u_half = 0.5 * u2;
        DoubleComplex const y_primary = (1.0 / z_full + 0.5 * y_shunt) * (sn_ / (u1 * u1));
        DoubleComplex const y_mutual = -0.5 / z_full * (sn_ / (u1 * u_half));
        DoubleComplex const denominator = 2.0 * z_half_branch * z_full;
        DoubleComplex const y_half =
            ((z_primary_branch + z_half_branch) / denominator + 0.25 * y_shunt) * (sn_ / (u_half * u_half));
        DoubleComplex const y_halves = -z_primary_branch / denominator * (sn_ / (u_half * u_half));
        Matrix3 y_prim;
        y_prim << y_primary, y_mutual, y_mutual, y_mutual, y_half, y_halves, y_mutual, y_halves, y_half;

        if (edge_status()) {
            return y_prim;
        }
        if (from_status()) {
            // the tolerance applies to the magnetising admittance in per unit of the to side, as in the two-winding
            // case
            if (cabs(y_shunt) * sn_ / u2 / u2 * base_u_to_ * base_u_to_ / base_power_1p < numerical_tolerance) {
                return Matrix3::Zero();
            }
            kron_eliminate(y_prim, 1);
            kron_eliminate(y_prim, 2);
            return y_prim;
        }
        if (to_status()) {
            kron_eliminate(y_prim, 0);
            return y_prim;
        }
        return Matrix3::Zero();
    }

    // a non-zero earthing impedance is only possible on a phase-to-earth or centre-tapped winding, whose earthed
    // terminal is then the electrode
    void add_electrodes(Matrix8& y) const {
        if (z_grounding_from_ != DoubleComplex{}) {
            y(electrode_from, electrode_from) += 1.0 / z_grounding_from_;
            kron_eliminate(y, electrode_from);
        }
        if (z_grounding_to_ != DoubleComplex{}) {
            y(electrode_to, electrode_to) += 1.0 / z_grounding_to_;
            kron_eliminate(y, electrode_to);
        }
    }

    // nodal admittance of a centre-tapped to winding: half 1 from slot a to the centre tap, half 2 from the centre tap
    // to slot b; clock 6 reverses both halves
    Matrix8 nodal_admittance_center_tapped() const {
        Matrix3 const y_prim = primitive_admittance_center_tapped();
        Idx const centre_tap = z_grounding_to_ == DoubleComplex{} ? earth : electrode_to;
        std::array<Terminals, 3> terminals{terminals_of(from_, 0, electrode_from, z_grounding_from_),
                                           Terminals{.positive = to_offset + to_.positive, .negative = centre_tap},
                                           Terminals{.positive = centre_tap, .negative = to_offset + to_.negative}};
        if (clock_ == IntS{6}) {
            std::swap(terminals[1].positive, terminals[1].negative);
            std::swap(terminals[2].positive, terminals[2].negative);
        }

        Matrix8 y = Matrix8::Zero();
        for (size_t row = 0; row != terminals.size(); ++row) {
            for (size_t col = 0; col != terminals.size(); ++col) {
                stamp(y, terminals[row], terminals[col], y_prim(static_cast<Idx>(row), static_cast<Idx>(col)));
            }
        }
        add_electrodes(y);
        return y;
    }

    // nodal admittance in siemens of the from and to slots: the windings stamped onto their terminals, the earthing
    // impedances added and their terminals eliminated
    Matrix8 nodal_admittance() const {
        if (to_.center_tapped()) {
            return nodal_admittance_center_tapped();
        }
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
        add_electrodes(y);
        return y;
    }

    // a 3x3 block of the nodal admittance in per unit: Y_pu(s, t) = Y(s, t) U_b,s U_b,t / S_b,1p
    static ComplexTensor<asymmetric_t> per_unit_block(Matrix8 const& y, Idx row_offset, Idx col_offset,
                                                      double base_u_row, double base_u_col) {
        double const scale = base_u_row * base_u_col / base_power_1p;
        return ComplexTensor<asymmetric_t>{(y.block<3, 3>(row_offset, col_offset) * scale).array()};
    }

    // the common mode of a phase-to-phase winding has no earth reference: add a small admittance to earth of
    // y_low / 2 on each entry of its two slots, the analogue of the zero-sequence low susceptance of transformer; a
    // centre-tapped winding is earthed at its centre tap and gets none
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
