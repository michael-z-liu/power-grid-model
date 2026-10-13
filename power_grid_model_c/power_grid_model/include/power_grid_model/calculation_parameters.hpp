// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "common/common.hpp"
#include "common/enum.hpp"
#include "common/grouped_index_vector.hpp"
#include "common/maybe_owning_view.hpp"
#include "common/small_vector.hpp"
#include "common/statistics.hpp"
#include "common/three_phase_tensor.hpp"

#include <array>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace power_grid_model {
constexpr Idx disconnected = -1;
constexpr IntS status_on = 1;
constexpr IntS status_off = 0;

// Entry of YBus, node addmittance matrix
struct YBusElement {
    YBusElementType element_type{};
    Idx idx{}; // index of the component
};

// everything here always per unit

// branch math calculation parameter and math output
template <symmetry_tag sym_type> struct BranchCalcParam {
    using sym = sym_type;

    std::array<ComplexTensor<sym>, 4> value;
    // getter
    ComplexTensor<sym>& yff() { return value[0]; }
    ComplexTensor<sym> const& yff() const { return value[0]; }
    ComplexTensor<sym>& yft() { return value[1]; }
    ComplexTensor<sym> const& yft() const { return value[1]; }
    ComplexTensor<sym>& ytf() { return value[2]; }
    ComplexTensor<sym> const& ytf() const { return value[2]; }
    ComplexTensor<sym>& ytt() { return value[3]; }
    ComplexTensor<sym> const& ytt() const { return value[3]; }
};

struct BusSolverOutput {
    LimitViolation q_limit_violated{LimitViolation::none};
};

// TODO(mgovers): cleanup v2: branch solver output should always be in current domain; conversion to power domain should
// be done in main_core/output.hpp
template <symmetry_tag sym_type> struct BranchSolverOutput {
    using sym = sym_type;

    ComplexValue<sym> s_f{};
    ComplexValue<sym> s_t{};
    ComplexValue<sym> i_f{};
    ComplexValue<sym> i_t{};
};

template <symmetry_tag sym_type> struct BranchShortCircuitSolverOutput {
    using sym = sym_type;

    ComplexValue<sym> i_f{};
    ComplexValue<sym> i_t{};
};

// fault math calculation parameters and math output
struct FaultCalcParam {
    DoubleComplex y_fault;
    FaultType fault_type{FaultType::nan};
    FaultPhase fault_phase{FaultPhase::nan};
};

template <symmetry_tag sym_type> struct FaultShortCircuitSolverOutput {
    using sym = sym_type;

    ComplexValue<sym> i_fault{};
};

// appliance solver math output, always injection direction
// s > 0, energy appliance -> node
template <symmetry_tag sym_type> struct ApplianceSolverOutput {
    using sym = sym_type;

    ComplexValue<sym> s{}; // TODO(mgovers): cleanup v2: appliance solver output should always be in current domain;
                           // conversion to power domain should be done in main_core/output.hpp
    ComplexValue<sym> i{};
};
template <symmetry_tag sym_type> struct ApplianceShortCircuitSolverOutput {
    using sym = sym_type;

    ComplexValue<sym> i{};
};

struct VoltageRegulatorSolverOutput {
    LimitViolation limit_violated{};

    // provide generator info, as the regulator component has no other access to it
    ID generator_id{};
    IntS generator_status{};
};

// voltage sensor calculation parameters for state estimation
// The value is the complex voltage
// If the imaginary part is NaN, it means the angle calculation is not correct
template <symmetry_tag sym> using VoltageSensorCalcParam = UniformComplexRandVar<sym>;

// power sensor calculation parameters for state estimation
// The value is the complex power
//   * for appliances, it is always in injection direction
//   * for branches, the direction is node -> branch
template <symmetry_tag sym> using PowerSensorCalcParam = DecomposedComplexRandVar<sym>;

// current sensor calculation parameters for state estimation
// The value is the complex current
//   * for branches, the direction is node -> branch
template <symmetry_tag sym_type> struct CurrentSensorCalcParam {
    using sym = sym_type;

    static constexpr bool symmetric{is_symmetric_v<sym>};

    AngleMeasurementType angle_measurement_type{};
    DecomposedComplexRandVar<sym> measurement;
};

template <typename T>
concept sensor_calc_param_type = std::derived_from<T, VoltageSensorCalcParam<symmetric_t>> ||
                                 std::derived_from<T, VoltageSensorCalcParam<asymmetric_t>> ||
                                 std::derived_from<T, PowerSensorCalcParam<symmetric_t>> ||
                                 std::derived_from<T, PowerSensorCalcParam<asymmetric_t>> ||
                                 std::derived_from<T, CurrentSensorCalcParam<symmetric_t>> ||
                                 std::derived_from<T, CurrentSensorCalcParam<asymmetric_t>>;

static_assert(sensor_calc_param_type<VoltageSensorCalcParam<symmetric_t>>);
static_assert(sensor_calc_param_type<VoltageSensorCalcParam<asymmetric_t>>);
static_assert(sensor_calc_param_type<PowerSensorCalcParam<symmetric_t>>);
static_assert(sensor_calc_param_type<PowerSensorCalcParam<asymmetric_t>>);
static_assert(sensor_calc_param_type<CurrentSensorCalcParam<symmetric_t>>);
static_assert(sensor_calc_param_type<CurrentSensorCalcParam<asymmetric_t>>);

struct TransformerTapRegulatorCalcParam {
    double u_set{};
    double u_band{};
    DoubleComplex z_compensation;
    IntS status{};
};

struct TransformerTapPosition {
    ID transformer_id{};
    IntS tap_position{};
};
using TransformerTapPositionOutput = std::vector<TransformerTapPosition>;

// from side, to side
// in case of indices for math model, -1/disconnected means the branch is not connected to that side
using BranchIdx = std::array<Idx, 2>;
// node 0, 1, 2 side
// in case of indices for math model, -1/disconnected means the branch is not connected to that side
using Branch3Idx = std::array<Idx, 3>;

// phase slots of a bus or a branch terminal: a = bit 0, b = bit 1, c = bit 2
// bit operations are confined to this class so that they all operate on unsigned 8-bit values
class PhaseMask {
  public:
    using underlying_type = std::uint8_t;

    constexpr PhaseMask() = default; // no phases
    explicit constexpr PhaseMask(underlying_type value) : bits_{static_cast<underlying_type>(value & all_bits)} {}

    constexpr underlying_type bits() const { return bits_; }
    constexpr bool none() const { return bits_ == 0; }
    constexpr bool all() const { return bits_ == all_bits; }
    // phase is the phase slot 0, 1, 2 for a, b, c; any other value is not present
    constexpr bool has(Idx phase) const {
        if (phase < 0 || phase >= static_cast<Idx>(single_phase_bits.size())) {
            return false;
        }
        return (bits_ & single_phase_bits[static_cast<std::size_t>(phase)]) != 0;
    }
    // all phases of other are present in this mask
    constexpr bool contains(PhaseMask other) const { return (bits_ & other.bits_) == other.bits_; }

    friend constexpr PhaseMask operator|(PhaseMask lhs, PhaseMask rhs) {
        return PhaseMask{static_cast<underlying_type>(lhs.bits_ | rhs.bits_)};
    }
    friend constexpr PhaseMask operator&(PhaseMask lhs, PhaseMask rhs) {
        return PhaseMask{static_cast<underlying_type>(lhs.bits_ & rhs.bits_)};
    }
    friend constexpr bool operator==(PhaseMask const& lhs, PhaseMask const& rhs) = default;

  private:
    static constexpr underlying_type all_bits{0b111};
    static constexpr std::array<underlying_type, 3> single_phase_bits{0b001, 0b010, 0b100};

    underlying_type bits_{};
};
static_assert(sizeof(PhaseMask) == 1);

inline constexpr PhaseMask phase_mask_none{};
inline constexpr PhaseMask phase_mask_a{0b001};
inline constexpr PhaseMask phase_mask_b{0b010};
inline constexpr PhaseMask phase_mask_c{0b100};
inline constexpr PhaseMask phase_mask_ab{0b011};
inline constexpr PhaseMask phase_mask_ac{0b101};
inline constexpr PhaseMask phase_mask_bc{0b110};
inline constexpr PhaseMask phase_mask_abc{0b111};

constexpr std::string_view phase_letters(PhaseMask mask) {
    constexpr std::array<std::string_view, 8> letters{"-", "a", "b", "ab", "c", "ac", "bc", "abc"};
    return letters[mask.bits()];
}

// how the phases of the two terminals of a branch are coupled
//   per_phase: phase k at one terminal is connected to phase k at the other (line-like)
//   winding: all phases of a terminal are coupled through a winding (transformer-like)
enum class PhaseCoupling : IntS { per_phase = 0, winding = 1 };

// non-earth terminals of a branch at each side and how they couple
struct BranchPhaseIncidence {
    PhaseMask from{phase_mask_abc};
    PhaseMask to{phase_mask_abc};
    PhaseCoupling coupling{PhaseCoupling::per_phase}; // per_phase requires from == to

    constexpr bool all_phases() const { return from.all() && to.all(); }
    friend constexpr bool operator==(BranchPhaseIncidence const& lhs, BranchPhaseIncidence const& rhs) = default;
};
inline constexpr BranchPhaseIncidence three_phase_winding{
    .from = phase_mask_abc, .to = phase_mask_abc, .coupling = PhaseCoupling::winding};

struct MathModelTopology {
    Idx slack_bus{};
    bool is_radial{};
    // true if at least one bus has fewer than three phases; only then bus_phase_mask is filled
    bool has_reduced_phase{false};
    std::vector<double> phase_shift;
    // phases present at each bus; size n_bus() if has_reduced_phase, otherwise empty (all buses abc)
    std::vector<PhaseMask> bus_phase_mask;
    std::vector<BranchIdx> branch_bus_idx;
    std::vector<BranchIdx> fill_in;
    DenseGroupedIdxVector sources_per_bus;
    DenseGroupedIdxVector shunts_per_bus;
    SparseGroupedIdxVector load_gens_per_bus;
    std::vector<LoadGenType> load_gen_type;
    DenseGroupedIdxVector voltage_sensors_per_bus;
    DenseGroupedIdxVector power_sensors_per_source;
    DenseGroupedIdxVector power_sensors_per_load_gen;
    DenseGroupedIdxVector power_sensors_per_shunt;
    DenseGroupedIdxVector power_sensors_per_branch_from;
    DenseGroupedIdxVector power_sensors_per_branch_to;
    DenseGroupedIdxVector power_sensors_per_bus;
    DenseGroupedIdxVector current_sensors_per_branch_from;
    DenseGroupedIdxVector current_sensors_per_branch_to;
    DenseGroupedIdxVector tap_regulators_per_branch;
    DenseGroupedIdxVector voltage_regulators_per_load_gen;

    Idx n_bus() const { return static_cast<Idx>(phase_shift.size()); }

    PhaseMask bus_phase(Idx bus) const {
        assert(bus >= 0 && bus < n_bus());
        assert(!has_reduced_phase || std::ssize(bus_phase_mask) == n_bus());
        return has_reduced_phase ? bus_phase_mask[static_cast<std::size_t>(bus)] : phase_mask_abc;
    }

    Idx n_branch() const { return static_cast<Idx>(branch_bus_idx.size()); }

    Idx n_source() const { return sources_per_bus.element_size(); }

    Idx n_shunt() const { return shunts_per_bus.element_size(); }

    Idx n_load_gen() const { return load_gens_per_bus.element_size(); }

    Idx n_voltage_sensor() const { return voltage_sensors_per_bus.element_size(); }

    Idx n_source_power_sensor() const { return power_sensors_per_source.element_size(); }

    Idx n_load_gen_power_sensor() const { return power_sensors_per_load_gen.element_size(); }

    Idx n_shunt_power_power_sensor() const { return power_sensors_per_shunt.element_size(); }

    Idx n_branch_from_power_sensor() const { return power_sensors_per_branch_from.element_size(); }

    Idx n_branch_to_power_sensor() const { return power_sensors_per_branch_to.element_size(); }

    Idx n_bus_power_sensor() const { return power_sensors_per_bus.element_size(); }

    Idx n_branch_from_current_sensor() const { return current_sensors_per_branch_from.element_size(); }

    Idx n_branch_to_current_sensor() const { return current_sensors_per_branch_to.element_size(); }

    Idx n_transformer_tap_regulator() const { return tap_regulators_per_branch.element_size(); }

    Idx n_load_gen_voltage_regulator() const { return voltage_regulators_per_load_gen.element_size(); }
};

struct SourceCalcParam {
    DoubleComplex y1;
    DoubleComplex y0;

    template <symmetry_tag sym> ComplexTensor<sym> y_ref() const {
        if constexpr (is_symmetric_v<sym>) {
            return y1;
        } else {
            return ComplexTensor<asymmetric_t>{(2.0 * y1 + y0) / 3.0, (y0 - y1) / 3.0};
        }
    }
};

template <symmetry_tag sym_type> struct VoltageRegulatorCalcParam {
    using sym = sym_type;

    IntS status{};
    DoubleComplex u_ref;
    double q_min{};
    double q_max{};

    // add generator id for later lookup
    ID generator_id{};
};

template <symmetry_tag sym_type> struct MathModelParam {
    using sym = sym_type;

    std::vector<BranchCalcParam<sym>> branch_param;
    ComplexTensorVector<sym> shunt_param;
    std::vector<SourceCalcParam> source_param;
};

template <symmetry_tag sym_type> struct MathModelParamIncrement {
    using sym = sym_type;

    std::vector<BranchCalcParam<sym>> branch_param;
    ComplexTensorVector<sym> shunt_param;
    std::vector<SourceCalcParam> source_param;

    IdxVector branch_param_to_change; // indices of changed branch_param
    IdxVector shunt_param_to_change;  // indices of changed shunt_param
    IdxVector source_param_to_change; // indices of changed source_param
};

template <typename T>
concept math_model_param_c =
    std::derived_from<T, MathModelParam<symmetric_t>> || std::derived_from<T, MathModelParamIncrement<symmetric_t>> ||
    std::derived_from<T, MathModelParam<asymmetric_t>> || std::derived_from<T, MathModelParamIncrement<asymmetric_t>>;

static_assert(math_model_param_c<MathModelParam<symmetric_t>>);
static_assert(math_model_param_c<MathModelParam<asymmetric_t>>);
static_assert(math_model_param_c<MathModelParamIncrement<symmetric_t>>);
static_assert(math_model_param_c<MathModelParamIncrement<asymmetric_t>>);

template <symmetry_tag sym_type> struct PowerFlowInput {
    using sym = sym_type;

    ComplexVector source;                // Complex u_ref of each source
    ComplexValueVector<sym> s_injection; // Specified injection power of each load_gen
    std::vector<VoltageRegulatorCalcParam<sym>> voltage_regulator;
    IntSVector load_gen_status;
    CalculationInitialization initialization{CalculationInitialization::default_initialization};
};

template <symmetry_tag sym_type> struct StateEstimationInput {
    using sym = sym_type;

    // connection status of shunt, load_gen, source
    // this is needed to determine if a measurement is relevant
    // if the shunt/load_gen/source is disconnected, all its measurements are discarded
    IntSVector shunt_status;
    IntSVector load_gen_status;
    IntSVector source_status;
    // measured value
    std::vector<VoltageSensorCalcParam<sym>> measured_voltage;
    std::vector<PowerSensorCalcParam<sym>> measured_source_power;
    std::vector<PowerSensorCalcParam<sym>> measured_load_gen_power;
    std::vector<PowerSensorCalcParam<sym>> measured_shunt_power;
    std::vector<PowerSensorCalcParam<sym>> measured_branch_from_power;
    std::vector<PowerSensorCalcParam<sym>> measured_branch_to_power;
    std::vector<PowerSensorCalcParam<sym>> measured_bus_injection;
    std::vector<CurrentSensorCalcParam<sym>> measured_branch_from_current;
    std::vector<CurrentSensorCalcParam<sym>> measured_branch_to_current;
};

struct ShortCircuitInput {
    DenseGroupedIdxVector fault_buses;
    std::vector<FaultCalcParam> faults;
    ComplexVector source; // Complex u_ref of each source
};

template <typename T>
concept symmetric_calculation_input_type =
    std::derived_from<T, PowerFlowInput<symmetric_t>> || std::derived_from<T, StateEstimationInput<symmetric_t>>;

static_assert(symmetric_calculation_input_type<PowerFlowInput<symmetric_t>>);
static_assert(symmetric_calculation_input_type<StateEstimationInput<symmetric_t>>);
static_assert(!symmetric_calculation_input_type<PowerFlowInput<asymmetric_t>>);
static_assert(!symmetric_calculation_input_type<StateEstimationInput<asymmetric_t>>);
static_assert(!symmetric_calculation_input_type<ShortCircuitInput>);

template <typename T>
concept asymmetric_calculation_input_type =
    std::derived_from<T, PowerFlowInput<asymmetric_t>> || std::derived_from<T, StateEstimationInput<asymmetric_t>> ||
    std::derived_from<T, ShortCircuitInput>;

static_assert(!asymmetric_calculation_input_type<PowerFlowInput<symmetric_t>>);
static_assert(!asymmetric_calculation_input_type<StateEstimationInput<symmetric_t>>);
static_assert(asymmetric_calculation_input_type<PowerFlowInput<asymmetric_t>>);
static_assert(asymmetric_calculation_input_type<StateEstimationInput<asymmetric_t>>);
static_assert(asymmetric_calculation_input_type<ShortCircuitInput>);

template <typename T>
concept calculation_input_type = symmetric_calculation_input_type<T> || asymmetric_calculation_input_type<T>;

static_assert(calculation_input_type<PowerFlowInput<symmetric_t>>);
static_assert(calculation_input_type<StateEstimationInput<symmetric_t>>);
static_assert(calculation_input_type<PowerFlowInput<asymmetric_t>>);
static_assert(calculation_input_type<StateEstimationInput<asymmetric_t>>);
static_assert(calculation_input_type<ShortCircuitInput>);

struct solver_output_t {};

template <symmetry_tag sym_type> struct SolverOutput {
    using type = solver_output_t;
    using sym = sym_type;

    std::vector<ComplexValue<sym>> u;
    std::vector<ComplexValue<sym>> bus_injection; // TODO(mgovers): remove this for v2
    std::vector<BusSolverOutput> bus;
    std::vector<BranchSolverOutput<sym>> branch;
    std::vector<ApplianceSolverOutput<sym>> source;
    std::vector<ApplianceSolverOutput<sym>> shunt;
    std::vector<ApplianceSolverOutput<sym>> load_gen;
    std::vector<VoltageRegulatorSolverOutput> voltage_regulator;
};

template <symmetry_tag sym_type> struct ShortCircuitSolverOutput {
    using type = solver_output_t;
    using sym = sym_type;

    std::vector<ComplexValue<sym>> u_bus;
    std::vector<FaultShortCircuitSolverOutput<sym>> fault;
    std::vector<BranchShortCircuitSolverOutput<sym>> branch;
    std::vector<ApplianceShortCircuitSolverOutput<sym>> source;
    std::vector<ApplianceShortCircuitSolverOutput<sym>> shunt;
};

template <typename T>
concept solver_output_type = std::derived_from<typename T::type, solver_output_t>;

static_assert(solver_output_type<SolverOutput<symmetric_t>>);
static_assert(solver_output_type<SolverOutput<asymmetric_t>>);
static_assert(solver_output_type<ShortCircuitSolverOutput<symmetric_t>>);
static_assert(solver_output_type<ShortCircuitSolverOutput<asymmetric_t>>);

template <typename T>
concept symmetric_solver_output_type = solver_output_type<T> && is_symmetric_v<typename T::sym>;

static_assert(symmetric_solver_output_type<SolverOutput<symmetric_t>>);
static_assert(!symmetric_solver_output_type<SolverOutput<asymmetric_t>>);
static_assert(symmetric_solver_output_type<ShortCircuitSolverOutput<symmetric_t>>);
static_assert(!symmetric_solver_output_type<ShortCircuitSolverOutput<asymmetric_t>>);

template <typename T>
concept asymmetric_solver_output_type = solver_output_type<T> && is_asymmetric_v<typename T::sym>;

static_assert(!asymmetric_solver_output_type<SolverOutput<symmetric_t>>);
static_assert(asymmetric_solver_output_type<SolverOutput<asymmetric_t>>);
static_assert(!asymmetric_solver_output_type<ShortCircuitSolverOutput<symmetric_t>>);
static_assert(asymmetric_solver_output_type<ShortCircuitSolverOutput<asymmetric_t>>);

template <typename T>
concept steady_state_solver_output_type = solver_output_type<T> && std::derived_from<T, SolverOutput<typename T::sym>>;

static_assert(steady_state_solver_output_type<SolverOutput<symmetric_t>>);
static_assert(steady_state_solver_output_type<SolverOutput<asymmetric_t>>);
static_assert(!steady_state_solver_output_type<ShortCircuitSolverOutput<symmetric_t>>);
static_assert(!steady_state_solver_output_type<ShortCircuitSolverOutput<asymmetric_t>>);

template <typename T>
concept short_circuit_solver_output_type =
    solver_output_type<T> && std::derived_from<T, ShortCircuitSolverOutput<typename T::sym>>;

static_assert(!short_circuit_solver_output_type<SolverOutput<symmetric_t>>);
static_assert(!short_circuit_solver_output_type<SolverOutput<asymmetric_t>>);
static_assert(short_circuit_solver_output_type<ShortCircuitSolverOutput<symmetric_t>>);
static_assert(short_circuit_solver_output_type<ShortCircuitSolverOutput<asymmetric_t>>);

struct OptimizerOutput {
    TransformerTapPositionOutput transformer_tap_positions;
};

template <typename T> struct SupernodeOutput;

// One entry per user node of a topological node. Only link-merged supernodes hold more than one, so
// inline room for a single element keeps the ordinary node from allocating at all.
template <symmetry_tag sym> using UserNodeValueVector = SmallVector<ComplexValue<sym>, 1>;

template <steady_state_solver_output_type SolverOutputType> struct SupernodeOutput<SolverOutputType> {
    using sym = decode_symmetry_v<SolverOutputType>;

    UserNodeValueVector<sym> bus_injection;    // user bus output
    std::vector<BranchSolverOutput<sym>> link; // user link
};
template <short_circuit_solver_output_type SolverOutputType> struct SupernodeOutput<SolverOutputType> {
    using sym = decode_symmetry_v<SolverOutputType>;

    std::vector<BranchShortCircuitSolverOutput<sym>> link; // user link
};

template <typename T> struct MathOutput {
    using SolverOutputType = T;
    using UnderlyingSolverOutputType = underlying_value_t<SolverOutputType>;
    using sym = decode_symmetry_v<SolverOutputType>;

    SolverOutputType solver_output;
    OptimizerOutput optimizer_output;
    std::vector<SupernodeOutput<UnderlyingSolverOutputType>> supernode_output;
};

// component indices at physical model side
// from, to node indices for branches
// node1, node2, node3 indices for 3-way branches
// node indices for source, load_gen, shunt
struct ComponentTopology {
    Idx n_node{};
    std::vector<BranchIdx> branch_node_idx;
    std::vector<Branch3Idx> branch3_node_idx;
    std::vector<BranchIdx> link_node_idx;
    IdxVector shunt_node_idx;
    IdxVector source_node_idx;
    IdxVector load_gen_node_idx;
    std::vector<LoadGenType> load_gen_type;
    IdxVector voltage_sensor_node_idx;
    IdxVector power_sensor_object_idx; // the index is relative to branch, source, shunt or load_gen
    std::vector<MeasuredTerminalType> power_sensor_terminal_type;
    IdxVector current_sensor_object_idx; // the index is relative to branch
    std::vector<MeasuredTerminalType> current_sensor_terminal_type;
    std::vector<ComponentType> regulator_type;
    IdxVector regulated_object_idx; // the index is relative to branch or branch3
    std::vector<ComponentType> regulated_object_type;

    constexpr Idx n_node_total() const { return n_node + std::ssize(branch3_node_idx); }
};

struct ReducedComponentTopology {
    Idx n_node{};
    MaybeOwningConstVector<BranchIdx> branch_node_idx;
    MaybeOwningConstVector<Branch3Idx> branch3_node_idx;
    MaybeOwningConstVector<Idx> shunt_node_idx;
    MaybeOwningConstVector<Idx> source_node_idx;
    MaybeOwningConstVector<Idx> load_gen_node_idx;
    MaybeOwningConstVector<LoadGenType> load_gen_type;
    MaybeOwningConstVector<Idx> voltage_sensor_node_idx;
    MaybeOwningConstVector<Idx> power_sensor_object_idx; // the index is relative to branch, source, shunt or load_gen
    MaybeOwningConstVector<MeasuredTerminalType> power_sensor_terminal_type;
    MaybeOwningConstVector<Idx> current_sensor_object_idx; // the index is relative to branch
    MaybeOwningConstVector<MeasuredTerminalType> current_sensor_terminal_type;
    MaybeOwningConstVector<ComponentType> regulator_type;
    MaybeOwningConstVector<Idx> regulated_object_idx; // the index is relative to branch or branch3
    MaybeOwningConstVector<ComponentType> regulated_object_type;

    constexpr Idx n_node_total() const { return n_node + std::ssize(branch3_node_idx); }

    static ReducedComponentTopology from_component_topology(ComponentTopology const& comp_topo) {
        assert(comp_topo.link_node_idx.empty() && "link is not supported in reduced topology");
        return ReducedComponentTopology{
            .n_node = comp_topo.n_node,
            .branch_node_idx = std::span{comp_topo.branch_node_idx},
            .branch3_node_idx = std::span{comp_topo.branch3_node_idx},
            .shunt_node_idx = std::span{comp_topo.shunt_node_idx},
            .source_node_idx = std::span{comp_topo.source_node_idx},
            .load_gen_node_idx = std::span{comp_topo.load_gen_node_idx},
            .load_gen_type = std::span{comp_topo.load_gen_type},
            .voltage_sensor_node_idx = std::span{comp_topo.voltage_sensor_node_idx},
            .power_sensor_object_idx = std::span{comp_topo.power_sensor_object_idx},
            .power_sensor_terminal_type = std::span{comp_topo.power_sensor_terminal_type},
            .current_sensor_object_idx = std::span{comp_topo.current_sensor_object_idx},
            .current_sensor_terminal_type = std::span{comp_topo.current_sensor_terminal_type},
            .regulator_type = std::span{comp_topo.regulator_type},
            .regulated_object_idx = std::span{comp_topo.regulated_object_idx},
            .regulated_object_type = std::span{comp_topo.regulated_object_type},
        };
    }
};

struct TopologicalNode {
    IdxVector user_nodes;
    std::vector<BranchIdx> internal_links;

    constexpr auto is_supernode() const noexcept -> bool { return user_nodes.size() > 1 && !internal_links.empty(); }
};

struct ComponentToTopoNodeCoupling {
    Idx n_topo_nodes{};
    // for every user node: which topo node it belongs to and which index it has within that topo node
    std::vector<Idx2D> user_nodes_to_topo_nodes;
    // for every user link: which topo node it belongs to and which index it has within that topo node
    std::vector<Idx2D> user_links_to_topo_nodes;
};

struct TopologicalNodesAndCoupling {
    std::vector<TopologicalNode> topo_nodes;
    ComponentToTopoNodeCoupling coupling;
};

struct ReducedTopology {
    ReducedComponentTopology reduced_comp_topo;
    TopologicalNodesAndCoupling topo_node_coup;
};

// connection property
using BranchConnected = std::array<IntS, 2>;
using Branch3Connected = std::array<IntS, 3>;

// component connection property at model side
struct ComponentConnections {
    std::vector<BranchConnected> branch_connected;
    std::vector<Branch3Connected> branch3_connected;
    std::vector<BranchConnected> link_connected;
    DoubleVector branch_phase_shift;
    // one entry per branch_connected entry; empty means all branches are abc -> abc, per phase
    std::vector<BranchPhaseIncidence> branch_phase_incidence;
    // 3-way branch, phase shift = phase_node_x - phase_internal_node
    std::vector<std::array<double, 3>> branch3_phase_shift;
    IntSVector source_connected;
};

// To couple 3-way branch, to math model, 3 virtual branches are created
// they are always in the same math model group
// sequence numbers of 3 virtual branches are saved in pos
struct Idx2DBranch3 {
    Idx group{};
    // 0: node 0 -> internal node
    // 1: node 1 -> internal node
    // 2: node 2 -> internal node
    std::array<Idx, 3> pos{};

    friend constexpr bool operator==(Idx2DBranch3 const& x, Idx2DBranch3 const& y) = default;
};

// couple component to math model
// use Idx2D to map component to math model
//		group = math model sequence number,
//		group = -1/disconnected means isolated component
//		pos = sequence number in math model,
//		pos = -1/disconnected means not connected at that side, only applicable for branches
struct ComponentToMathCoupling {
    std::vector<Idx2D> fault;
};

// couple component to math model
// like ComponentToMathCoupling but for components that are immutable after the topology is fixed
// use Idx2D to map component to math model
//		group = math model sequence number,
//		group = -1/disconnected means isolated component
//		pos = sequence number in math model,
//		pos = -1/disconnected means not connected at that side, only applicable for branches
struct TopologicalComponentToMathCoupling {
    std::vector<Idx2D> node;
    std::vector<Idx2D> branch;
    std::vector<Idx2DBranch3> branch3;
    std::vector<Idx2D> shunt;
    std::vector<Idx2D> load_gen;
    std::vector<Idx2D> source;
    std::vector<Idx2D> voltage_sensor;
    std::vector<Idx2D> power_sensor;   // can be coupled to branch-from/to, source, load_gen, or shunt sensor
    std::vector<Idx2D> current_sensor; // can be coupled to branch-from/to
    std::vector<Idx2D> regulator;
    std::vector<Idx2D> voltage_regulator;
};

} // namespace power_grid_model
