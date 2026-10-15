// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#pragma once

#include "../calculation_parameters.hpp"
#include "../common/three_phase_tensor.hpp"
#include "power_grid_model/common/common.hpp"

namespace power_grid_model {

inline ComplexTensor<asymmetric_t> kron_reduction(const ComplexTensor4& matrix_to_reduce) {
    ComplexTensor4 const& Y = matrix_to_reduce;
    auto const Y_aa = ComplexTensor<asymmetric_t>(Y(0, 0), Y(1, 1), Y(2, 2), Y(1, 0), Y(2, 0), Y(2, 1));
    ComplexValue<asymmetric_t> const Y_ab(Y(0, 3), Y(1, 3), Y(2, 3));
    ComplexValue<asymmetric_t> const Y_ba(Y(3, 0), Y(3, 1), Y(3, 2));
    DoubleComplex const Y_bb_inv = 1.0 / Y(3, 3);
    return Y_aa - vector_outer_product(Y_ba, Y_ab) * Y_bb_inv;
}

// inverse of the block of the present phases: the absent rows and columns of matrix are ignored (they may be
// non-finite), those of the result are exactly zero (by assignment, never by multiplication)
inline ComplexTensor<asymmetric_t> masked_inv(ComplexTensor<asymmetric_t> matrix, PhaseMask phases) {
    if (phases.all()) {
        return inv(matrix);
    }
    mask_absent_phases(matrix, phases, phases);
    set_absent_phase_identity(matrix, phases);
    ComplexTensor<asymmetric_t> result = inv(matrix);
    mask_absent_phases(result, phases, phases);
    return result;
}

// replace the rows and columns of the absent phase conductors (0, 1, 2 for a, b, c; the neutral 3 is kept) by the
// identity, so that kron_reduction keeps them decoupled with an identity block
inline ComplexTensor4 mask_absent_phase_conductors(ComplexTensor4 matrix, PhaseMask phases) {
    for (Idx phase = 0; phase != 3; ++phase) {
        if (!phases.has(phase)) {
            matrix.row(phase).setZero();
            matrix.col(phase).setZero();
            matrix(phase, phase) = DoubleComplex{1.0};
        }
    }
    return matrix;
}

} // namespace power_grid_model
