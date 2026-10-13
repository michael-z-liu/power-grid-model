// SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
//
// SPDX-License-Identifier: MPL-2.0

#include <power_grid_model/common/matrix_utils.hpp>

#include <power_grid_model/common/common.hpp>
#include <power_grid_model/common/three_phase_tensor.hpp>

#include <doctest/doctest.h>

#include <complex>

namespace power_grid_model {

using namespace std::complex_literals;

TEST_CASE("Test matrix utils") {
    SUBCASE("Distinct complex entries") {
        // M(i, j) = (3i + j + 1) + (9 - 3i - j) * 1j
        ComplexTensor<asymmetric_t> matrix;
        matrix << DoubleComplex{1.0, 9.0}, DoubleComplex{2.0, 8.0}, DoubleComplex{3.0, 7.0}, // row 0
            DoubleComplex{4.0, 6.0}, DoubleComplex{5.0, 5.0}, DoubleComplex{6.0, 4.0},       // row 1
            DoubleComplex{7.0, 3.0}, DoubleComplex{8.0, 2.0}, DoubleComplex{9.0, 1.0};       // row 2

        // diagonal: (1+9j, 5+5j, 9+1j); off-diagonal: the remaining six elements
        DoubleComplex const expected_diagonal_average{5.0, 5.0};
        DoubleComplex const expected_off_diagonal_average{5.0, 5.0};

        CHECK(cabs(average_of_diagonal_of_matrix(matrix) - expected_diagonal_average) < numerical_tolerance);
        CHECK(cabs(average_of_off_diagonal_of_matrix(matrix) - expected_off_diagonal_average) < numerical_tolerance);
    }

    SUBCASE("Each off-diagonal element contributes with weight 1/6") {
        for (Idx row = 0; row < 3; ++row) {
            for (Idx col = 0; col < 3; ++col) {
                if (row == col) {
                    continue;
                }
                CAPTURE(row);
                CAPTURE(col);
                ComplexTensor<asymmetric_t> matrix = ComplexTensor<asymmetric_t>::Zero();
                matrix(row, col) = 6.0;

                CHECK(cabs(average_of_off_diagonal_of_matrix(matrix) - 1.0) < numerical_tolerance);
                CHECK(cabs(average_of_diagonal_of_matrix(matrix)) < numerical_tolerance);
            }
        }
    }

    SUBCASE("Each diagonal element contributes with weight 1/3") {
        for (Idx i = 0; i < 3; ++i) {
            CAPTURE(i);
            ComplexTensor<asymmetric_t> matrix = ComplexTensor<asymmetric_t>::Zero();
            matrix(i, i) = 3.0;

            CHECK(cabs(average_of_diagonal_of_matrix(matrix) - 1.0) < numerical_tolerance);
            CHECK(cabs(average_of_off_diagonal_of_matrix(matrix)) < numerical_tolerance);
        }
    }
}

} // namespace power_grid_model
