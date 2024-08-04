#include <check.h>

#include "s21_matrix.h"

START_TEST(s21_create_matrix_1) {
  matrix_t A = {};
  ck_assert_int_eq(s21_create_matrix(5, 5, &A), OK);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_create_matrix_2) {
  ck_assert_int_eq(s21_create_matrix(5, 5, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_create_matrix_3) {
  matrix_t A = {};
  ck_assert_int_eq(s21_create_matrix(5, 0, &A), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_remove_matrix_1) {
  matrix_t A = {};
  s21_create_matrix(5, 5, &A);
  s21_remove_matrix(&A);
  ck_assert_ptr_null(A.matrix);
}
END_TEST

START_TEST(s21_remove_matrix_2) {
  s21_remove_matrix(NULL);
  ck_assert_int_eq(1, 1);
}
END_TEST

START_TEST(s21_remove_matrix_3) {
  matrix_t A = {};
  s21_remove_matrix(&A);
  ck_assert_int_eq(1, 1);
}
END_TEST

START_TEST(s21_eq_matrix_1) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_initialize_matrix(&A, 21, 21);
  s21_initialize_matrix(&B, 21, 21);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_2) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_3) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_initialize_matrix(&A, 1, 0.000000001);
  s21_initialize_matrix(&B, 1, 0.000000002);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_4) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(3, 4, &B);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_eq_matrix_5) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_eq_matrix_6) {
  matrix_t A = {};
  matrix_t B = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_initialize_matrix(&A, 21, 1);
  s21_initialize_matrix(&B, 42, 1);
  ck_assert_int_eq(s21_eq_matrix(&A, &B), FAILURE);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_sum_matrix_1) {
  matrix_t A = {};
  matrix_t B = {};
  ck_assert_int_eq(s21_sum_matrix(&A, &B, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_sum_matrix_2) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(5, 5, &B);

  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_sum_matrix_4) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &eq_matrix);
  s21_initialize_matrix(&A, 1, 1);
  s21_initialize_matrix(&B, 1, 1);
  s21_initialize_matrix(&eq_matrix, 2, 2);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_sum_matrix_5) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &eq_matrix);
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_sum_matrix_6) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 5;
  A.matrix[2][0] = 0, A.matrix[2][1] = 0, A.matrix[2][2] = 6;
  B.matrix[0][0] = 1, B.matrix[0][1] = 0, B.matrix[0][2] = 0;
  B.matrix[1][0] = 2, B.matrix[1][1] = 0, B.matrix[1][2] = 0;
  B.matrix[2][0] = 3, B.matrix[2][1] = 4, B.matrix[2][2] = 1;
  ck_assert_int_eq(s21_sum_matrix(&A, &B, &result), OK);
  eq_matrix.matrix[0][0] = 2, eq_matrix.matrix[0][1] = 2,
  eq_matrix.matrix[0][2] = 3;
  eq_matrix.matrix[1][0] = 2, eq_matrix.matrix[1][1] = 4,
  eq_matrix.matrix[1][2] = 5;
  eq_matrix.matrix[2][0] = 3, eq_matrix.matrix[2][1] = 4,
  eq_matrix.matrix[2][2] = 7;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_sub_matrix_1) {
  matrix_t A = {};
  matrix_t B = {};
  ck_assert_int_eq(s21_sub_matrix(&A, &B, NULL), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_sub_matrix_2) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  s21_create_matrix(4, 4, &A);
  s21_create_matrix(5, 5, &B);
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_sub_matrix_4) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &eq_matrix);
  s21_initialize_matrix(&A, 1, 1);
  s21_initialize_matrix(&B, 1, 1);
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_sub_matrix_5) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(5, 5, &A);
  s21_create_matrix(5, 5, &B);
  s21_create_matrix(5, 5, &eq_matrix);
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_sub_matrix_6) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 5;
  A.matrix[2][0] = 0, A.matrix[2][1] = 0, A.matrix[2][2] = 6;
  B.matrix[0][0] = 1, B.matrix[0][1] = 0, B.matrix[0][2] = 0;
  B.matrix[1][0] = 2, B.matrix[1][1] = 0, B.matrix[1][2] = 0;
  B.matrix[2][0] = 3, B.matrix[2][1] = 4, B.matrix[2][2] = 1;
  ck_assert_int_eq(s21_sub_matrix(&A, &B, &result), OK);
  eq_matrix.matrix[0][0] = 0, eq_matrix.matrix[0][1] = 2,
  eq_matrix.matrix[0][2] = 3;
  eq_matrix.matrix[1][0] = -2, eq_matrix.matrix[1][1] = 4,
  eq_matrix.matrix[1][2] = 5;
  eq_matrix.matrix[2][0] = -3, eq_matrix.matrix[2][1] = -4,
  eq_matrix.matrix[2][2] = 5;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_number_1) {
  matrix_t A = {};
  double number = 3.14;
  s21_create_matrix(3, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_mult_number(&A, number, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_mult_number_4) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  double number = 3;
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  ck_assert_int_eq(s21_mult_number(&A, number, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_number_5) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  double number = 3;
  s21_create_matrix(3, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  s21_create_matrix(3, 3, &eq_matrix);
  s21_initialize_matrix(&eq_matrix, 3, 3);
  ck_assert_int_eq(s21_mult_number(&A, number, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_number_6) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  double number = 2;
  s21_create_matrix(3, 3, &A);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 2;
  A.matrix[2][0] = 2, A.matrix[2][1] = 3, A.matrix[2][2] = 4;
  s21_create_matrix(3, 3, &eq_matrix);
  eq_matrix.matrix[0][0] = 2, eq_matrix.matrix[0][1] = 4,
  eq_matrix.matrix[0][2] = 6;
  eq_matrix.matrix[1][0] = 0, eq_matrix.matrix[1][1] = 8,
  eq_matrix.matrix[1][2] = 4;
  eq_matrix.matrix[2][0] = 4, eq_matrix.matrix[2][1] = 6,
  eq_matrix.matrix[2][2] = 8;
  ck_assert_int_eq(s21_mult_number(&A, number, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_matrix_1) {
  matrix_t A = {};
  s21_create_matrix(3, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_mult_matrix(&A, NULL, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_mult_matrix_3) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(3, 2, &B);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}
END_TEST

START_TEST(s21_mult_matrix_4) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &B);
  s21_create_matrix(3, 3, &eq_matrix);
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_matrix_5) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  s21_create_matrix(3, 3, &B);
  s21_initialize_matrix(&B, 1, 1);
  s21_create_matrix(3, 3, &eq_matrix);
  eq_matrix.matrix[0][0] = 30, eq_matrix.matrix[0][1] = 36,
  eq_matrix.matrix[0][2] = 42;
  eq_matrix.matrix[1][0] = 66, eq_matrix.matrix[1][1] = 81,
  eq_matrix.matrix[1][2] = 96;
  eq_matrix.matrix[2][0] = 102, eq_matrix.matrix[2][1] = 126,
  eq_matrix.matrix[2][2] = 150;
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_mult_matrix_6) {
  matrix_t A = {};
  matrix_t B = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 2, &A);
  A.matrix[0][0] = 1, A.matrix[0][1] = 4;
  A.matrix[1][0] = 2, A.matrix[1][1] = 5;
  A.matrix[2][0] = 3, A.matrix[2][1] = 6;
  s21_create_matrix(2, 3, &B);
  B.matrix[0][0] = 1, B.matrix[0][1] = -1, B.matrix[0][2] = 1;
  B.matrix[1][0] = 2, B.matrix[1][1] = 3, B.matrix[1][2] = 4;
  s21_create_matrix(3, 3, &eq_matrix);
  eq_matrix.matrix[0][0] = 9, eq_matrix.matrix[0][1] = 11,
  eq_matrix.matrix[0][2] = 17;
  eq_matrix.matrix[1][0] = 12, eq_matrix.matrix[1][1] = 13,
  eq_matrix.matrix[1][2] = 22;
  eq_matrix.matrix[2][0] = 15, eq_matrix.matrix[2][1] = 15,
  eq_matrix.matrix[2][2] = 27;
  ck_assert_int_eq(s21_mult_matrix(&A, &B, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_transpose_1) {
  matrix_t A = {};
  matrix_t result = {};
  ck_assert_int_eq(s21_transpose(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_transpose_2) {
  matrix_t A = {};
  matrix_t result = {};
  s21_create_matrix(3, 2, &A);
  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 2);
  ck_assert_int_eq(result.columns, 3);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_transpose_3) {
  matrix_t A = {};
  matrix_t result = {};
  s21_create_matrix(2, 3, &A);
  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_transpose_4) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(2, 3, &A);
  s21_create_matrix(3, 2, &eq_matrix);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  eq_matrix.matrix[0][0] = 1, eq_matrix.matrix[0][1] = 4;
  eq_matrix.matrix[1][0] = 2, eq_matrix.matrix[1][1] = 5;
  eq_matrix.matrix[2][0] = 3, eq_matrix.matrix[2][1] = 6;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  ck_assert_int_eq(result.rows, 3);
  ck_assert_int_eq(result.columns, 2);
  s21_remove_matrix(&A);
  s21_remove_matrix(&eq_matrix);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_transpose_5) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 2, &A);
  s21_create_matrix(2, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 4;
  A.matrix[1][0] = 2, A.matrix[1][1] = 5;
  A.matrix[2][0] = 3, A.matrix[2][1] = 6;
  ck_assert_int_eq(s21_transpose(&A, &result), OK);
  eq_matrix.matrix[0][0] = 1, eq_matrix.matrix[0][1] = 2,
  eq_matrix.matrix[0][2] = 3;
  eq_matrix.matrix[1][0] = 4, eq_matrix.matrix[1][1] = 5,
  eq_matrix.matrix[1][2] = 6;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  ck_assert_int_eq(result.rows, 2);
  ck_assert_int_eq(result.columns, 3);
  s21_remove_matrix(&A);
  s21_remove_matrix(&eq_matrix);
  s21_remove_matrix(&result);
}
END_TEST

START_TEST(s21_calc_complements_1) {
  matrix_t A = {};
  s21_create_matrix(3, 3, &A);
  ck_assert_int_eq(s21_calc_complements(&A, NULL), INCORRECT_MATRIX);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_calc_complements_2) {
  matrix_t A = {};
  matrix_t result = {};
  s21_create_matrix(1, 3, &A);
  s21_initialize_matrix(&A, 1, 3);
  ck_assert_int_eq(s21_calc_complements(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_calc_complements_3) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 1, A.matrix[0][1] = 2, A.matrix[0][2] = 3;
  A.matrix[1][0] = 0, A.matrix[1][1] = 4, A.matrix[1][2] = 2;
  A.matrix[2][0] = 5, A.matrix[2][1] = 2, A.matrix[2][2] = 1;
  ck_assert_int_eq(s21_calc_complements(&A, &result), OK);
  eq_matrix.matrix[0][0] = 0, eq_matrix.matrix[0][1] = 10,
  eq_matrix.matrix[0][2] = -20;
  eq_matrix.matrix[1][0] = 4, eq_matrix.matrix[1][1] = -14,
  eq_matrix.matrix[1][2] = 8;
  eq_matrix.matrix[2][0] = -8, eq_matrix.matrix[2][1] = -2,
  eq_matrix.matrix[2][2] = 4;
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_determinant_1) {
  matrix_t A = {};
  double det = 0;
  ck_assert_int_eq(s21_determinant(&A, &det), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_determinant_2) {
  matrix_t A = {};
  double det = 0;
  s21_create_matrix(3, 2, &A);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_determinant(&A, &det), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_3) {
  matrix_t A = {};
  double det = 0;
  s21_create_matrix(3, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_determinant(&A, &det), OK);
  ck_assert_double_eq(det, 0);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_4) {
  matrix_t A = {};
  double det = 0;
  s21_create_matrix(2, 2, &A);
  s21_initialize_matrix(&A, 3, 3);
  ck_assert_int_eq(s21_determinant(&A, &det), OK);
  ck_assert_double_eq(det, -18);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_determinant_5) {
  matrix_t A = {};
  double det = 0;
  s21_create_matrix(1, 1, &A);
  A.matrix[0][0] = 21;
  ck_assert_int_eq(s21_determinant(&A, &det), OK);
  ck_assert_double_eq(det, 21);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_1) {
  matrix_t A = {};
  matrix_t result = {};
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), INCORRECT_MATRIX);
}
END_TEST

START_TEST(s21_inverse_matrix_3) {
  matrix_t A = {};
  matrix_t result = {};
  s21_create_matrix(5, 3, &A);
  s21_initialize_matrix(&A, 1, 1);
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), CALCULATION_ERROR);
  s21_remove_matrix(&A);
}
END_TEST

START_TEST(s21_inverse_matrix_4) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(1, 1, &A);
  s21_create_matrix(1, 1, &eq_matrix);
  A.matrix[0][0] = 21;
  eq_matrix.matrix[0][0] = 1.0 / 21.0;
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

START_TEST(s21_inverse_matrix_5) {
  matrix_t A = {};
  matrix_t result = {};
  matrix_t eq_matrix = {};
  s21_create_matrix(3, 3, &A);
  s21_create_matrix(3, 3, &eq_matrix);
  A.matrix[0][0] = 2, A.matrix[0][1] = 5, A.matrix[0][2] = 7;
  A.matrix[1][0] = 6, A.matrix[1][1] = 3, A.matrix[1][2] = 4;
  A.matrix[2][0] = 5, A.matrix[2][1] = -2, A.matrix[2][2] = -3;
  eq_matrix.matrix[0][0] = 1, eq_matrix.matrix[0][1] = -1,
  eq_matrix.matrix[0][2] = 1;
  eq_matrix.matrix[1][0] = -38, eq_matrix.matrix[1][1] = 41,
  eq_matrix.matrix[1][2] = -34;
  eq_matrix.matrix[2][0] = 27, eq_matrix.matrix[2][1] = -29,
  eq_matrix.matrix[2][2] = 24;
  ck_assert_int_eq(s21_inverse_matrix(&A, &result), OK);
  ck_assert_int_eq(s21_eq_matrix(&result, &eq_matrix), SUCCESS);
  s21_remove_matrix(&A);
  s21_remove_matrix(&result);
  s21_remove_matrix(&eq_matrix);
}
END_TEST

Suite *s21_matrix_suite(void) {
  Suite *s;
  TCase *tc_core;
  s = suite_create("s21_matrix");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, s21_create_matrix_1);
  tcase_add_test(tc_core, s21_create_matrix_2);
  tcase_add_test(tc_core, s21_create_matrix_3);

  tcase_add_test(tc_core, s21_remove_matrix_1);
  tcase_add_test(tc_core, s21_remove_matrix_2);
  tcase_add_test(tc_core, s21_remove_matrix_3);

  tcase_add_test(tc_core, s21_eq_matrix_1);
  tcase_add_test(tc_core, s21_eq_matrix_2);
  tcase_add_test(tc_core, s21_eq_matrix_3);
  tcase_add_test(tc_core, s21_eq_matrix_4);
  tcase_add_test(tc_core, s21_eq_matrix_5);
  tcase_add_test(tc_core, s21_eq_matrix_6);

  tcase_add_test(tc_core, s21_sum_matrix_1);
  tcase_add_test(tc_core, s21_sum_matrix_2);
  tcase_add_test(tc_core, s21_sum_matrix_4);
  tcase_add_test(tc_core, s21_sum_matrix_5);
  tcase_add_test(tc_core, s21_sum_matrix_6);

  tcase_add_test(tc_core, s21_sub_matrix_1);
  tcase_add_test(tc_core, s21_sub_matrix_2);
  tcase_add_test(tc_core, s21_sub_matrix_4);
  tcase_add_test(tc_core, s21_sub_matrix_5);
  tcase_add_test(tc_core, s21_sub_matrix_6);

  tcase_add_test(tc_core, s21_mult_number_1);
  tcase_add_test(tc_core, s21_mult_number_4);
  tcase_add_test(tc_core, s21_mult_number_5);
  tcase_add_test(tc_core, s21_mult_number_6);

  tcase_add_test(tc_core, s21_mult_matrix_1);
  tcase_add_test(tc_core, s21_mult_matrix_3);
  tcase_add_test(tc_core, s21_mult_matrix_4);
  tcase_add_test(tc_core, s21_mult_matrix_5);
  tcase_add_test(tc_core, s21_mult_matrix_6);

  tcase_add_test(tc_core, s21_transpose_1);
  tcase_add_test(tc_core, s21_transpose_2);
  tcase_add_test(tc_core, s21_transpose_3);
  tcase_add_test(tc_core, s21_transpose_4);
  tcase_add_test(tc_core, s21_transpose_5);

  tcase_add_test(tc_core, s21_calc_complements_1);
  tcase_add_test(tc_core, s21_calc_complements_2);
  tcase_add_test(tc_core, s21_calc_complements_3);

  tcase_add_test(tc_core, s21_determinant_1);
  tcase_add_test(tc_core, s21_determinant_2);
  tcase_add_test(tc_core, s21_determinant_3);
  tcase_add_test(tc_core, s21_determinant_4);
  tcase_add_test(tc_core, s21_determinant_5);

  tcase_add_test(tc_core, s21_inverse_matrix_1);
  tcase_add_test(tc_core, s21_inverse_matrix_3);
  tcase_add_test(tc_core, s21_inverse_matrix_4);
  tcase_add_test(tc_core, s21_inverse_matrix_5);

  suite_add_tcase(s, tc_core);
  return s;
}

/*  TESTS MAIN */

int main(void) {
  int number_failed;
  Suite *s;
  SRunner *sr;

  s = s21_matrix_suite();
  sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? 0 : 1;
}