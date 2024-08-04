#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int matrix_status = OK;
  if (result == NULL || rows < 1 || columns < 1) {
    matrix_status = INCORRECT_MATRIX;
  } else {
    result->rows = rows;
    result->columns = columns;
    result->matrix = (double **)calloc(rows + rows * columns, sizeof(double));
    if (result->matrix == NULL) {
      matrix_status = INCORRECT_MATRIX;
    } else {
      double *starting_point = (double *)(result->matrix + rows);
      int iteration_limit = rows <= columns ? columns : rows;
      for (int i = 0; i < iteration_limit; i++)
        result->matrix[i] = starting_point + i * columns;
    }
  }
  return matrix_status;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    free(A->matrix);
    A->columns = 0;
    A->rows = 0;
    A->matrix = NULL;
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int matrix_status = SUCCESS;
  matrix_status = s21_check_matrix(A, B);
  if (matrix_status) matrix_status = s21_check_matrix_equal_size(A, B);
  if (matrix_status) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-07)
          matrix_status = FAILURE;
      }
    }
  }
  return matrix_status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int code = 1;
  return sum_sub_func(A, B, result, code);
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int code = 2;
  return sum_sub_func(A, B, result, code);
}

int sum_sub_func(matrix_t *A, matrix_t *B, matrix_t *result, int code) {
  int matrix_status = OK;
  if (s21_check_matrix(A, B) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;
  if (matrix_status == OK && (A->rows != B->rows || A->columns != B->columns))
    matrix_status = CALCULATION_ERROR;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->rows, B->columns, result);
  if (!matrix_status) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < B->columns; j++) {
        if (code == 1) result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        if (code == 2) result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }
  return matrix_status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int matrix_status = OK;
  if (s21_check_matrix(A, B) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;
  if (matrix_status == OK && s21_check_rows_columns(A, B) == FAILURE)
    matrix_status = CALCULATION_ERROR;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->rows, B->columns, result);
  if (!matrix_status) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        for (int k = 0; k < B->rows; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }
  return matrix_status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int matrix_status = OK;
  if (s21_check_single_matrix(A) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->rows, A->columns, result);
  if (!matrix_status) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return matrix_status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int matrix_status = OK;
  if (s21_check_single_matrix(A) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->columns, A->rows, result);
  if (!matrix_status) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return matrix_status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int matrix_status = OK;
  if (s21_check_single_matrix(A) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;
  if (!matrix_status &&
      (s21_matrix_is_square(A) == FAILURE || (A->rows == 1 || A->columns == 1)))
    matrix_status = CALCULATION_ERROR;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->rows, A->columns, result);
  if (!matrix_status) {
    for (int i = 0; i < A->rows && !matrix_status; i++) {
      int sign = i % 2 == 0 ? 1 : -1;
      for (int j = 0; j < A->columns && !matrix_status; j++) {
        matrix_t *minor = s21_minor_create(i, j, A);
        if (minor == NULL) {
          matrix_status = INCORRECT_MATRIX;
        } else {
          double det = 0;
          matrix_status = s21_determinant(minor, &det);
          if (!matrix_status) {
            result->matrix[i][j] = sign * det;
            sign = -sign;
          }
          s21_remove_matrix(minor);
          free(minor);
          minor = NULL;
        }
      }
    }
  }
  return matrix_status;
}

matrix_t *s21_minor_create(int excluded_row, int excluded_column, matrix_t *A) {
  int matrix_status = OK;
  matrix_t *minor = calloc(1, sizeof(matrix_t));
  if (s21_check_single_matrix(A) == FAILURE) matrix_status = INCORRECT_MATRIX;

  if (!matrix_status)
    matrix_status = s21_create_matrix(A->rows - 1, A->columns - 1, minor);
  if (!matrix_status && minor != NULL) {
    for (int i = 0, row = 0; i < A->rows; i++) {
      for (int j = 0, column = 0; j < A->columns; j++) {
        if (i != excluded_row && j != excluded_column)
          minor->matrix[row][column++] = A->matrix[i][j];
      }
      if (i != excluded_row) row++;
    }
  }
  return matrix_status == INCORRECT_MATRIX ? NULL : minor;
}

int s21_determinant(matrix_t *A, double *result) {
  int matrix_status = OK;
  int sign = 1;
  if (s21_check_single_matrix(A) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;
  if (!matrix_status && s21_matrix_is_square(A) == FAILURE)
    matrix_status = CALCULATION_ERROR;

  if (!matrix_status) {
    if (A->rows == 1) *result = A->matrix[0][0];
    if (A->rows == 2)
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
    if (A->rows > 2) {
      *result = 0;
      for (int i = 0; i < A->columns && !matrix_status; i++) {
        matrix_t *minor = s21_minor_create(0, i, A);
        if (minor == NULL) matrix_status = INCORRECT_MATRIX;
        if (minor != NULL) {
          double minor_det = 0;
          matrix_status = s21_determinant(minor, &minor_det);
          if (!matrix_status) {
            *result += sign * A->matrix[0][i] * minor_det;
            sign = -sign;
          }
          s21_remove_matrix(minor);
          free(minor);
          minor = NULL;
        }
      }
    }
  }
  return matrix_status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int matrix_status = OK;
  double det = 0;
  if (s21_check_single_matrix(A) == FAILURE || result == NULL)
    matrix_status = INCORRECT_MATRIX;
  if (!matrix_status && (s21_matrix_is_square(A) == FAILURE))
    matrix_status = CALCULATION_ERROR;
  if (!matrix_status && s21_determinant(A, &det) && det == 0)
    matrix_status = CALCULATION_ERROR;

  if (!matrix_status && det != 0) {
    if (A->rows == 1) {
      matrix_status = s21_create_matrix(A->rows, A->columns, result);
      if (!matrix_status) result->matrix[0][0] = 1 / A->matrix[0][0];
    } else {
      matrix_t tmp_calc_complement = {0}, tmp_transpose = {0};
      matrix_status = s21_calc_complements(A, &tmp_calc_complement);
      if (!matrix_status)
        matrix_status = s21_transpose(&tmp_calc_complement, &tmp_transpose);
      if (!matrix_status)
        matrix_status = s21_mult_number(&tmp_transpose, 1 / det, result);
      s21_remove_matrix(&tmp_calc_complement);
      s21_remove_matrix(&tmp_transpose);
    }
  }
  if (!matrix_status && det == 0) matrix_status = CALCULATION_ERROR;
  return matrix_status;
}

int s21_check_matrix(matrix_t *A, matrix_t *B) {
  int matrix_status = SUCCESS;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0)
    matrix_status = FAILURE;
  if (B == NULL || B->matrix == NULL || B->rows <= 0 || B->columns <= 0)
    matrix_status = FAILURE;
  return matrix_status;
}

int s21_check_matrix_equal_size(matrix_t *A, matrix_t *B) {
  return (A->rows == B->rows && A->columns == B->columns) ? SUCCESS : FAILURE;
}

int s21_check_rows_columns(matrix_t *A, matrix_t *B) {
  return A->columns == B->rows ? SUCCESS : FAILURE;
}

int s21_check_single_matrix(matrix_t *A) {
  int matrix_status = SUCCESS;
  if (A == NULL || A->matrix == NULL || A->rows <= 0 || A->columns <= 0)
    matrix_status = FAILURE;
  return matrix_status;
}

int s21_matrix_is_square(matrix_t *A) {
  int matrix_status;
  if (A != NULL) {
    matrix_status = A->rows == A->columns ? SUCCESS : FAILURE;
  }
  return matrix_status;
}

void s21_initialize_matrix(matrix_t *A, double start_value,
                           double iteration_step) {
  if (A != NULL && A->matrix != NULL) {
    double value = start_value;
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        A->matrix[i][j] = value;
        value += iteration_step;
      }
    }
  }
}