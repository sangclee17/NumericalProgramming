/***************************************************************************
 *   File        : code_3.c
 *   Name        : Sangchul Lee (Scott)
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#define MATRIX_SIZE 5

void gauss_elimination(float A[][5], float* C, int n);
void solve_for_X(float A[][5], float* C, float* X, int n);

int main()
{

    float A[5][5] = { { 1, -2, 0, 0, 0 }, { 2, 4, 5, 0, 0 }, { 0, 8, -9, 2, 0 }, { 0, 0, 6, 3, 4 }, { 0, 0, 0, 6, 3 } },
          C[5] = { 1, 2, -3, 4, 5 }, X[5];
    int n = MATRIX_SIZE;

    gauss_elimination(A, C, n);
    printf("\nAfter Gauss Elimination\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.0f\t", A[i][j]);
        }
        printf("|  %.0f\n", C[i]);
    }
    printf("\nsolution\n");
    solve_for_X(A, C, X, n);
    for (int i = 0; i < n; i++) {
        printf("X%d = %.04f\n", i, X[i]);
    }

    return 0;
}
void gauss_elimination(float A[][5], float* C, int n)
{
    for (int i = 1; i < n; i++) {
        A[i][i] = A[i][i] - (A[i][i - 1] * A[i - 1][i]) / A[i - 1][i - 1];
        C[i] = C[i] - (A[i][i - 1] * C[i - 1]) / A[i - 1][i - 1];
        A[i][i - 1] = 0;
    }
}

void solve_for_X(float A[][5], float* C, float* X, int n)
{
    for (int i = n - 1; i >= 0; i--) {
        if (i == n - 1) {
            X[i] = C[i] / A[i][i];
        } else {
            X[i] = (C[i] - A[i][i + 1] * X[i + 1]) / A[i][i];
        }
    }
}
