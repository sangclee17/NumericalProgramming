/***************************************************************************
 *   File        : code_4.c
 *   Name        : Sangchul Lee (Scott)
 ***************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define NUM_OF_DATA 4

float solve_lagrange(float* x, float* f_x, float p);
void interval_h(float* h, float* x, int n);
void compute_c1(float c1[][4], float* h, int n);
void compute_c2(float* c2, float* h, float* y, int n);
void gauss_elimination(float A[][4], float* C, int n);
void thomas_algorithm(float* X, float A[][4], float* C, int n);

int main()
{
    float x[4] = { 0, 1, 3, 8 }, f_x[4] = { 0, 11, 17, 17 }, p = 5, h[4], c1[4][4], c2[4], c0[4], sum;
    int n = NUM_OF_DATA;

    float f = solve_lagrange(x, f_x, p);
    printf("\nusing second order Lagrange interpolating polynomials\n");
    printf("f(x=%.0f) = %.04f\n", p, f);

    interval_h(h, x, n);
    compute_c1(c1, h, n);
    compute_c2(c2, h, f_x, n);
    gauss_elimination(c1, c2, n);
    thomas_algorithm(c0, c1, c2, n);

    for (int i = 0; i < n; i++) {
        if (x[i] <= p && p <= x[i + 1]) {
            float a = f_x[i];
            float b = ((f_x[i + 1] - f_x[i]) / h[i]) - (h[i] / 3 * (2 * c0[i] + c0[i + 1]));
            float c = c0[i];
            float d = (c0[i + 1] - c0[i]) / (3 * h[i]);
            sum = a + b * (p - x[i]) + c * pow((p - x[i]), 2) + d * pow((p - x[i]), 3);
        }
    }
    printf("\nusing cubic spline\n");
    printf("f(x=%.0f) = %.04f\n\n", p, sum);

    return 0;
}

float solve_lagrange(float* x, float* f_x, float p)
{
    float l, f = 0;
    for (int i = 0; i < 4; i++) {
        l = 1;
        for (int j = 0; j < 4; j++) {
            if (j != i) {
                l = l * ((p - x[j]) / (x[i] - x[j]));
            }
        }
        f = f + l * f_x[i];
    }
    return f;
}
void interval_h(float* h, float* x, int n)
{
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
    }
}
void compute_c1(float c1[][4], float* h, int n)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                if (i == 0) {
                    c1[i][j] = 1;
                } else if (i == n - 1) {
                    c1[i][j] = 1;
                } else {
                    c1[i][j] = 2 * (h[i - 1] + h[i]);
                    c1[i][i - 1] = h[i - 1];
                    c1[i][i + 1] = h[i];
                }
            }
        }
    }
}
void compute_c2(float* c2, float* h, float* y, int n)
{
    for (int i = 0; i < n; i++) {
        if (i == 0 || i == n - 1) {
            c2[i] = 0;
        } else {
            c2[i] = (3 / h[i] * (y[i + 1] - y[i])) + (3 / h[i - 1] * (y[i - 1] - y[i]));
        }
    }
}
void gauss_elimination(float A[][4], float* C, int n)
{
    for (int i = 1; i < n; i++) {
        A[i][i] = A[i][i] - (A[i][i - 1] * A[i - 1][i]) / A[i - 1][i - 1];
        C[i] = C[i] - (A[i][i - 1] * C[i - 1]) / A[i - 1][i - 1];
        A[i][i - 1] = 0;
    }
}
void thomas_algorithm(float* X, float A[][4], float* C, int n)
{
    for (int i = n - 1; i >= 0; i--) {
        if (i == n - 1) {
            X[i] = C[i] / A[i][i];
        } else {
            X[i] = (C[i] - A[i][i + 1] * X[i + 1]) / A[i][i];
        }
    }
}
