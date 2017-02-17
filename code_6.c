/***************************************************************************
 *   File        : code_6.c
 *   Name        : Sangchul Lee (Scott)
 ***************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265359

void compute_fi_n(double* fi_n, int n);
FILE* safe_open_file(const char* file_name, const char* mode);
double RHS(double* f, int n, int x);
void compute_fi_1(double* fi_1, double* fi_n, int t, int n);
void compute_new_fi_n(double* fi_1, double* fi_n, int t, int n);

int main()
{
    int n, t;

    printf("resolution of x : ");
    scanf("%d", &n);
    printf("resolution of t: ");
    scanf("%d", &t);

    double fi_n[n + 1], fi_1[n + 1];

    FILE* fp = safe_open_file("wave.csv", "w");
    compute_fi_n(fi_n, n);

    for (int i = 0; i <= n; i++) {
        fprintf(fp, "%f,%f\n", (double)i / n, fi_n[i]);
    }

    for (int j = 1; j <= t; j++) {
        compute_fi_1(fi_1, fi_n, t, n);
        compute_new_fi_n(fi_1, fi_n, t, n);
        for (int i = 0; i <= n; i++) {
            fprintf(fp, "%f,%f\n", (double)i / n, fi_n[i]);
        }
    }
    fclose(fp);
    return 0;
}

void compute_fi_n(double* fi_n, int n)
{
    for (int i = 0; i <= n; i++) {
        double del_x_i = (double)1 / n * i;

        if (del_x_i < 0.125) {
            fi_n[i] = 0.0;
        } else if (del_x_i <= 0.375) {
            fi_n[i] = 0.5 * (1 - cos(8 * PI * (del_x_i - 0.125)));
        } else if (del_x_i <= 1.0) {
            fi_n[i] = 0.0;
        }
    }
}

FILE* safe_open_file(const char* file_name, const char* mode)
{
    FILE* fp = fopen(file_name, mode);
    if (fp == NULL) {
        perror("file open error.");
        exit(EXIT_FAILURE);
    }
    return fp;
}

double RHS(double* f, int n, int x)
{
    double x_resolution = (double)1 / n;
    /* upwind scheme*/
    if (x == 0) {
        return (-1.0) * (f[1] - f[0]) / x_resolution;
    } else {
        return (-1.0) * (f[x] - f[x - 1]) / x_resolution;
    }
    /* central scheme 
     if (x == 0) {
        return (-1.0) * (f[1] - f[0]) / x_resolution;
    } else if (x == n) {
        return (-1.0) * (f[n] - f[n - 1]) / x_resolution;
    } else {
        return (-1.0) * (f[x + 1] - f[x - 1]) / 2.0 / x_resolution;
    }*/  
}

void compute_fi_1(double* fi_1, double* fi_n, int t, int n)
{
    double t_resolution = (double)0.2 / t;

    for (int i = 0; i <= n; i++) {
        fi_1[i] = fi_n[i] + t_resolution * RHS(fi_n, n, i);
    }
}
void compute_new_fi_n(double* fi_1, double* fi_n, int t, int n)
{
    double t_resolution = (double)0.2 / t;

    for (int i = 0; i <= n; i++) {
        fi_n[i] = fi_n[i] + 1.0 / 2.0 * t_resolution * (RHS(fi_n, n, i) + RHS(fi_1, n, i));
    }
}
