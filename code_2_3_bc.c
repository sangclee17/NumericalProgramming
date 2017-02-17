/***************************************************************************
 *   File        : code_2_3_bc.c
 *   Name        : Sangchul Lee (Scott)
 ***************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOT_CONVERGING (-1)
#define PI 3.14159265359

double f(double beta, double theta, double M);
double fprime(double beta, double M);
double newtonraphs(double beta, double theta, double M);
FILE* safe_open_file(const char* file_name, const char* mode);

int main(int argc, char* argv[])
{
    double beta;
    double LEFT;
    double M;

    FILE* fp = safe_open_file("newton_raphs.csv", "a");

    while (1) {

        printf("\nEnter 0 to terminate!!\n");
        printf("Enter Mach Number M you want to evaluate: ");
        scanf("%la", &M);

        if (M == 0) {
            break;
        }
        printf("Enter Initial Guess: ");
        scanf("%la", &LEFT);

        for (int theta = 0; theta < 45; theta++) {
            beta = newtonraphs(LEFT, theta, M);
            if ((beta <= PI / 2.0) && (beta > 0)) {
                fprintf(fp, "%d,%.10f\n", theta, beta * 180 / PI);
            }
        }
    }

    fclose(fp);

    return 0;
}

double f(double beta, double theta, double M)
{
    double gamma = 1.4;
    double thetarad = theta * (PI / 180);

    return (2 * cos(beta) / sin(beta) * ((pow(M * sin(beta), 2) - 1) / (pow(M, 2) * (gamma + cos(2 * beta)) + 2))) - tan(thetarad);
}

double fprime(double beta, double M)
{
    double gamma = 1.4;

    return (4 * pow(M, 2) * cos(beta) * cos(beta) / sin(beta) * sin(beta)) / ((gamma + cos(2 * beta)) * pow(M, 2) + 2) - (2 * (pow((cos(beta) / sin(beta)), 2) + 1) * (pow(M, 2) * pow(sin(beta), 2) - 1)) / ((gamma + cos(2 * beta)) * pow(M, 2) + 2) + (4 * pow(M, 2) * sin(2 * beta) * cos(beta) / sin(beta) * (pow(M, 2) * pow(sin(beta), 2) - 1)) / pow((pow(M, 2) * (gamma + cos(2 * beta)) + 2), 2);
}

double newtonraphs(double beta, double theta, double M)
{
    double fx1, fpx1;
    int iterations = 0;
    double eps = 1e-4;
    int limit = 10000;

    fx1 = f(beta, theta, M);
    fpx1 = fprime(beta, M);

    while (fabs(fx1) > eps) {
        iterations = iterations + 1;
        if (iterations == limit) {
            exit(NOT_CONVERGING);
        }
        beta = beta - fx1 / fpx1;
        fx1 = f(beta, theta, M);
        fpx1 = fprime(beta, M);
    }
    return beta = beta - fx1 / fpx1;
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
