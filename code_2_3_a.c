/***************************************************************************
 *   File        : code_2_3_a.c
 *   Name        : Sangchul Lee (Scott)
 ***************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOT_CONVERGING (-1)
#define PI 3.14159265359

double f(double beta, double theta);
double fprime(double beta);
double newtonraphs(double beta, double theta, double eps, int limit);

int main(int argc, char* argv[])
{
    double beta;
    double LEFT;
    int theta = 20;

    printf("Enter Initial Guess: ");
    scanf("%lf", &LEFT);

    beta = newtonraphs(LEFT, theta, 1e-4, 1000);
    printf("(main) beta = %.10f\n", beta * 180 / PI);

    return 0;
}

double f(double beta, double theta)
{
    double gamma = 1.4;
    double M = 5.0;
    double thetarad = theta * (PI / 180);

    return (2 * cos(beta) / sin(beta) * ((pow(M * sin(beta), 2) - 1) / (pow(M, 2) * (gamma + cos(2 * beta)) + 2))) - tan(thetarad);
}

double fprime(double beta)
{
    double gamma = 1.4;
    double M = 5.0;

    return (4 * pow(M, 2) * cos(beta) * cos(beta) / sin(beta) * sin(beta)) / ((gamma + cos(2 * beta)) * pow(M, 2) + 2) - (2 * (pow((cos(beta) / sin(beta)), 2) + 1) * (pow(M, 2) * pow(sin(beta), 2) - 1)) / ((gamma + cos(2 * beta)) * pow(M, 2) + 2) + (4 * pow(M, 2) * sin(2 * beta) * cos(beta) / sin(beta) * (pow(M, 2) * pow(sin(beta), 2) - 1)) / pow((pow(M, 2) * (gamma + cos(2 * beta)) + 2), 2);
}

double newtonraphs(double beta, double theta, double eps, int limit)
{
    double fx1, fpx1;
    int iterations = 0;

    fx1 = f(beta, theta);
    fpx1 = fprime(beta);

    while (fabs(fx1) > eps) {
        iterations = iterations + 1;
        if (iterations == limit) {
            exit(NOT_CONVERGING);
        }
        beta = beta - fx1 / fpx1;
        fx1 = f(beta, theta);
        fpx1 = fprime(beta);
    }

    return beta = beta - fx1 / fpx1;
}
