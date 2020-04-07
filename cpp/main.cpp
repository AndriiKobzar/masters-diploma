#include "cmath"
#include "ctime"
#include "integrand_params.h"
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_specfunc.h>
#include <iostream>
#include <mpi.h>
#include <omp.h>

#define REAL(z, i) ((z)[2 * (i)])
#define IMAG(z, i) ((z)[2 * (i) + 1])
#define ABS_ERR 0
#define REL_ERR 0.001
#define WORKSPACE_SIZE 20
#define EPSILON 0.000001
#define NUM_THREADS 8

using namespace std;

double get_random(gsl_rng *r, double sigma) {
    double generatedValue = gsl_ran_gaussian(r, sigma);
    unsigned long min = gsl_rng_min(r);
    unsigned long max = gsl_rng_max(r);
    return (generatedValue);
}

double fbc(double g, int x) {
    return (pow((x + 1), g) + pow(abs(x - 1), g) - 2 * pow(x, g)) / 2;
}

double *get_lambda(double h, int n) {
    int m = 2 * n - 2;
    double *c = new double[2 * m];
    double g = 2 * h;
    for (int i = 0; i < n; ++i) {
        REAL(c, i) = fbc(g, i);
        IMAG(c, i) = 0;
    }
    for (int i = 1; i < n; ++i) {
        REAL(c, m - i) = REAL(c, i);
        IMAG(c, m - i) = 0;
    }

    gsl_fft_complex_radix2_forward(c, 1, m);
    auto result = new double[m];
    for (int i = 0; i < m; ++i) {
        result[i] = sqrt(REAL(c, i));
    }
    delete[] c;
    return result;
}

double *get_fgn(gsl_rng *r, double *lambda, int lambda_len, double t,
                double h) {
    auto random = new double[2 * lambda_len];
    for (int i = 0; i < lambda_len; ++i) {
        REAL(random, i) = get_random(r, 1);
        IMAG(random, i) = 0;
    }
    gsl_fft_complex_radix2_inverse(random, 1, lambda_len);
    for (int i = 0; i < lambda_len; ++i) {
        REAL(random, i) = REAL(random, i) * lambda[i];
    }
    gsl_fft_complex_radix2_forward(random, 1, lambda_len);
    auto result = new double[lambda_len];
    for (int i = 0; i < lambda_len / 2; ++i) {
        result[i] = REAL(random, i) * pow(t / lambda_len, h);
    }
    return result;
}

double sigma(double x) { return pow(sin(x), 2) + 0.05; }

double sigma_der(double x) { return sin(2 * x); }

double sigma_second_der(double x) { return 2 * cos(2 * x); }

double step(double *array, int length, double interval, double x) {
    int index = (int)floor(x / interval);
    if (index >= length) {
        return array[length - 1];
    }
    return array[index];
}

double integral_of_sigma_square(double *y, int length) {
    double sum = 0;
    for (int i = 0; i < length; ++i) {
        sum += pow(sigma(y[i]), 2);
    }
    return sum;
}

double *get_euler_trajectory(double *fbm, int length, double t, double alpha) {
    auto result = new double[length];
    double currentValue = 0;
    double increment = t / length;
    result[0] = currentValue;
    for (int i = 1; i < length; ++i) {
        currentValue =
            currentValue - alpha * currentValue * increment + fbm[i - 1];
        result[i] = currentValue;
    }
    return result;
}

double c(double h) {
    double dividend = 2 * h * gsl_sf_gamma(1.5 - h);
    double divisor = gsl_sf_gamma(h + 0.5) * gsl_sf_gamma(2 - 2 / h);
    return (h - 0.5) * sqrt(dividend / divisor);
}

double dvbyt_integrand(double s, void *p) {
    struct dvbyt_integrand_params *params = (struct dvbyt_integrand_params *)p;
    double h = params->h;
    double u = params->u;
    double alpha = params->alpha;
    double t = params->t;
    double s_minus_u = abs(s - u) <= EPSILON ? 0.001 : s - u;
    double r = c(h) * exp(-alpha * t) * pow(u, 0.5 - h) * pow(s, h - 0.5) *
               exp(alpha * s) * pow(s_minus_u, h - 1.5);
    return r;
}

double dvbyt(double h, double alpha, double u, double t) {
    if (u >= t) {
        return 0;
    }
    gsl_function F;
    F.function = &dvbyt_integrand;
    struct dvbyt_integrand_params params = {h, alpha, u, t};
    F.params = &params;
    gsl_integration_cquad_workspace *w =
        gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    double r;
    gsl_integration_cquad(&F, u, t, ABS_ERR, REL_ERR, w, &r, nullptr, nullptr);
    gsl_integration_cquad_workspace_free(w);
    return r;
}

double s_array_integrand(double x, void *p) {
    struct integrand_params *params = (struct integrand_params *)p;
    double step_y_x = step(params->y, params->n, params->t / params->n, x);
    return sigma(step_y_x) * sigma_der(step_y_x) *
           dvbyt(params->h, params->alpha, params->k * params->t / params->n,
                 x);
}

double *s_array(double *y, int n, double t, double alpha, double h) {
    double *result = new double[n];
    gsl_function F;
    F.function = &s_array_integrand;
    gsl_integration_cquad_workspace *w =
        gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    for (int i = 0; i < n; ++i) {
        struct integrand_params params = {alpha, t, n, h, (i + 1), y};
        F.params = &params;
        double localResult;
        gsl_integration_cquad(&F, 0.001, t, 0, 0.001, w, &localResult, nullptr,
                              nullptr);
        result[i] = 2 * localResult;
        //printf("s[%i] = %f\n", i, result[i]);
    }
    gsl_integration_cquad_workspace_free(w);
    return result;
}

double density_inner_integrand(double v, void *p) {
    struct inner_density_integrand_params *params =
        (struct inner_density_integrand_params *)p;
    return exp(-params->alpha * (params->s - v)) * pow(v, params->h - 0.5) *
           pow(v - params->k * params->t / params->n, params->h - 1.5);
}

double density_outer_integrand(double s, void *p) {
    struct outer_density_integrand_params *params =
        (struct outer_density_integrand_params *)p;
    double tau = params->k * params->t / params->n;
    if (abs(tau - s) < 0.000001) {
        return 0;
    }
    double step_x = step(params->y, params->n, params->t / params->n, s);
    double integrationResult;
    gsl_function integrand;
    integrand.function = &density_inner_integrand;
    inner_density_integrand_params innerDensityIntegrandParams = {
        params->n, params->t, params->alpha, params->h, params->k, s};
    integrand.params = &innerDensityIntegrandParams;

    gsl_integration_cquad_workspace *w =
        gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    gsl_integration_cquad(&integrand, params->k * params->t / params->n, s,
                          ABS_ERR, REL_ERR, w, &integrationResult, nullptr,
                          nullptr);
    gsl_integration_cquad_workspace_free(w);
    return sigma(step_x) * sigma_der(step_x) * integrationResult;
}

double dvbeta_inner_integrand(double q, void *p) {
    struct sndsn_integrand_params *params = (struct sndsn_integrand_params *)p;
    double step_s_q = step(params->s, params->n, params->t / params->n, q);
    double step_y_tau =
        step(params->y, params->n, params->t / params->n, params->tau);
    return -4 * pow(params->eta, 2) * step_s_q *
           (pow(sigma_der(step_y_tau), 2) +
            sigma(step_y_tau) * sigma_second_der(step_y_tau)) *
           dvbyt(params->h, params->alpha, params->v, params->tau) *
           dvbyt(params->h, params->alpha, q, params->tau);
}

double dvbeta_outer_integrand(double tau, void *p) {
    struct sndsn_integrand_params *params = (struct sndsn_integrand_params *)p;
    params->tau = tau;
    // gsl_integration_cquad_workspace *w =
    // gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    double result;
    gsl_function integrand;
    integrand.function = &dvbeta_inner_integrand;
    integrand.params = params;
    double err;
    size_t neval;
    gsl_integration_qng(&integrand, 0.001, params->t, ABS_ERR, REL_ERR, &result,
                        &err, &neval);
    // gsl_integration_cquad_workspace_free(w);
    return result;
}

double dvbeta(double v, sndsn_integrand_params *params) {
    params->v = v;
    // gsl_integration_cquad_workspace *w =
    // gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    double result;
    gsl_function integrand;
    integrand.function = &dvbeta_outer_integrand;
    integrand.params = params;
    double err;
    size_t neval;
    gsl_integration_qng(&integrand, 0.001, params->t, ABS_ERR, REL_ERR, &result,
                        &err, &neval);
    // gsl_integration_cquad_workspace_free(w);
    return result;
}

double sndsn_integrand(double tau, void *p) {
    struct sndsn_integrand_params *params = (struct sndsn_integrand_params *)p;
    return step(params->s, params->n, params->t / params->n, tau) *
           dvbeta(tau, params);
}

double density(double *y, double *s, double *sbm_incr, int n, double t,
               double h, double alpha) {
    double eta = 0;
    for (int i = 0; i < n; ++i) {
        eta += s[i] * s[i];
    }
    eta = 1 / eta;
    double first = 0;
    gsl_integration_cquad_workspace *w =
        gsl_integration_cquad_workspace_alloc(WORKSPACE_SIZE);
    for (int k = 1; k < n; ++k) {
        //printf("first[%i]\n", k);
        double integrationResult;
        gsl_function integrand;
        integrand.function = &density_outer_integrand;
        outer_density_integrand_params params = {y, n, t, alpha, h, k};
        integrand.params = &params;
        gsl_integration_cquad(&integrand, k * t / n, t, ABS_ERR, REL_ERR, w,
                              &integrationResult, nullptr, nullptr);
        first +=
            pow((k + 1) * t / n, 0.5 - h) * integrationResult * sbm_incr[k];
    }
    first *= 2 * eta * c(h);
    double sndsn;
    gsl_function integrand;
    integrand.function = &sndsn_integrand;
    struct sndsn_integrand_params params;
    params.t = t;
    params.alpha = alpha;
    params.n = n;
    params.h = h;
    params.s = s;
    params.y = y;
    params.eta = eta;
    integrand.params = &params;
    double err;
    size_t neval;
    gsl_integration_qng(&integrand, 0.001, t, 2, 2, &sndsn, &err, &neval);
    gsl_integration_cquad_workspace_free(w);
    return first - sndsn;
}

double *sbm_increments(gsl_rng *r, double d, int n) {
    double *result = new double[n];
    for (int i = 0; i < n; ++i) {
        result[i] = get_random(r, d);
    }
    return result;
}

int get_number_of_points_for_node(int numberOfPoints, int worldSize, int rank) {
    int numberOfLocalPoints = (int)(numberOfPoints / worldSize);
    int remainder = numberOfPoints % worldSize;
    double start = rank * (numberOfLocalPoints);
    double end = start + numberOfLocalPoints;
    if (remainder > 0 && rank < remainder) {
        numberOfLocalPoints++;
    }
    return numberOfLocalPoints;
}

int main(int argc, char *argv[]) {
    const int root = 0;
    double h = 0.6;
    int n = 33;
    double t = 1.;
    double alpha = 0.6;
    double globalStart = 0;
    double globalEnd = 2;
    int numberOfPoints = 20;
    const int calculationsPerPoint = 100;
    int errCode;

    if ((errCode = MPI_Init(&argc, &argv)) != 0) {
        return errCode;
    }

    int myRank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    printf("world size: %i\n", worldSize);
    if (myRank == root) {
        printf("MAIN PROCESS\n");
    }

    gsl_set_error_handler_off();
    int quotient = (int)(numberOfPoints / worldSize);
    int remainder = numberOfPoints % worldSize;
    int rangeStart = 0;
    int rangeEnd = 0;
    if (remainder > 0) {
        if (myRank <= remainder) {
            rangeStart = (myRank * (quotient + 1));
            rangeEnd = rangeStart + quotient + 1;
        } else {
            rangeStart = myRank * quotient + remainder;
            rangeEnd = rangeStart + quotient;
        }
    } else {
        rangeStart = myRank * quotient;
        rangeEnd = (myRank + 1) * quotient - 1;
    }

    const long double sysTime = time(nullptr);
    const unsigned long seed = (unsigned long)sysTime + myRank*1000;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(r, seed);
    omp_set_num_threads(NUM_THREADS);
    double* result = new double[rangeEnd - rangeStart + 1];
    for (int i = 0; i <= rangeEnd-rangeStart; i++) {
        double u = globalStart + i * (globalEnd - globalStart) / numberOfPoints;
        double average = 0;
	printf("rank %i point %i\n", myRank, i);
#pragma omp parallel for shared(h, n, t, alpha, u, r)
        for (int j = 0; j < calculationsPerPoint; j++) {
            double *lambda = get_lambda(h, n);
            double *fbm = get_fgn(r, lambda, 2 * n - 2, t, h);
            double *y = get_euler_trajectory(fbm, n, t, alpha);
            double integral = integral_of_sigma_square(y, n);
            if (integral <= u) {
                printf("density(%f)=%f\n", u, .0);
                continue;
            }

            double *s = s_array(y, n, t, alpha, h);
            double *sbm_incr = sbm_increments(r, t / n, n);
            double densityValue = density(y, s, sbm_incr, n, t, h, alpha);
            printf("density(%f)=%f", u, densityValue);
            delete[] lambda;
            delete[] fbm;
            delete[] y;
            delete[] s;
            delete[] sbm_incr;
            #pragma omp atomic
            average += densityValue / numberOfPoints;
        }
        result[i] = average;
    }
    gsl_rng_free(r);
    int* displacements;
    double* recv;
    int* recvCounts;
    if(myRank == 0) {
        recv = new double[numberOfPoints];
        recvCounts = new int[worldSize];
        displacements = new int[worldSize];
        for (int i = 0; i < worldSize; i++)
        {
            if(i==0){
                displacements[i] = 0;
		recvCounts[i] = quotient + (remainder>0 ? 1 : 0);
		continue;
            } else {
                displacements[i] = displacements[i-1] + quotient;
            }
            recvCounts[i] = quotient;
            if(i<remainder) {
                displacements[i]++;
                recvCounts[i]++;
            } else if(i==remainder) {
                displacements[i]++;
            }
        }
        
    }
    MPI_Gatherv(result, rangeEnd - rangeStart + 1, MPI_DOUBLE, recv, recvCounts, displacements, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if(myRank == root) {
        for (int i = 0; i < numberOfPoints; i++)
        {
            printf("delta[%i]=%f",i, recv[i]);
        }
    }
    MPI_Finalize();
    return 0;
}
