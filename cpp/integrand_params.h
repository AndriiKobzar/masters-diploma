//
// Created by Andriy Kobzar on 22.03.2020.
//

#ifndef CPP_INTEGRAND_PARAMS_H
#define CPP_INTEGRAND_PARAMS_H


struct integrand_params {
    double alpha;
    double t;
    int n;
    double h;
    int k;
    double* y;
};

struct dvbyt_integrand_params {
    double h;
    double alpha;
    double u;
    double t;
};

struct outer_density_integrand_params {
    double* y;
    int n;
    double t;
    double alpha;
    double h;
    int k;
};

struct inner_density_integrand_params {
    int n;
    double t;
    double alpha;
    double h;
    int k;
    double s;
};

struct sndsn_integrand_params {
    double* s;
    double* y;
    double t;
    int n;
    double h;
    double alpha;
    double eta;
    double tau;
    double v;
};

#endif //CPP_INTEGRAND_PARAMS_H
