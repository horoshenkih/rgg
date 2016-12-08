//
// Created by serkh on 12/6/16.
//
#include <exception>

#include "hyperbolic.h"

double cosh_distance(double r1, double phi1, double r2, double phi2) {
    return cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(phi1 - phi2);
}

double distance(double r1, double phi1, double r2, double phi2) {
    return acosh(cosh_distance(r1, phi1, r2, phi2));
}

std::vector<double> grad_distance(double r1, double phi1, double r2, double phi2) {
    double cd = cosh_distance(r1, phi1, r2, phi2);
    std::vector<double> gradients{0.,0.,0.,0.};
    if (fabs(cd - 1) < 1e-15) {
        return gradients;
    }
    double normalizer = 1. / sqrt(cd-1) / sqrt(cd+1);
    gradients[0] = sinh(r1) * cosh(r2) - cosh(r1) * sinh(r2) * cos(phi1 - phi2);
    gradients[1] = sinh(r1) * sinh(r2) * sin(phi1 - phi2);
    gradients[2] = cosh(r1) * sinh(r2) - sinh(r1) * cosh(r2) * cos(phi1 - phi2);
    gradients[3] = -gradients[1];
    for (int i = 0; i < gradients.size(); ++i) {
        gradients[i] *= normalizer;
    }

    return gradients;
}

Coordinates poincare_average(vector<Coordinates>& coords) {
    double A = 0;
    double B = 0;
    double rnorm = 0;
    for (Coordinates c : coords) {
        double r_i = c[0];
        double phi_i = c[1];
        A += sinh(r_i) * sin(phi_i);
        B += sinh(r_i) * cos(phi_i);
        rnorm += cosh(r_i);
    }
    if (fabs(A) < 1e-10 and fabs(B) < 1e-10) {
        return Coordinates{0,0};
    }
    double phi = atan2(A, B);
    double r = atanh((A * sin(phi) + B * cos(phi)) / rnorm);
    return Coordinates{r, phi};
}

Coordinates poincare_average(vector<Coordinates>& coords, vector<double>& weights) {
    if (coords.size() != weights.size()) {
        throw std::length_error("number of points and weights does not coinside");
    }
    // TODO copypaste
    double A = 0;
    double B = 0;
    double rnorm = 0;
    for (int i = 0; i < coords.size(); ++i) {
        Coordinates c = coords[i];
        double w_i = weights[i];
        double r_i = c[0];
        double phi_i = c[1];
        A += sinh(r_i) * sin(phi_i) * w_i;
        B += sinh(r_i) * cos(phi_i) * w_i;
        rnorm += cosh(r_i) * w_i;
    }
    if (fabs(A) < 1e-10 and fabs(B) < 1e-10) {
        return Coordinates{0,0};
    }
    double phi = atan2(A, B);
    double r = atanh((A * sin(phi) + B * cos(phi)) / rnorm);
    return Coordinates{r, phi};
}