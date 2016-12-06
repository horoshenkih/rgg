//
// Created by serkh on 12/6/16.
//


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

    return gradients;
}