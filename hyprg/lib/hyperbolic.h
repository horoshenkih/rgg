//
// Created by serkh on 12/6/16.
//

#ifndef HYPRG_HYPERBOLIC_H
#define HYPRG_HYPERBOLIC_H

#include <cmath>
#include <vector>
#include "graph.h"

double distance(double, double, double, double);
std::vector<double> grad_distance(double, double, double, double);

Coordinates poincare_average(vector<Coordinates>&);
Coordinates poincare_average(vector<Coordinates>&, vector<double>&);
#endif //HYPRG_HYPERBOLIC_H
