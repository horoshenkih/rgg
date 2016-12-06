//
// Created by serkh on 12/6/16.
//

#ifndef HYPRG_LOSS_FUNCTION_H
#define HYPRG_LOSS_FUNCTION_H

#include <cmath>
#include <algorithm>
#include "embedding_model.h"

using std::max;
using std::min;

class LossFunction {
public:
    virtual double loss(const Edge&, const DistanceInfo&) const = 0;
    virtual Gradient loss_gradient(const Edge&, const DistanceInfo&) const = 0;
};

class SmoothEdgePredictor {
public:
    double predict(double d, double r, double beta) const;
    double prediction_gradient_d(double d, double r, double beta) const;
    double prediction_gradient_r(double d, double r, double beta) const;
};

class SmoothEdgeLoss : public LossFunction {
protected:
    double beta;
    SmoothEdgePredictor edge_predictor;
    double EPS;
public:
    SmoothEdgeLoss();
    virtual double elementary_loss(double t, double p) const = 0;
    virtual double elementary_loss_gradient(double t, double p) const = 0;
    double loss(const Edge&, const DistanceInfo&) const;
    Gradient loss_gradient(const Edge&, const DistanceInfo&) const;
};

class MSE : public SmoothEdgeLoss {
public:
    double elementary_loss(double t, double p) const;
    double elementary_loss_gradient(double t, double p) const;
};

class LogLoss : public SmoothEdgeLoss {
public:
    double elementary_loss(double t, double p) const;
    double elementary_loss_gradient(double t, double p) const;
};

#endif //HYPRG_LOSS_FUNCTION_H
