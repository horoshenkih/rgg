//
// Created by serkh on 12/6/16.
//

#include "loss_function.h"

double SmoothEdgePredictor::predict(double d, double r, double beta) const {
    return 1. / (1. + exp(beta * (d - r)));
}

double SmoothEdgePredictor::prediction_gradient_d(double d, double r, double beta) const {
    double p = predict(d, r, beta);
    return -exp((d - r) * beta) * beta * pow(p,2);
}

double SmoothEdgePredictor::prediction_gradient_r(double d, double r, double beta) const {
    return -prediction_gradient_d(d, r, beta);
}

SmoothEdgeLoss::SmoothEdgeLoss() {
    beta = 1.;
    EPS = 1e-15;
    edge_predictor = SmoothEdgePredictor();
}

double SmoothEdgeLoss::loss(const Edge& edge, const DistanceInfo& distance_info) const {
    double radius = distance_info.get_radius();
    double distance = distance_info.get_distance();
    double t = edge.get_value();
    double p = edge_predictor.predict(distance, radius, beta);
    double edge_weight = edge.get_weight();

    return edge_weight * elementary_loss(t, p);
}

Gradient SmoothEdgeLoss::loss_gradient(const Edge& edge, const DistanceInfo& distance_info) const {
    double radius = distance_info.get_radius();
    double distance = distance_info.get_distance();
    double t = edge.get_value();
    double p = edge_predictor.predict(distance, radius, beta);
    double edge_weight = edge.get_weight();

    double dp = edge_predictor.prediction_gradient_d(distance, radius, beta);
    DistanceGradients dg = distance_info.get_distance_gradients();
    vector<double> dd{0,0,0,0};
    pair<Node, Node> nodepair = edge.get_node_pair();
    Node v1 = nodepair.first;
    Node v2 = nodepair.second;
    dd[0] = dg[v1][0];
    dd[1] = dg[v1][1];
    dd[2] = dg[v2][0];
    dd[3] = dg[v2][1];
    vector<double> d_grad;
    double multiplier = dp * edge_weight * elementary_loss_gradient(t, p);
    for (auto x : dd) {
        d_grad.push_back(x * multiplier);
    }

    return d_grad;
}

double MSE::elementary_loss(double t, double p) const { return pow(p-t, 2); }
double MSE::elementary_loss_gradient(double t, double p) const { return 2. * (p-t); }

double LogLoss::elementary_loss(double t, double p) const {
    p = min(max(EPS, p), 1 - EPS);
    return -t * log(p) - (1 - t) * log(1 - p);
}
double LogLoss::elementary_loss_gradient(double t, double p) const {
    p = min(max(EPS, p), 1 - EPS);
    return (p - t) / p / (1 - p);
}
