#pragma once

#include "EquazioneDifferenzialeBase.hpp"

// integratore concreto, metodo di Runge-Kutta

class RungeKutta : public EquazioneDifferenzialeBase
{

public:
    virtual std::vector<double> Passo(double t,
                                      const std::vector<double> &x,
                                      double h,
                                      const FunzioneVettorialeBase &f) const override
    {
        std::vector<double> k1, k2, k3, k4;
        k1 = f(t, x);
        k2 = f(t + h / 2., x + k1 * h / 2.);
        k3 = f(t + h / 2., x + k2 * h / 2.);
        k4 = f(t + h, x + k3 * h);

        return x + (k1 + 2. * k2 + 2. * k3 + k4) * (h / 6.);
    };
};
