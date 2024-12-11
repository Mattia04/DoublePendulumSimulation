#pragma once

#include "EquazioneDifferenzialeBase.hpp"

// integratore concreto, metodo di Eulero

class Eulero : public EquazioneDifferenzialeBase
{

public:
    virtual std::vector<double> Passo(double t,
                                      const std::vector<double> &x,
                                      double h,
                                      const FunzioneVettorialeBase &f) const override
    {
        return x + h * f(t, x);
    };
};
