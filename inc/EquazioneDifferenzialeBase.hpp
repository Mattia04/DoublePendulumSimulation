#pragma once

#include "VectorOperations.hpp"
#include "FunzioneVettorialeBase.hpp"
#include <cmath>

// ===========================================================================
// classe astratta per un integratore di equazioni differenziali
// ===========================================================================

class EquazioneDifferenzialeBase
{

public:
    virtual std::vector<double> Passo(double t,
                                      const std::vector<double> &x,
                                      double h,
                                      const FunzioneVettorialeBase &f) const = 0;
};
