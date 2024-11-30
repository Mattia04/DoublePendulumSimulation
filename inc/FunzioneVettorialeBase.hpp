#pragma once

#include "VectorOperations.hpp"

// ===========================================================================
// classe astratta, restituisce la derivata valutata nel punto x
// ===========================================================================

class FunzioneVettorialeBase
{

public:
    virtual std::vector<double> Eval(double t, const std::vector<double> &x) const = 0;

    std::vector<double> operator()(double t, const std::vector<double> &x) const { return Eval(t, x); };
};
