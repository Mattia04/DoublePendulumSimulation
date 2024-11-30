#pragma once

#include <cmath>

#include "FunzioneVettorialeBase.hpp"

class DoublePendulum : public FunzioneVettorialeBase {
public:
    DoublePendulum() = default;
    DoublePendulum(double g) : m_g(g) {};
    DoublePendulum(double g, double l) : m_g(g), m_l(l) {};

    ~DoublePendulum() = default;

    virtual std::vector<double> Eval(double t, const std::vector<double> &x) const{
        double m_cos = std::cos(x[1]-x[0]);
        double m_sin = std::sin(x[1]-x[0]);

        double ddPhi = - (x[3]*x[3] * m_sin * m_cos / 2 + x[2] * x[2] * m_sin + m_g/m_l * (std::sin(x[1]) - std::sin(x[0]) * m_cos)) / (1 - m_cos*m_cos /2);

        return std::vector<double> {x[2],
          x[3],
            - ddPhi * m_cos/2 + x[3]*x[3] * m_sin / 2 - m_g/m_l * sin(x[0]),
            ddPhi
      };
    };
private:
    double m_g = 9.806; // m/s^2
    double m_l = 1.0; // m
};
