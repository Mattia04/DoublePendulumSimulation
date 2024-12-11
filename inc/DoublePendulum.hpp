#pragma once

#include <cmath>

#include "FunzioneVettorialeBase.hpp"

class DoublePendulum : public FunzioneVettorialeBase {
public:
    DoublePendulum() = default;
    DoublePendulum(double g) : m_g(g) {};
    DoublePendulum(double g, double l) : m_g(g), m_l(l) {};

    virtual ~DoublePendulum() = default;

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

class DoublePendulumDamped : public FunzioneVettorialeBase {
    public:
        DoublePendulumDamped() = default;
        DoublePendulumDamped(double l1, double l2, double m1, double m2) : m_l1(l1), m_l2(l2), m_m1(m1), m_m2(m2) {
                // TODO add variable check
        };
        DoublePendulumDamped(double g, double l1, double l2, double m1, double m2, double gamma) : m_g(g), m_l1(l1), m_l2(l2), m_m1(m1), m_m2(m2), m_gamma1(gamma), m_gamma2(gamma) {
                // TODO add variable check
        };
        DoublePendulumDamped(double g, double l1, double l2, double m1, double m2, double gamma1, double gamma2) : m_g(g), m_l1(l1), m_l2(l2), m_m1(m1), m_m2(m2), m_gamma1(gamma1), m_gamma2(gamma2) {
                // TODO add variable check
        };

        virtual ~DoublePendulumDamped() = default;

		virtual std::vector<double> Eval(double t, const std::vector<double> &x) const{
            double m_cos = std::cos(x[0]-x[1]);
            double m_sin = std::sin(x[0]-x[1]);
            double m_sin1 = std::sin(x[0]);
            double m_sin2 = std::sin(x[1]);
            double m_m = m_m1 + m_m2;

            double ddtheta2 = -( - m_m2 * m_l2 * x[3]*x[3] * m_sin * m_cos / m_m - m_l1 * x[2] * x[2] * m_sin // parte cinetica
                               + m_g * (m_sin2 - m_cos * m_sin1) // parte gravitazionale
                               + m_gamma2 * x[3] / (m_l2 * m_m2) - m_gamma1 * x[2] * m_cos / (m_l1 * m_m) // parte viscosa
                               )/(m_l2 * (1 - m_m2 * m_cos * m_cos / m_m)); // coefficiente

            double ddtheta1 = -( m_m2 * m_l2 * m_cos * ddtheta2 // soluzione di theta 2
                                + m_m2 * m_l2 * x[3] * x[3] * m_sin // parte cinetica
                                + m_g * m_sin1 * m_m // parte gravitazionale
                                + m_gamma1 * x[2] / m_l1 // parte viscosa
            )/ (m_l1 * m_m); // coefficiente

            return std::vector<double> {x[2], x[3], ddtheta1, ddtheta2};
		}

		double GetL1() const { return m_l1; }
        double GetL2() const { return m_l2; }

        void SetG(double g) { m_g = g; }
        void SetM1(double m) { m_m1 = m; }
        void SetM2(double m) { m_m2 = m; }
        void SetL1(double l) { m_l1 = l; }
        void SetL2(double l) { m_l2 = l; }
        void SetGamma1(double gamma) { m_gamma1 = gamma; }
        void SetGamma2(double gamma) { m_gamma2 = gamma; }

    private:
        double m_g = 9.806;
        double m_l1 = 1.0;
        double m_l2 = 1.0;
        double m_m1 = 1.0;
        double m_m2 = 1.0;
        double m_gamma1 = 0.00;
        double m_gamma2 = 0.00;
};