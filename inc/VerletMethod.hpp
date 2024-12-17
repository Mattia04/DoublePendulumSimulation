#pragma once

#include <span>
#include <algorithm>

#include "EquazioneDifferenzialeBase.hpp"

class Verlet : public EquazioneDifferenzialeBase
{

public:
  virtual std::vector<double> Passo(double t,
                                    const std::vector<double> &x,
                                    double h,
                                    const FunzioneVettorialeBase &f) const override
  {
    size_t n = x.size() / 2; // Number of particles

    std::vector<double> pos(x.begin(), x.begin() + n); // Extract positions
    std::vector<double> vel(x.begin() + n, x.end());   // Extract velocities

    // Step 1: Compute acceleration using current position and velocity
    std::vector<double> acc(n);
    std::vector<double> temp = f(t, x);
    std::copy(temp.begin() + n, temp.end(), acc.begin());

    // Step 2: Half-step velocity update
    std::vector<double> vel_half(n);
    for (size_t i = 0; i < n; ++i)
      vel_half[i] = vel[i] + 0.5 * acc[i] * h;

    // Step 3: Position update
    std::vector<double> new_pos(n);
    for (size_t i = 0; i < n; ++i)
      new_pos[i] = pos[i] + vel_half[i] * h;

    // Step 4: Recompute acceleration with updated position and half-step velocity
    std::vector<double> temp_state(2 * n);
    std::copy(new_pos.begin(), new_pos.end(), temp_state.begin());
    std::copy(vel_half.begin(), vel_half.end(), temp_state.begin() + n);

    std::vector<double> new_acc = f(t + h, temp_state);

    // Step 5: Full-step velocity update
    std::vector<double> new_vel(n);
    for (size_t i = 0; i < n; ++i)
      new_vel[i] = vel_half[i] + 0.5 * new_acc[i] * h;

    // Combine new positions and velocities into the output vector
    std::vector<double> new_state(2 * n);
    std::copy(new_pos.begin(), new_pos.end(), new_state.begin());
    std::copy(new_vel.begin(), new_vel.end(), new_state.begin() + n);

    return new_state;
  };
};
