# Double Pendulum Simulation

This repo contains a simulator of a double pendulum.

The simulation is written in C++, I have used Lagrangian mechanics to calculate the differential equations for the classes, and I'm using Runge-Kutta of order 4 as integrator (it is also available Euler method for integration, and I'm working on fixing Verlet's method).

The current simulation is able to plot around 1000 pendulums without lagging (on my machine) in all possible configurations.

## Examples

Here are two videos simulating 500 pendulums, the pendulums witout trace and just the trace.

![video2](https://github.com/user-attachments/assets/ba59a9b0-ce9a-4156-b8d3-abbb8c160d2e)

![video1](https://github.com/user-attachments/assets/6e861552-fef6-4da5-b3b9-d9030a9eb070)

Here is another video demostrating that the period of the pendulum oscillations (for small oscillations) is in function to the length of the pendulum.

![video4](https://github.com/user-attachments/assets/e4abb7db-5750-466c-86bf-9dae4a4769a5)

## Configurations:

Since there are a lot of parameters you could want to change i decided to create a `SETTINGS.yaml` file to store your favourite configuration.

I won't list here all possible settings since they are too many, but I will list the most important ones here:

- With `show_pendulums`, `show_trace` you can show the pendulums and the trace left by the mass2, you can also change the trace length with `trace_length`.
- You can simulate a different `number_of_pendulums`, or change the `time_step` (a time step of 0.01 is recommended for "rk4"), or the `total_time`.
- Most importantly you can change the pendulums parameters in the `pendulums` section, you can also modify a single pendulum by creating a `pendulum_<number>` section with the modified parameters, you can modify every parameter already contained in `pendulums` except for `spread`.

## Pendulum Description

### `DoublePendulum` class

The `DoublePendulum` class creates a double pendulum with: the same rod lengths, the same masses and no damping.

### `DoublePendulumDamped` class

The `DoublePendulumDamped` class creates a double pendulum with: different rod lengths, different masses and damping (for each mass).

The class constructors do not accept non-physical values for the pendulum, but if you want to use a negative mass or inverted gravity. You can still achieve (at your risk) this by using the `Set<parameter>` methods (Note: NEVER use masses or lengths equal to zero). 

Note: there is an issue with the damping if using large coefficients (~ 1), but it works fine using smaller ones (< 0.5).

## TODOs

- Add phase space visualization
- Add marker size that goes like: `sqrt(mass)`
- Fix damping (for `damping1 << 1` and `damping2 > 1`).
- Verlet method of integration is not working as intended.
- Add energy logging to an external file to be displayed with matplotlib at the end of simulation.
- Add multithreading. Edit: I'm thinking of doing this in a different project combined with GPU computation. 
  - The project idea is to simulate a huge number of double pendulums (I'm thinking the order of a million) at the same time to make a video about the pseudo-equilibrium zones (with colors).
