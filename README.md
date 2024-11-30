# Double Pendulum Simulation

This repo contains a simulator of a double pendulum.

The simulation is written in C++, I'm using Runge-Kutta's method of integration of order 4 to calculate the steps of the pendulum.

The current simulation is able to plot around 1000 pendulums without lagging (on my machine) in all configurations.

## Configurations:

When launching the program there are multiple argument you can give:

`./main <number of pendulums, default: 1> <show trace (1 true/0 false, default: 1)> <show pendulums (options: [0: don't draw, 1: draw just the point, 2: draw all], default: 2)> <length of the trace (default: 100)> <NOT IMPLEMENTED YET configuration space or phase space (options: [0: configuration space, 1: phase space], default: 0)> <NOT IMPLEMENTED YET damping>`

I'm considering of adding a CONFIG file since there are too many arguments.

## Pendulum Description

For now the double pendulum has the same length, the same masses and no damping.

In future I will add the possibility of having different length and masses and damping.

The simulation only shows configuration space, I will add phase space in future.

