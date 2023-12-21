# REMUS - Rocket Equation Model Using Simulations

## Overview

REMUS is a C++ simulation program designed to numerically model the motion of a rocket using the Rocket Equation. The program solves the system from time 0 until the rocket returns to h = 0. It incorporates environmental and rocket-specific parameters provided in a file called `parameters.txt` and integrates the system using the 4th-order Runge-Kutta scheme.

## Project Background

This project initially aimed to model a concept rocket for the Altitude Record Team's APEX I. APEX I successfully broke the UKRA I-Class record, achieving an altitude of 9333 ft (2845 m) and reaching Mach 1.55.

## Features

- Reads environmental and rocket-specific parameters from `parameters.txt`.
- Integrates the system using the 4th-order Runge-Kutta scheme.
- Writes time `t` and values of `h`, `v`, and `m` at each time step to the file `output.txt`.
- Includes user-defined classes for rocket and rocket stage representation.
- Extensible to model multi-stage rockets.

## Usage

1. Ensure the `parameters.txt` file is in the current working directory with the specified format.
2. Run the REMUS program.
3. Check the `output.txt` file for the simulation results.

## Example `parameters.txt`

```plaintext
# Environmental parameters
1.22 9.81 0.2 0.0 0.0
# Rocket-specific parameters
2.0 0.0 0.001 0.005 100 0.01
```

## Contributing

If you are interested in contributing to REMUS, feel free to fork the repository and submit pull requests.

## License

This project is licensed under the MIT License.

## Acknowledgments

Special thanks to the Altitude Record Team for the real-world application of the Rocket Equation Model.

