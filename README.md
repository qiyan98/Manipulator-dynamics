# Manipulator-friction-dynamics-simulation
Project on the friction dynamics and control simulation for a manipulator, 2018-2020 at SJTU.

The paper *Tribo-Dynamic Analysis and Motion Control of a Rotating Manipulator Based on the Load and Temperature Dependent Friction Model*  is **accepted** by *Journal of Engineering Tribology*, coming online soon. 

## Abstract

Joint friction has a significant influence on the dynamics and motion control of manipulators. However, the friction effect is often omitted or simplified in previous studies. In this paper, we establish the dynamics model of a single joint in a rotating industrial manipulator taking detailed friction effects into consideration and propose a new control algorithm for the friction compensation purpose. Firstly, the manipulator dynamics modeling is carried out employing a recently-proposed extended static friction model, which depicts load and temperature influence on Coulomb, Stribeck and viscous terms. Moreover, based on the established dynamics model, the paper presents a new adaptive fast nonsingular terminal sliding mode (AFNTSM) controller. The proposed approach has the advantages of continuous control inputs, fast convergence rate, no singularity and great robustness against disturbances. Furthermore, its adaptive property does not require any prior knowledge of the upper bound of the uncertainties. Finally, the proposed controller is applied to the manipulator joint trajectory tracking problem with varying friction subject to load and temperature changes. The numerical simulation verifies the effectiveness of our proposed method and its advantages over other controllers.

## To get started

Run `SixDOF_runsim.m` to start simulation. Tune control parameters for each controller in `SixDOF_simulation.m`.

Please note that simulation with additive white noise is significantly slower, especially when `ode45` or similar high-precision solver is used. For a quick evaluation, it's recommended to use forward Euler discretization with simulated noise.