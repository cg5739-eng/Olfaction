# Olfaction

A 2D simulation of a robot locating an odor source using four directional air intakes and proportional control.

## Physics

The room is represented by a 5 × 5 m concentration field, not individual particles. Odor is continuously emitted by a stationary source and transported by advection and diffusion:

`∂C/∂t = −∇·(uC) + D∇²C + S − λC + aspiration exchange`

Airflow is the sum of prescribed room ventilation and aspiration-induced flow. Four equal intake sinks and a balanced central exhaust are modeled using a potential-flow Poisson equation. Aspirated odor is mixed and returned to the room.

The robot measures weighted concentrations in front of its four inlets. Its velocity is proportional to opposing sensor differences:

- vx = Kp (C+x − C−x)
- vy = Kp (C+y − C−y)

The controller does not know the source location.

## Assumptions and limitations

- Uniform temperature and background pressure; concentration is uniform along the vertical (z) axis.
- Odor is a passive scalar. Concentration units are arbitrary.
- Still, uniform, recirculating and time-varying synthetic ventilation are available; this is not full CFD.
- Sensors are instantaneous, identical and noiseless, with no saturation.
- The robot has a 30 cm visual body but does not obstruct airflow.
- The simulation stops within 20 cm of the source using ground-truth position, not an onboard detection algorithm.
- Gradient ascent can fail in weak or misleading concentration fields. Parameters are illustrative, not experimentally calibrated.

## Run locally

Install Python 3, download this repository and open a terminal in its folder. VS Code is not required.

```sh
python3 -m pip install -r requirements.txt
python3 odor_sim.py
```

On Windows, use `python` instead of `python3` if needed.

Left-click the room to move the source; right-click to move the robot. Adjust ventilation, gain, suction and simulation speed. Use pre-fill to disperse odor before moving the robot, and reset to start again.

Run numerical checks with `python3 odor_sim.py --test`.
