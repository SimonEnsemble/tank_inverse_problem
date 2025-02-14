# reconstructing the shape of a solid in a draining tank from its liquid level dynamics

this repo contains data and code to reproduce the results in our article:
> G. Fabusola, C. Simon. "Reconstructing the shape of a solid inside a draining tank from its liquid level dynamics." _Chemical Engineering Science. (under revision)_ (2024)

## raw data

the `data` folder contains:
* `tank_length_measurements.txt`: raw length measurements to characterize the tank geometry
* `level_sensor_calibration.csv`: calibration curve for the liquid level sensor
* `liq_level_data_empty_*.csv`: liquid level data for draining tank _without_ an exogenous object inside of it
* `liq_level_data_w_bottle.csv`: liquid level data for draining tank _with_ a glass bottle inside of it
* `bottle_area.csv`: length measurements of the circumference of the bottle (--> inferred cross-sectional area)

## code

* `tank_problem.jl`, a Pluto notebook, contains all Julia code to reproduce the results and draw the figures in our article. this notebook reads in the data in `data` and conducts Bayesian inference. the notebook is well-organized and documented.
* `Water_level_measurement_Audrino.ino` contains the Arduino code to read data from our liquid level sensor.

## drawings and photo

the folder `drawings_and_photo/` contains several sketches of the tank geometry and a photo of the tank during the experiment.
