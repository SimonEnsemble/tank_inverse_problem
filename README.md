# Inferring the cross-sectional area profile of an unseen solid in a draining tank from liquid level dynamics

this repo contains data and code to reproduce the results in our article:
>Gbenga Fabusola, Cory M. Simon, Inferring the cross-sectional area profile of an unseen solid in a draining tank from liquid level dynamics, Chemical Engineering Science, Volume 309, 2025, 121488, ISSN 0009-2509, https://doi.org/10.1016/j.ces.2025.121488.


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

## drawings and photos

the folder `drawings_and_photos/` contains several sketches of the tank geometry and photos of the draining tank during the experiment.

## simple activity for the classroom

(not discussed in the paper, but uses data from this project.)

the Pluto notebook `classroom/model_tank.jl` contains a simple learning activity for the classroom: model the liquid level in the draining tank and fit the coefficient of discharge to the data.
