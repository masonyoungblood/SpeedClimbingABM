Speed Climbing ABM
================

[![Static Badge](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.16905877-blue)](https://doi.org/10.5281/zenodo.16905877)
[![Static Badge](https://img.shields.io/badge/Paper-10.1098%2Frspb.2025.1433-red)](https://doi.org/10.1098/rspb.2025.1433)

![speed_climbing](https://helios-i.mashable.com/imagery/articles/03XKxYMwkNiZJGiHUu2FWnw/hero-image.fill.size_1248x702.v1628089817.png)  
&copy; Mashable

Data and code associated with ["Simulation-based inference with deep learning suggests speed climbers combine innovation and copying to improve performance"](https://doi.org/10.1098/rspb.2025.1433)&mdash;published in Proceedings of the Royal Society B.

This agent-based model simulates a dynamic population of professional speed climbers, and incorporates parameters for athletic improvement, innovation of “beta” (or route sequence), and copying of other climbers’ beta. Here are the locations of various files in the repository:

- `climbing_times/all_climbing_times.RData`: all climbing times for all competitors in each year.
- `climbing_times/best_climbing_times.RData`: best climbing time for each competitor in each year (used for simulation).
- `climbing_times/ifsc_world_championship` and `climbing_times/ifsc_world_championship.numbers`: manually-transcribed climbing sequences from IFSC events.
- `climbing_times/lfigil_ifsc_data` and `lfigil_processing.R`: climbing times from IFSC events (collected by [lfigil](https://github.com/lfigil/ifsc-data)) and a processing script to convert them to the RData files.
- `analysis`: full simulation-based inference pipeline.
- `SpeedClimbingABM.R`: source code for the agent-based model.
- `grid.csv` and `speed_wall.png`: coordinates and visualization of each hold on the route.

Files over 100 MB are ignored to avoid GitHub limit problems.
