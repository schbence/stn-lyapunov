# Lyapunov measure for State-Transition Networks

The repository provides a method to construct State-Transition Networks (STN) from discrete-time or quasi-continuous multivariate time-series.
Furthermore it contains an implementation of the Lyapunov measure defined for STNs (as introduced in [article link?]). 

## Constructing State-Transition Networks
State-transition networks can be constructed using the function `STN` from `stn.py`.
Given the time-series data `data` and the discretization resolution `b` one obatins a weighted and directed `igraph.Graph` object.

### Example 1 - STN from Lorenz dynamics

<img src="./plots/ex-1-lorenz.png" width="300">

<img src="./plots/ex-1-stn.png" width="250">

## Calculating the Lyapunov network measure

### Example 2 - Lyapunov measure of STNs from the Henon Map

### Example 3 - Lyapunov measure of STNs from the Lorenz system

