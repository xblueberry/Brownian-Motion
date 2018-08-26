# Brownian-Motion

This project is part of the Computational Physics class taken in my Bachelor degree. 

* The first program allows the user to generate a random variable with normal distribution of parameters mean and variance chosen by the user with the Box-Muller method. The output is a dataset with the probability density function for each x value; a histogram can be easily generated from that. 

* In the second part, the random movement of N independent 2-dimensional particles is simulated. It is considered that all particles are situated in the origin at time = 0. In each step the particle position is defined by:

xn(t+Δt)=xn(t)+Δx

yn(t+Δt)=yn(t)+Δy

where Δx and Δy are random gaussian numbers with mean = 0 and variance = Δt(2kBT/λ). 

* In the third part, we add a boundary condition to the second program. We add a parameter to choose the position of an orthogonal wall to the x-axis. This wall reflects the particles on the x-axis, such that when they touch the wall their Δx is inverted.

