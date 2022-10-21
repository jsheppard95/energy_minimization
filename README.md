Jackson Sheppard
CH E 210D, Exercise 2
10/20/22

Here we present an implementation of Conjugate Gradient Descent to find
minimum energy stable packings of a collection of spherical particles. We
consider the dimensionless Lennard-Jones interaction potential and initialize
random configurations including two to 26 stable particles, then perform
conjugate gradient descent to find minimum energy configurations. We repeat
this procedure `K` times for each system particle number `N`, with `K` set to
`100`, `1000`, and `10000`, plotting the minimum, average, and maximum result
over the `K` trials of the conjugate gradient minimization search. Finally,
we plot the minimum and average results found, and compare to both the
numerical results taken from
[Leary, J. Global Optimization 11, 35 (1997)](https://link.springer.com/article/10.1023/A:1026500301312)
and a predicted model for the global energy minimum as a function of `N` based
on surface area (via surface tension) and volume (via bulk energy density),
fitting the model using least-squares minimization. We then look for dips in
the residual plot for this model, indicating a lower stable energy minima than
expected from macroscopic scaling, thus indicating a stable packing
configuration.