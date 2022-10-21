# energy_minimization

Jackson Sheppard\
CH E 210D, Exercise 2\
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

## Installation and Usage
Note the following instructions have been tested using Linux Ubuntu 20.04.5
only. Compiling Fortran may bring additional challenges on other operating
systems.

Clone [this repository](https://github.com/jsheppard95/energy_minimization)
and navigate to its root. Install dependencies defined in the
`environment.yml` file using `conda`:

```
$ conda env create --name envname --file=environment.yml
```

Compile Fortran code `ex2lib.f90` as a Python module as follows:

```
$ f2py -c -m ex2lib ex2lib.f90
```

This generates the
["Shared Object" library file](https://superuser.com/questions/71404/what-is-an-so-file)
`ex2lib.cpython-39-x86_64-linux-gnu.so`. With this file in the working
directory, one can then `import` functions defined in `ex2lib.f90` as Python
functions, noting that casing is lowered when used in Python. For example, the
Fortran function `CalcForces` is available in Python as

```
import ex2lib
ex2lib.calcforces(*args)
```

making the force and energy calculations far more efficient than if
implemented directly in Python due to the fact that Fortran is a compiled
language.

Generate energy minimization data with the following minimization script. Note
this process can take on the order of days. `K`-specific files (see below) are
however written sequentially, allowing for analysis in parallel.

```
$ python exercise2.py
```

This will write three files to the output directory `data/`:

```
K100_energy_min.txt
K1000_energy_min.txt
K10000_energy_min.txt
```

Each of the following format:

```
N: {}, Minimum P.E: {}, Average P.E: {}, Maximum P.E: {}
```

where `N` ranges fropm `2` to `26` inclusive, and the minimum, average, and maximum
are over the `K` conjugate gradient searches for that `N`.

Generate plots and global energy minimum model with the following analysis
script:

```
$ python analysis.py
```

This script reads the minimization data contained in `data` (along with the
minimizations obtained by Leary and defined in `data/LearyData.csv`), and
plots the results, writing figures toi the output directory `plots`.

## Results
