`landsim` : Efficient simulation of landscape-scale populations

This is an R package with tools to run efficient, forwards-time simulation
of populations on habitats derived from GIS raster layers.

Install it by installing [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) and running
```
devtools::install_github("petrelharp/landsim")
```

Here are some vignettes showing what it can do:

1.  [Object structure](http://petrelharp.github.io/landsim/vignettes/object-structure.html) : Describes the modeling strategy and how information is stored in *population* and *demography* objects.
2.  [Density dependence](http://petrelharp.github.io/landsim/vignettes/density-dependence.html) : Demonstrates a few different ways of implementing density-dependent population size regulation.
3.  [Types of selection](http://petrelharp.github.io/landsim/vignettes/types-of-selection.html) : Demonstrates a few different ways of implementing genotype-dependent selection.
4.  [Spatial sweeps](http://petrelharp.github.io/landsim/vignettes/spatial-sweeps.html) : Simulates and analyzes an expanding wave of an advantageous allele.
5.  [A selective sweep](http://petrelharp.github.io/landsim/vignettes/selective-sweeps.html) : Works through the modeling choices of a model of selection.
6.  [Steps in a generation](http://petrelharp.github.io/landsim/vignettes/simulation-methods.html) : Works through the guts of the steps in a simulation.
7.  [Migration methods](http://petrelharp.github.io/landsim/vignettes/migration-methods.html) : Explains how to implement various types of migration: conservative or not, and different boundary types.
8.  [Migration method run times](http://petrelharp.github.io/landsim/vignettes/run-times.html) : Examines the details of two methods of migration, using precomputed sparse Matrix objects, and using `raster` methods.
9.  [Sierpinski habitat](http://petrelharp.github.io/landsim/vignettes/sierpinski.html) : Creates a random (Sierpinski gasket) habitat.
