---
layout: default
title: ablate::particles::initializers::Initializer
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::particles::initializers::Initializer](./ablate::particles::initializers::Initializer.html)
## ablate::particles::initializers::CellInitializer
simple cell initializer that puts particles in every element

particlesPerCellPerDim (req) 
: (int) particles per cell per dimension

## ablate::particles::initializers::BoxInitializer
simple box initializer that puts particles in a defined box

lower (req) 
: (double list) the lower bound of the box

upper (req) 
: (double list) the upper bound of the box

particlesPerDim (req) 
: (int) the particles per box dimension

