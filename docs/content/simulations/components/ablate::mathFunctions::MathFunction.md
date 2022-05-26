---
layout: default
title: ablate::mathFunctions::MathFunction
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)
## ablate::mathFunctions::SimpleFormula*
a string based function to be parsed with muparser. The (string) formula that may accept x, y, z, t as variables

## ablate::mathFunctions::ConstantValue
sets a constant value to all values in field

## ablate::mathFunctions::ParsedSeries
 computes a series result from a string function with variables x, y, z, t, and i where i index of summation. $$\sum_{i = m}^n formula(x, y, z, t, n)$$

formula (req) 
: (string) see ParsedFunction for details on the string formatting.

lowerBound (req) 
: (int) the inclusive lower bound of summation (m)

upperBound (req) 
: (int) the inclusive upper bound of summation (n)

## ablate::mathFunctions::LinearTable
A table that is built from a spreadsheet that allows linear interpolation of variables based on monotonically increasing independent variables

file (req) 
: (file path or url) a file with csv data and header

independent (req) 
: (string) the name of the independent column name as defined in the header

dependent (req) 
: (string list) the names of the dependent column in the order in which to apply them

mappingFunction (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html))  the function that maps from the physical x,y,z, and t space to the table independent variable

## ablate::mathFunctions::Formula
 computes string function with variables x, y, z, and t where additional variables can be specified using other functions

formula (req) 
: (string) see ParsedFunction for details on the string formatting.

nested
: (string,[ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html) map) a map of nested MathFunctions.  These functions are assumed to compute a single scalar value

constants
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) constants that can be used in the formula

## ablate::mathFunctions::geom::Sphere
assigns a uniform value to all points inside the sphere

center (req) 
: (double list) the sphere center

radius
: (double) the sphere radius

insideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the values for inside the sphere, defaults to 1

outsideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the outside values, defaults to zero

## ablate::mathFunctions::geom::Box
assigns a uniform value to all points inside the box

lower (req) 
: (double list) the box lower corner

upper (req) 
: (double list) the box upper corner

insideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the values for inside the sphere, defaults to 1

outsideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the outside values, defaults to zero

## ablate::mathFunctions::geom::Surface
Assigned a unified number to all points inside of cad geometry file.

path (req) 
: (file path or url) the path to the step/stp file

insideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the values for inside the sphere, defaults to 1

outsideValues
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the outside values, defaults to zero

egadsVerboseLevel
: (int) the egads verbose level for output (default is 0, max is 3)

