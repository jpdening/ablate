---
layout: default
title: ablate::solver::TimeStepper
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::solver::TimeStepper](./ablate::solver::TimeStepper.html)
## ablate::solver::TimeStepper*
the basic stepper

name (req) 
: (string) the time stepper name

domain (req) 
: ([ablate::domain::Domain](./ablate::domain::Domain.html)) the mesh used for the simulation

arguments (req) 
: (string,string map) arguments to be passed to petsc

io
: ([ablate::io::Serializer](./ablate::io::Serializer.html)) the serializer used with this timestepper

initialization
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) initialization field functions

exactSolution
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) optional exact solutions that can be used for error calculations

absoluteTolerances
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) optional absolute tolerances for a field

relativeTolerances
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) optional relative tolerances for a field


