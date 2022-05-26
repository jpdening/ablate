---
layout: default
title: ablate::finiteElement::boundaryConditions::BoundaryCondition
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteElement::boundaryConditions::BoundaryCondition](./ablate::finiteElement::boundaryConditions::BoundaryCondition.html)
## ablate::finiteElement::boundaryConditions::Essential
essential (Dirichlet condition) for FE based problems

boundaryName (req) 
: (string) the name for this boundary condition

labelIds (req) 
: (int list) the ids on the mesh to apply the boundary condition

boundaryValue (req) 
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html)) the field function used to describe the boundary

labelName
: (string) the mesh label holding the boundary ids (default marker)

