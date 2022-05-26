---
layout: default
title: ablate::finiteVolume::boundaryConditions::BoundaryCondition
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteVolume::boundaryConditions::BoundaryCondition](./ablate::finiteVolume::boundaryConditions::BoundaryCondition.html)
## ablate::finiteVolume::boundaryConditions::EssentialGhost
essential (Dirichlet condition) for ghost cell based boundaries

boundaryName (req) 
: (string) the name for this boundary condition

labelIds (req) 
: (int list) the ids on the mesh to apply the boundary condition

boundaryValue (req) 
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html)) the field function used to describe the boundary

labelName
: (string) the mesh label holding the boundary ids (default Face Sets)

enforceAtFace
: (bool) optionally update the boundary to enforce the value at the face instead of the cell using linear interpolation (default falue)

