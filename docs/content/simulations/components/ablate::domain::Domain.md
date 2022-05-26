---
layout: default
title: ablate::domain::Domain
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::domain::Domain](./ablate::domain::Domain.html)
## ablate::domain::BoxMesh
simple uniform box mesh

name (req) 
: (string) the name of the domain/mesh object

fields
: ([ablate::domain::FieldDescriptor](./ablate::domain::FieldDescriptor.html) list) a list of fields/field descriptors

modifiers
: ([ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html) list) a list of domain modifier

faces (req) 
: (int list) the number of faces in each direction

lower (req) 
: (double list) the lower bound of the mesh

upper (req) 
: (double list) the upper bound of the mesh

boundary
: (string list) custom boundary types (NONE, GHOSTED, MIRROR, PERIODIC)

simplex
: (bool) sets if the elements/cells are simplex

## ablate::domain::DMPlex
DMPlex that can be set using PETSc options

fields
: ([ablate::domain::FieldDescriptor](./ablate::domain::FieldDescriptor.html) list) a list of fields/field descriptors

name (req) 
: (string) the mesh dm name

modifiers
: ([ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html) list) a list of domain modifier

## ablate::domain::FileMesh
read a DMPlex from a file

name (req) 
: (string) the name of the domain/mesh object

path (req) 
: (file path or url) the path to the mesh file

fields
: ([ablate::domain::FieldDescriptor](./ablate::domain::FieldDescriptor.html) list) a list of fields/field descriptors

modifiers
: ([ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html) list) a list of domain modifier

## ablate::domain::BoxMeshBoundaryCells
simple uniform box mesh with boundary solver cells.  Available labels are: interiorCells, domain (interior and boundary cells), boundaryFaces, boundaryCells, boundaryCellsLeft, boundaryCellsRight, boundaryCellsBottom, boundaryCellsTop, boundaryCellsFront, and boundaryCellsBack

name (req) 
: (string) the name of the domain/mesh object

fields
: ([ablate::domain::FieldDescriptor](./ablate::domain::FieldDescriptor.html) list) a list of fields/field descriptors

preModifiers
: ([ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html) list) a list of domain modifiers to apply before ghost labeling

postModifiers
: ([ablate::domain::modifiers::Modifier](./ablate::domain::modifiers::Modifier.html) list) a list of domain modifiers to apply after ghost labeling

faces (req) 
: (int list) the number of faces in each direction

lower (req) 
: (double list) the lower bound of the mesh

upper (req) 
: (double list) the upper bound of the mesh

simplex
: (bool) sets if the elements/cells are simplex

