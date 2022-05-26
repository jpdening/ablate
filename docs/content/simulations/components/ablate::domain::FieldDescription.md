---
layout: default
title: ablate::domain::FieldDescription
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::domain::FieldDescription](./ablate::domain::FieldDescription.html)
## ablate::domain::FieldDescription*
A single custom field description

name (req) 
: (string) the name of the field

prefix
: (string) optional prefix (defaults to name)

components
: (string list) the components in the field (defaults to 1)

location
: (ablate::domain::FieldLocation enum) if it is a solution (SOL) or auxiliary (aux) field

type (req) 
: (ablate::domain::FieldType enum) if it is a finite volume (FV) or finite element (FE) field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region in which this field lives

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) field specific options

