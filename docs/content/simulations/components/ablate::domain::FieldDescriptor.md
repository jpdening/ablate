---
layout: default
title: ablate::domain::FieldDescriptor
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::domain::FieldDescriptor](./ablate::domain::FieldDescriptor.html)
## ablate::finiteVolume::CompressibleFlowFields
FVM fields need for compressible flow

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state to be used for the flow

extraVariables
: (string list) Any extra variables to transport

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region for the compressible flow (defaults to entire domain)

conservedFieldOptions
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) petsc options used for the conserved fields.  Common options would be petscfv_type and petsclimiter_type

## ablate::finiteElement::LowMachFlowFields
FE fields need for incompressible/low-Mach flow

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region for the compressible flow (defaults to entire domain)

includeSourceTerms
: (bool) include aux fields for source terms (defaults to false)

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

