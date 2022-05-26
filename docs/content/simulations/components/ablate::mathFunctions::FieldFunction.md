---
layout: default
title: ablate::mathFunctions::FieldFunction
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html)
## ablate::finiteVolume::fieldFunctions::Euler
initializes the euler conserved field variables based upon a CompressibleFlowState

state (req) 
: ([ablate::finiteVolume::fieldFunctions::CompressibleFlowState](./ablate::finiteVolume::fieldFunctions::CompressibleFlowState.html)) The CompressibleFlowState used to initialize

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) A subset of the domain to apply the field function

## ablate::finiteVolume::fieldFunctions::MassFractions
initializes the yi field function variables based upon a the list of functions and eos

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) The eos with the list of species

values (req) 
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) The list of mass fraction functions

## ablate::finiteVolume::fieldFunctions::DensityMassFractions
initializes the densityYi conserved field variables based upon a CompressibleFlowState

state (req) 
: ([ablate::finiteVolume::fieldFunctions::CompressibleFlowState](./ablate::finiteVolume::fieldFunctions::CompressibleFlowState.html)) The CompressibleFlowState used to initialize

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) A subset of the domain to apply the field function

## ablate::finiteVolume::fieldFunctions::DensityExtraVariables
initializes the densityEV conserved field variables based upon a CompressibleFlowState and specified EV

state (req) 
: ([ablate::finiteVolume::fieldFunctions::CompressibleFlowState](./ablate::finiteVolume::fieldFunctions::CompressibleFlowState.html)) The CompressibleFlowState used to initialize

functions (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html) list) The EV values in order

## ablate::mathFunctions::FieldFunction*
a field description that can be used for initialization or exact solution 

fieldName (req) 
: (string) the field name

field (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the math function used to describe the field

timeDerivative
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the math function used to describe the field time derivative

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) A subset of the domain to apply the field function

