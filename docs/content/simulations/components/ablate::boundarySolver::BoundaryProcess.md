---
layout: default
title: ablate::boundarySolver::BoundaryProcess
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::boundarySolver::BoundaryProcess](./ablate::boundarySolver::BoundaryProcess.html)
## ablate::boundarySolver::lodi::IsothermalWall
Enforces a isothermal wall with fixed velocity/temperature

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) The EOS describing the flow field at the wall

pgs
: ([ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)) Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows

## ablate::boundarySolver::lodi::OpenBoundary
Treats boundary as open.

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) The EOS describing the flow field at the boundary

reflectFactor (req) 
: (double) boundary reflection factor

referencePressure (req) 
: (double) reference pressure

maxAcousticsLength (req) 
: (double) maximum length in the domain for acoustics to propagate 

pgs
: ([ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)) Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows

## ablate::boundarySolver::lodi::Inlet
Enforces an inlet with specified velocity

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) The EOS describing the flow field at the wall

pgs
: ([ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)) Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows

velocity
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) optional velocity function that can change over time

## ablate::boundarySolver::physics::Sublimation
Adds in the euler/yi sources for a sublimating material.  Should be used with a LODI boundary.

latentHeatOfFusion (req) 
: (double) the latent heat of fusion [J/kg]

transportModel (req) 
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the effective conductivity model to compute heat flux to the surface [W/(m⋅K)]

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the eos used to compute temperature on the boundary

massFractions
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html)) the species to deposit the off gas mass to (required if solving species)

additionalHeatFlux
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) additional normal heat flux into the solid function

pgs
: ([ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)) Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows

disablePressure
: (bool) disables the pressure contribution to the momentum equation. Should be true when advection is not solved. (Default is false)

