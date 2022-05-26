---
layout: default
title: ablate::finiteVolume::processes::Process
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteVolume::processes::Process](./ablate::finiteVolume::processes::Process.html)
## ablate::finiteVolume::processes::EulerTransport
build advection/diffusion for the euler field

parameters
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the parameters used by advection

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used to describe the flow

fluxCalculator
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) the flux calculator (default is no advection)

transport
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the diffusion transport model (default is no diffusion)

## ablate::finiteVolume::processes::SpeciesTransport
diffusion/advection for the species yi field

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used to describe the flow

fluxCalculator
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) the flux calculator (default is no advection)

transport
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the diffusion transport model (default is no diffusion)

## ablate::finiteVolume::processes::EVTransport
diffusion/advection for the specified EV

conserved (req) 
: (string) the name of the conserved (density*ev) of the variable

nonConserved (req) 
: (string) the name of the non-conserved (ev) of the variable

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used to describe the flow

fluxCalculator
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) the flux calculator (default is no advection)

transport
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the diffusion transport model (default is no diffusion)

## ablate::finiteVolume::processes::TChemReactions
reactions using the TChem v1 library

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the tChem v1 eos

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) any PETSc options for the chemistry ts

inertSpecies
: (string list) fix the Jacobian for any undetermined inertSpecies

massFractionBounds
: (double list) sets the minimum/maximum mass fraction passed to TChem Library. Must be a vector of size two [min,max] (default is no bounds)

## ablate::finiteVolume::processes::TwoPhaseEulerAdvection


eosGas (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) 

eosLiquid (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) 

fluxCalculatorGasGas (req) 
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) 

fluxCalculatorGasLiquid (req) 
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) 

fluxCalculatorLiquidGas (req) 
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) 

fluxCalculatorLiquidLiquid (req) 
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) 

## ablate::finiteVolume::processes::Gravity
build advection/diffusion for the euler field

vector (req) 
: (double list) gravitational acceleration vector

## [ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)
Rescales the thermodynamic pressure gradient scaling the acoustic propagation speeds to allow for a larger time step.

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used for the flow

alphaInit (req) 
: (double) the initial alpha

domainLength (req) 
: (double) the reference length of the domain

maxAlphaAllowed
: (double) the maximum allowed alpha during the simulation (default 100)

maxDeltaPressureFac
: (double) max variation from mean pressure (default 0.05)

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

## ablate::finiteVolume::processes::ArbitrarySource
uses math functions to add arbitrary sources to the fvm method

## ablate::finiteVolume::processes::Buoyancy
build advection/diffusion for the euler field

vector (req) 
: (double list) gravitational acceleration vector

## ablate::finiteVolume::processes::ConstantPressureFix
Reset the density based upon a constant pressure value

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used to describe the flow

pressure (req) 
: (double) the constant pressure value

