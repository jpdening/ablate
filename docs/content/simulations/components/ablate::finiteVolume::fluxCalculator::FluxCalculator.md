---
layout: default
title: ablate::finiteVolume::fluxCalculator::FluxCalculator
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)
## ablate::finiteVolume::fluxCalculator::Ausm
AUSM Flux Spliting: "A New Flux Splitting Scheme" Liou and Steffen, pg 26, Eqn (6), 1993

## ablate::finiteVolume::fluxCalculator::AverageFlux
Takes the average of the left/right faces.  Only useful for debugging.

## ablate::finiteVolume::fluxCalculator::OffFlux
Turns of convective flux through the face.

## ablate::finiteVolume::fluxCalculator::AusmpUp
A sequel to AUSM, Part II: AUSM+-up for all speeds, Meng-Sing Liou, Pages 137-170, 2006

mInf
: (double) the reference mach number

pgs
: ([ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)) Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows

## ablate::finiteVolume::fluxCalculator::Rieman
Exact Rieman Solution

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) only valid for perfect gas

## ablate::finiteVolume::fluxCalculator::Riemann2Gas
Exact Riemann Solution for 2 Perfect Gasses

eosL (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) only valid for perfect gas

eosR (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) only valid for perfect gas

## ablate::finiteVolume::fluxCalculator::RiemannStiff
Exact Riemann Solution for 2 Stiffened Gasses

eosL (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) only valid for perfect or stiffened gas

eosR (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) only valid for perfect or stiffened gas

