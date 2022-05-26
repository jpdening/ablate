---
layout: default
title: ablate::finiteVolume::processes::PressureGradientScaling
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteVolume::processes::PressureGradientScaling](./ablate::finiteVolume::processes::PressureGradientScaling.html)
## ablate::finiteVolume::processes::PressureGradientScaling*
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

