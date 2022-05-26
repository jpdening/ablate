---
layout: default
title: ablate::eos::transport::TransportModel
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)
## ablate::eos::transport::Constant*
constant value transport model (often used for testing)

k
: (double) thermal conductivity [W/(m K)]

mu
: (double) viscosity [Pa s]

diff
: (double) diffusivity [m2/s]

## ablate::eos::transport::Sutherland
Sutherland Transport model

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) The EOS used to compute Cp (needed for Conductivity)

enabledProperties
: (ablate::eos::transport::TransportProperty enum list) list of enabled properties.  When empty or default all properties are enabled.

