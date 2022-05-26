---
layout: default
title: ablate::eos::EOS
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::eos::EOS](./ablate::eos::EOS.html)
## ablate::eos::PerfectGas
perfect gas eos

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) parameters for the perfect gas eos

species
: (string list) species to track.  Note: species mass fractions do not change eos

## ablate::eos::TChem
TChem ideal gas eos

mechFile (req) 
: (file path or url) the mech file (CHEMKIN Format)

thermoFile (req) 
: (file path or url) the thermo file (CHEMKIN Format)

## ablate::eos::StiffenedGas
stiffened gas eos

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) parameters for the stiffened gas eos

species
: (string list) species to track.  Note: species mass fractions do not change eos

