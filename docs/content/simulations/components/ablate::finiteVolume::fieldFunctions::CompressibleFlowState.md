---
layout: default
title: ablate::finiteVolume::fieldFunctions::CompressibleFlowState
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::finiteVolume::fieldFunctions::CompressibleFlowState](./ablate::finiteVolume::fieldFunctions::CompressibleFlowState.html)
## ablate::finiteVolume::fieldFunctions::CompressibleFlowState*
a simple structure used to describe a compressible flow field using an EOS, T, pressure, vel, Yi

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the eos used for the flow field

temperature (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the temperature field (K)

pressure (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the pressure field (Pa)

velocity (req) 
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the velocity field (m/2)

massFractions
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html)) a fieldFunctions used to describe all mass fractions

