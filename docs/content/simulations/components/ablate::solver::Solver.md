---
layout: default
title: ablate::solver::Solver
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::solver::Solver](./ablate::solver::Solver.html)
## ablate::finiteVolume::FiniteVolumeSolver
finite volume solver

id (req) 
: (string) the name of the flow field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the options passed to PETSC for the flow

processes (req) 
: ([ablate::finiteVolume::processes::Process](./ablate::finiteVolume::processes::Process.html) list) the processes used to describe the flow

boundaryConditions
: ([ablate::finiteVolume::boundaryConditions::BoundaryCondition](./ablate::finiteVolume::boundaryConditions::BoundaryCondition.html) list) the boundary conditions for the flow field

computePhysicsTimeStep
: (bool) determines if a physics based time step is used to control the FVM time stepping (default is false)

## ablate::finiteVolume::CompressibleFlowSolver
compressible finite volume flow

id (req) 
: (string) the name of the flow field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the options passed to PETSc

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the equation of state used to describe the flow

parameters
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the parameters used for field values

transport
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the diffusion transport model

fluxCalculator
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) the flux calculators (defaults to none)

additionalProcesses
: ([ablate::finiteVolume::processes::Process](./ablate::finiteVolume::processes::Process.html) list) any additional processes besides euler/yi/ev transport

boundaryConditions
: ([ablate::finiteVolume::boundaryConditions::BoundaryCondition](./ablate::finiteVolume::boundaryConditions::BoundaryCondition.html) list) the boundary conditions for the flow field

computePhysicsTimeStep
: (bool) determines if a physics based time step is used to control the FVM time stepping (default is false)

## ablate::finiteVolume::ReactingCompressibleFlowSolver
reacting compressible finite volume flow

id (req) 
: (string) the name of the flow field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the options passed to PETSc

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the TChem v1 equation of state used to describe the flow

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the compressible flow parameters cfl, gamma, etc.

transport
: ([ablate::eos::transport::TransportModel](./ablate::eos::transport::TransportModel.html)) the diffusion transport model

fluxCalculator
: ([ablate::finiteVolume::fluxCalculator::FluxCalculator](./ablate::finiteVolume::fluxCalculator::FluxCalculator.html)) the flux calculator (defaults to AUSM)

boundaryConditions
: ([ablate::finiteVolume::boundaryConditions::BoundaryCondition](./ablate::finiteVolume::boundaryConditions::BoundaryCondition.html) list) the boundary conditions for the flow field

additionalProcesses
: ([ablate::finiteVolume::processes::Process](./ablate::finiteVolume::processes::Process.html) list) any additional processes besides euler/yi/ev transport and rxn

computePhysicsTimeStep
: (bool) determines if a physics based time step is used to control the FVM time stepping (default is false)

## ablate::finiteElement::LowMachFlowSolver
incompressible FE flow

id (req) 
: (string) the name of the flow field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) options for the flow passed directly to PETSc

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the flow field parameters

boundaryConditions (req) 
: ([ablate::finiteElement::boundaryConditions::BoundaryCondition](./ablate::finiteElement::boundaryConditions::BoundaryCondition.html) list) the boundary conditions for the flow field

auxFields (req) 
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) enables and sets the update functions for the auxFields

## ablate::finiteElement::IncompressibleFlowSolver
incompressible FE flow

id (req) 
: (string) the name of the flow field

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) options for the flow passed directly to PETSc

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the flow field parameters

boundaryConditions (req) 
: ([ablate::finiteElement::boundaryConditions::BoundaryCondition](./ablate::finiteElement::boundaryConditions::BoundaryCondition.html) list) the boundary conditions for the flow field

auxFields
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) enables and sets the update functions for the auxFields

## ablate::boundarySolver::BoundarySolver
A solver used to compute boundary values in boundary cells

id (req) 
: (string) the name of the flow field

region (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.

fieldBoundary (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region describing the faces between the boundary and field

processes (req) 
: ([ablate::boundarySolver::BoundaryProcess](./ablate::boundarySolver::BoundaryProcess.html) list) a list of boundary processes

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the options passed to PETSC for the flow

mergeFaces
: (bool) determine if multiple faces should be merged for a single cell, default if false

## ablate::particles::Tracer
massless particles that advect with the flow

id (req) 
: (string) the name of this particle solver

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) options for the flow passed directly to PETSc

ndims (req) 
: (int) the number of dimensions for the particle

initializer (req) 
: ([ablate::particles::initializers::Initializer](./ablate::particles::initializers::Initializer.html)) the initial particle setup methods

exactSolution
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the particle location exact solution

## ablate::particles::Inertial
particles (with mass) that advect with the flow

id (req) 
: (string) the name of this particle solver

region
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.  Default is entire domain

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) options for the flow passed directly to PETSc

ndims (req) 
: (int) the number of dimensions for the particle

parameters (req) 
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) fluid parameters for the particles (fluidDensity, fluidViscosity, gravityField)

initializer (req) 
: ([ablate::particles::initializers::Initializer](./ablate::particles::initializers::Initializer.html)) the initial particle setup methods

fieldInitialization (req) 
: ([ablate::mathFunctions::FieldFunction](./ablate::mathFunctions::FieldFunction.html) list) the initial particle fields setup methods

exactSolution
: ([ablate::mathFunctions::MathFunction](./ablate::mathFunctions::MathFunction.html)) the particle location/velocity exact solution

## ablate::radiation::Radiation
A solver for radiative heat transfer in participating media

id (req) 
: (string) the name of the flow field

region (req) 
: ([ablate::domain::Region](./ablate::domain::Region.html)) the region to apply this solver.

rays (req) 
: (int) number of rays used by the solver

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) the options passed to PETSC for the flow

properties (req) 
: ([ablate::eos::radiationProperties::RadiationModel](./ablate::eos::radiationProperties::RadiationModel.html)) the radiation properties model

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

