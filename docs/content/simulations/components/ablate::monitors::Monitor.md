---
layout: default
title: ablate::monitors::Monitor
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::monitors::Monitor](./ablate::monitors::Monitor.html)
## ablate::monitors::FieldErrorMonitor
Computes and reports the error every time step

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

## ablate::monitors::SolutionErrorMonitor
Computes and reports the error every time step

scope (req) 
: (ablate::monitors::SolutionErrorMonitor::Scope enum) how the error should be calculated ('vector', 'component')

type (req) 
: (ablate::monitors::SolutionErrorMonitor::Norm enum) norm type ('l1','l1_norm','l2', 'linf', 'l2_norm')

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

## ablate::monitors::TimeStepMonitor
Reports the current step, time, and dt

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

interval
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) report interval object, defaults to every

## ablate::monitors::IgnitionDelayPeakYi
Compute the ignition time based upon peak mass fraction

species (req) 
: (string) the species used to determine the peak Yi

location (req) 
: (double list) the monitor location

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record the final ignition time (default is stdout)

historyLog
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record the time and yi history (default is none)

## ablate::monitors::IgnitionDelayTemperature
Compute the ignition time based upon temperature change

eos (req) 
: ([ablate::eos::EOS](./ablate::eos::EOS.html)) the eos used to compute temperature

location (req) 
: (double list) the monitor location

thresholdTemperature (req) 
: (double) the temperature used to define ignition delay

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record the final ignition time (default is stdout)

historyLog
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record the time and yi history (default is none)

## ablate::monitors::ExtractLineMonitor
Outputs the results along a line as a curve file (beta)

interval (req) 
: (int) output interval

prefix (req) 
: (string) the file prefix

start (req) 
: (double list) the line start location

end (req) 
: (double list) the line end location

outputFields (req) 
: (string list) a list of fields to write to the curve

outputAuxFields (req) 
: (string list) a list of aux fields to write to the curve 

## ablate::monitors::DmViewFromOptions
replicates the DMViewFromOptions function in PETSC

scope (req) 
: (ablate::monitors::DmViewFromOptions::Scope enum) determines if DMViewFromOptions is called initially (initial) or every time step (monitor)

options
: (string) if provided these options are used for the DMView call, otherwise global options is used

optionName
: (string) if provided the optionsName is used for DMViewFromOptions.  Needed if using global options.

## ablate::monitors::ParticleCount
Outputs the total number of particles in the domain

interval
: (int) output interval

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

## ablate::monitors::ParticleAverage
Outputs the average particle location in the domain

interval
: (int) output interval

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

## ablate::monitors::CurveMonitor
Write 1D results to a curve file

interval
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) output interval

prefix
: (string) the file prefix

## ablate::monitors::MaxMinAverage
Prints the min/max/average for a field

field (req) 
: (string) the name of the field

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

interval
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) report interval object, defaults to every

## ablate::monitors::PhysicsTimeStep
Reports the physics based time stepping from the FVM without enforcing it

log
: ([ablate::monitors::logs::Log](./ablate::monitors::logs::Log.html)) where to record log (default is stdout)

interval
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) report interval object, defaults to every

## ablate::monitors::Probes
Records the values of the specified variables at a specific point in space

probes (req) 
: ([ablate::monitors::probes::ProbeInitializer](./ablate::monitors::probes::ProbeInitializer.html)) where to record log (default is stdout)

variables (req) 
: (string list) list of variables to output

interval
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) report interval object, defaults to every

bufferSize
: (int) how often the probe file is written (default is 100, must be > 0)

