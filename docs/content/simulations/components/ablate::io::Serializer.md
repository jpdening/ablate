---
layout: default
title: ablate::io::Serializer
parent: Component List
grand_parent: Running Simulations
nav_exclude: true
---
# [ablate::io::Serializer](./ablate::io::Serializer.html)
## ablate::io::Hdf5Serializer*
default serializer for IO

interval (req) 
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) The interval object used to determine write interval.

## ablate::io::Hdf5MultiFileSerializer
serializer for IO that writes each time to a separate hdf5 file

interval (req) 
: ([ablate::io::interval::Interval](./ablate::io::interval::Interval.html)) The interval object used to determine write interval.

options
: ([ablate::parameters::Parameters](./ablate::parameters::Parameters.html)) options for the viewer passed directly to PETSc including (hdf5ViewerView, viewer_hdf5_collective, viewer_hdf5_sp_output

