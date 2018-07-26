# GlobalCorrelator
Firmware for Level-1 Particle Reconstruction


## Track-muon linking

This setup is designed for the TAMU cmstca machine.  
To run the 'test' setup for track-muon linking, follow these instructions:

```
# Move to the right directory and setup the environment
cd pulsar_devkit/simple_examples/HLSIPs/
source setup.sh
```

There is one macro needed to run the algorithm, `run_hls_tkmu.tcl`.  
Feasibly this could all be run from the GUI as well.  
This processes both a C++ (floating-point) and HLS algorithm.

```
vivado_hls run_hls_tkmu.tcl
```

The output will contain information from HLS (Vivado) and `cout` statements from the framework.
The relevant scripts:

- `tk-mu_simple_test.cpp`: Main script run by the framework
- `tk-mu_simple_ref.cpp`: C++ floating-point algorithm
- `tk-mu_simple.h`: Contains struct definitions, LUTs, etc.
- `tk-mu_simple.cpp`: HLS implementation of the C++ track propagation algorithm
