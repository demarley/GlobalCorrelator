############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2016 Xilinx, Inc. All Rights Reserved.
##
## source setup.sh
## vivado_hls run_hls_tanh.tcl 
############################################################

# open the project, don't forget to reset
open_project -reset hls_match

set_top match_hw
add_files src/tk-mu_simple.cpp
add_files src/dataformats.h

## test-bench files
add_files -tb tk-mu_simple_test.cpp
add_files -tb tk-mu_simple_ref.h
add_files -tb src/dataformatsExtended.h

# reset the solution
open_solution -reset "solution1"
 #set_part {xcku9p-ffve900-2-i-EVAL}
set_part {xc7vx690tffg1927-2}
 #set_part {xcku5p-sfvb784-3-e}
 #set_part {xcku115-flvf1924-2-i}
create_clock -period 5 -name default
 #source "./nb1/solution1/directives.tcl"

# do stuff
#csim_design
csynth_design
 #cosim_design -trace_level all
export_design -format ip_catalog  -vendor "cern-cms"

# exit Vivado HLS
exit
