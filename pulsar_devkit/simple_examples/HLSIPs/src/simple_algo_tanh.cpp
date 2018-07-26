/*
HLS implementation of TANH function via LUT

https://github.com/p2l1pfp/GlobalCorrelator_HLS/blob/dev/firmware/data.h
 typedef ap_int<16> pt_t;
 typedef ap_int<10>  etaphi_t;
 typedef ap_int<5>  vtx_t;
 typedef ap_uint<3>  particleid_t;
 typedef ap_int<10> z0_t;  // 40cm / 0.1
*/
#include "simple_algo_tanh.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif

// https://github.com/Xilinx/RFNoC-HLS-NeuralNet/blob/master/rfnoc/hls/test_activations/test_activations.h/cpp
void simple_algo_tanh_hw(val_t data, result_t& res){ 
//void simple_algo_tanh_hw(etaphi_t data, etaphi_t& res){ 
    tanh<val_t, result_t>(data, res);
//    tanh<etaphi_t, etaphi_t>(data, res);
    return;
}

// THE END
