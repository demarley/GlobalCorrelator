/*
// *************************************************
//       TanH Activation
// Implemented following:
//  https://github.com/Xilinx/RFNoC-HLS-NeuralNet/blob/master/rfnoc/hls/nnet_lib/nnet_activation.h#L111-L153
//  -- remove references to NN layers
//  -- Make the range for tanh (0,4) [antisymmetric function]
//     -- Use +/- in the function call below
// *************************************************

Update bit values for FPGA outputs 
Twiki page:     https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1TriggerPhase2InterfaceSpecifications
FPGA reference: https://indico.cern.ch/event/696147/contributions/2855708/
                        attachments/1606757/2550062/TAMUHardware_SvenDildick_20180226.pdf

if z0<0: z0 += 1024 then multiply by INVZ_CONVERSION
if sinhEta<0: sinhEta+=8192 then multiply by INVETA_CONVERSION
if rinv<0: rinv+=16384 then multiply by INVRINV_CONVERSION
*/
#ifndef TK_MU_SIMPLE_H
#define TK_MU_SIMPLE_H

#include "ap_fixed.h"
#include "ap_int.h"
#include "dataformats.h"

#define DEBUG 0

namespace {

fphi_t phiOffSetValues[27] = {
  -0.0387851, // 1
  0.193925, // 2
  0.426636, // 3
  0.659347, // 4
  0.892057, // 5
  1.124768, // 6
  1.357478, // 7
  1.590189, // 8
  1.822899, // 9
  2.055610, // 10
  2.288321, // 11
  2.521031, // 12
  2.753742, // 13
  2.986452, // 14
  -3.064022, // 15
  -2.831312, // 16
  -2.598601, // 17
  -2.365891, // 18
  -2.133180, // 19
  -1.900470, // 20
  -1.667759, // 21
  -1.435049, // 22
  -1.202338, // 23
  -0.969627, // 24
  -0.736917, // 25
  -0.504206, // 26
  -0.271496  // 27
};

}

// reference and hardware functions
SwPropTrack tkmu_simple_ref( const SwTrack& in );
HwPropTrack tkmu_simple_hw (       HwTrack& in );

HwTrackMuon match_hw(const HwTrack&, const HwMuon&);
HwTrackMuon match_prop_hw(const HwPropTrack&, const HwMuon&);
SwTrackMuon match_sw(const SwTrack&, const SwMuon&);
SwTrackMuon match_prop_sw(const SwPropTrack&, const SwMuon&);

// decode track eta and phi
feta_t decode_track_eta(const HwTrack&);
fphi_t decode_track_phi(const HwTrack&);

// calculate deltaR
float deltaR(float eta1, float phi1, float eta2, float phi2);

// template functions
template<class data_T, int N_TABLE>
void init_deta_table(data_T table_out[N_TABLE]){
    /* deta_LUT  = track.hwZ0 * (1/550)*/
    for (int ii = 0; ii < N_TABLE; ii++) {
        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table 
        data_T real_val = in_val/550.;  // convert to proper type
        if (DEBUG) std::cout << "deta_LUT:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void deta_LUT(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T deta_table[TABLE_SIZE];
    init_deta_table<res_T,TABLE_SIZE>(deta_table);

    #pragma HLS PIPELINE
    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_Z0_RANGE;

    if (index<0) res = deta_table[0];
    else if (index>TABLE_SIZE-1) res = deta_table[TABLE_SIZE-1];
    else res = deta_table[index];

    return;
}

template<class data_T, class res_T>
void deta_LUT(data_T &data, res_T &res) { 
    /* Gateway to deta_LUT (z0/550) : checks boundaries */
    res = 0;
    deta_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 

    return;
}


///////////////////////////
///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_minus_LUT(data_T table_out[N_TABLE]){
    /* delta_minus_LUT  z0 / (z0-850) */
    for (int ii = 0; ii < N_TABLE; ii++) {
        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table 
        float numerator   = in_val;    // just repeat the calculation from delta_LUT
        float denominator = 850. - in_val;
        data_T real_val   = numerator/denominator; // convert to proper type
        if (DEBUG) std::cout << "delta_minus_LUT (z0/(850-z0)):  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }
}

template<class data_T, class res_T, int TABLE_SIZE>
void delta_minus_LUT(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T delta_minus_table[TABLE_SIZE];
    init_delta_minus_LUT<res_T,TABLE_SIZE>(delta_minus_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_Z0_RANGE;

    if (index<0) res = delta_minus_table[0];
    else if (index>TABLE_SIZE-1) res = delta_minus_table[TABLE_SIZE-1];
    else res = delta_minus_table[index];

    return;
}

template<class data_T, class res_T>
void delta_minus_LUT(data_T &data, res_T &res) { 
    /* Gateway to delta_minus_LUT (z0/(850-z0)) */
    res = 0;
    delta_minus_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 

    return;
}



///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_plus_LUT(data_T table_out[N_TABLE]){
    /* delta_plus_LUT  z0 / (z0+850) */
    for (int ii = 0; ii < N_TABLE; ii++) {
        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table 
        float numerator   = in_val;    // just repeat the calculation from delta_LUT
        float denominator = 850. + in_val;
        data_T real_val    = numerator/denominator; // convert to proper type
        if (DEBUG) std::cout << "delta_plus_LUT (z0/(850+z0)):  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }
}

template<class data_T, class res_T, int TABLE_SIZE>
void delta_plus_LUT(data_T &data, res_T &res) {
    /* Initialize the lookup table */
    res_T delta_plus_table[TABLE_SIZE];
    init_delta_plus_LUT<res_T,TABLE_SIZE>(delta_plus_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_Z0_RANGE;

    if (index<0) res = delta_plus_table[0];
    else if (index>TABLE_SIZE-1) res = delta_plus_table[TABLE_SIZE-1];
    else res = delta_plus_table[index];

    return;
}

template<class data_T, class res_T>
void delta_plus_LUT(data_T &data, res_T &res) { 
    /* Gateway to delta_plus_LUT (z0/(850+z0)) */
    res = 0;
    delta_plus_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 

    return;
}


///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_LUT(data_T table_out[N_TABLE]){
    /* delta_LUT  = track.hwZ0 / 850 */
    for (int ii = 0; ii < N_TABLE; ii++) {
        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table 
        data_T real_val = in_val/850.;  // convert to proper type
        if (DEBUG) std::cout << "delta_LUT (z0/850):  Lookup table Index " <<  ii << ": In Value = " << in_val << "; Result = " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void delta_LUT(data_T &data, res_T &res) {
    /* Initialize the LUT */
    res_T delta_table[TABLE_SIZE];
    init_delta_LUT<res_T,TABLE_SIZE>(delta_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_Z0_RANGE;

    if (index<0) res = delta_table[0];
    else if (index>TABLE_SIZE-1) res = delta_table[TABLE_SIZE-1];
    else res = delta_table[index];

    return;
}

template<class data_T, class res_T>
void delta_LUT(data_T &data, res_T &res) { 
    /* Gateway to delta_LUT (z0/850) : checks boundaries */
    res = 0;
    delta_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 

    return;
}


///////////////////////////
// -- TANH FUNCTION 
template<class data_T, int N_TABLE>
void init_tanh_table(data_T table_out[N_TABLE]) {
    /* Implement tanh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value (unsigned 4-bit, range 0 to +4)
        float in_val = (ETA_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = tanh(in_val);
        if (DEBUG) std::cout << "Tanh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE/*=1024*/>
void tanh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T tanh_table[TABLE_SIZE];
    init_tanh_table<res_T, TABLE_SIZE>(tanh_table);

    #pragma HLS PIPELINE

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_ETA_RANGE;

    if (index<0) res = tanh_table[0];
    else if (index>TABLE_SIZE-1) res = tanh_table[TABLE_SIZE-1];
    else res = tanh_table[index];

    return;
}

// Gateway to calling tanh(x)
template<class data_T, class res_T>
void tanh(data_T &data, res_T &res) { 
    /* Get the tanh value from the LUT -- symmetric function */
    res = 0;
    tanh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 

    return;
}


///////////////////////////
// -- 1/cosh(x) LUT (follow tanh example)
template<class data_T, int N_TABLE>
void init_cosh_table(data_T table_out[N_TABLE]) {
    /* Implement cosh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value
        float in_val = (COSH_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = 1./cosh(in_val);
        if (DEBUG) std::cout << "1/cosh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE/*=1024*/>
void invCosh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T cosh_table[TABLE_SIZE];
    init_cosh_table<res_T, TABLE_SIZE>(cosh_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_COSH_RANGE;

    if (index<0) res = cosh_table[0];
    else if (index>TABLE_SIZE-1) res = cosh_table[TABLE_SIZE-1];
    else res = cosh_table[index];

    return;
}

template<class data_T, class res_T>
void invCosh(data_T &data, res_T &res) { 
    /* Get the tanh value from the LUT -- symmetric function */
    // should only get positive values!
    res = 0;
    invCosh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 

    return;
}


///////////////////////////
// -- ARCSINH FUNCTION 
//    Tracker provides sinh(eta) and we need eta!
//    http://mathworld.wolfram.com/InverseHyperbolicSine.html
template<class data_T, int N_TABLE>
void init_arcsinh_table(data_T table_out[N_TABLE]) {
    /* Implement arcsinh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value (unsigned 4-bit, range 0 to +4)
        float in_val = (SINHETA_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = log(in_val + sqrt(1+pow(in_val,2)));
        if (DEBUG) std::cout << "Arcsinh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE/*=1024*/>
void arcsinh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T arcsinh_table[TABLE_SIZE];
    init_arcsinh_table<res_T, TABLE_SIZE>(arcsinh_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    int index = TABLE_SIZE - data * TABLE_SIZE * INV_SINHETA_RANGE;

    if (index<0) res = arcsinh_table[0];
    else if (index>TABLE_SIZE-1) res = arcsinh_table[TABLE_SIZE-1];
    else res = arcsinh_table[index];

    return;
}

// Gateway to calling arcsinh(x)
template<class data_T, class res_T>
void arcsinh(data_T &data, res_T &res) { 
    /* Get the arcsinh value from the LUT -- anti-symmetric function */
    res = 0;
    arcsinh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 

    return;
}

// delta R matching
// choose eta precision of the muon as the dR precision of the match
template<class data_T, class data_S, class data_U, class data_V>
data_S dr2_int(data_T eta1, data_S phi1, data_U eta2, data_V phi2) {
  // eta1, phi1: track properties
  // eta2, phi2: muon properties
  data_T deta = (eta1-eta2);
  data_S dphi = (phi1-phi2);
  data_S dR = deta*deta + dphi*dphi;
  bool debug(false);
  if (debug) {
    std::cout << "eta1 " << eta1 << std::endl;
    std::cout << "eta2 " << eta2 << std::endl;
    std::cout << "phi1 " << phi1 << std::endl;
    std::cout << "phi2 " << phi2 << std::endl;
    std::cout << "dR " << dR << std::endl << std::endl;
  }
  return dR;
}

#endif
