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
#include <stdlib.h>

#define DEBUG 0

typedef ap_int<15> invpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_int<12> pt_t;     // convert from RINV
typedef ap_int<14> eta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_int<19> phi_t;    // phi (50 micro-rad)
typedef ap_int<10> chisq_t;  // chi^2 (0 - 100; 0.1 steps)
typedef ap_int<1> q_t;       // charge
typedef ap_int<11> z0_t;     // z0  (1 mm over +/-14.9 cm)
typedef ap_int<5> sector_t;

typedef ap_fixed<15,2> finvpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_fixed<12,9> fpt_t;
typedef ap_fixed<14,4> feta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<19,3> fphi_t;    // phi (50 micro-rad)
typedef ap_fixed<10,7> fchisq_t;  // chi^2 (0 - 100; 0.1 steps) 
typedef ap_fixed<11,5> fz0_t;     // z0  (1 mm over +/-14.9 cm) 


// size of the LUTs
#define ETA_TABLE_SIZE 8192  // 13 unsigned bits
#define Z0_TABLE_SIZE 1024   // 10 unsigned bits

/*
// range for LUTs
#define SINHETA_RANGE 6
#define ETA_RANGE 3
#define COSH_RANGE 3
#define Z0_RANGE 15
#define INV_SINHETA_RANGE 1/SINHETA_RANGE
#define INV_ETA_RANGE 1/ETA_RANGE
#define INV_COSH_RANGE 1/COSH_RANGE
#define INV_Z0_RANGE 1/Z0_RANGE

// Conversions between binary and floating point (using example file to derive)
#define RINV_CONVERSION 792055              //1314229             // 1/0.000000760902077
#define PT_CONVERSION 87719298E-6           // 1/(0.01*0.3*3.8); 87719298E-6
#define ETA_CONVERSION 512                  //855   // 1/0.0011698 = 854.84698
#define PHI_CONVERSION 211233               //original: 219037
#define Z_CONVERSION 17                     //original: 18 // 1/0.05615 = 17.81 -> 18
#define INVRINV_CONVERSION 1262538462E-15   //0.000001262538462  //original: 760902077E-15    // 0.000000760902077
#define INVETA_CONVERSION 19531261E-10      //original: 11698E-7
#define INVPHI_CONVERSION 4734119709E-15    //0.000004734119709  // original: 456544E-11
#define INVZ_CONVERSION 5859375E-8          //0.05859375 //original: 56152375E-9         //0.056152375 -- shift of 1024 for negative values!
*/

// -- Define structs for physics objects in software
struct TrackObj_tkmu {
    float pt;
    float eta;
    float phi;
    float z0;
    int q;
    int VALID;
    int BX;
};

struct PropTrackObj_tkmu : TrackObj_tkmu {
    float propEta;
    float propPhi;
};

// -- Define structs for physics objects in hardware
struct TkObj_tkmu {
    invpt_t hwRinv;
    invpt_t hwInvPt;
    pt_t hwPt;
    eta_t hwSinhEta;
    eta_t hwEta;
    phi_t hwPhi;
    eta_t hwZ0;  // same precision at eta_t
    q_t hwQ;
    chisq_t hwX2;
    q_t VALID;   // VALID bit
    fpt_t BX;    // bunch crossing
    sector_t hwSector;
};

struct PropTkObj_tkmu : TkObj_tkmu {
    eta_t hwPropEta;
    phi_t hwPropPhi;
};

inline void clear(TkObj_tkmu & c) {
    c.hwRinv  = 0; 
    c.hwInvPt = 0; 
    c.hwEta = 0; 
    c.hwPhi = 0; 
    c.hwSinhEta = 0; 
    c.hwPt  = 0;
    c.hwZ0  = 0;
    c.hwQ   = 0;
    c.hwX2  = 0;
    c.VALID = 0;
    c.BX = 0;
}
inline void init(PropTkObj_tkmu & p, const TkObj_tkmu & i){
    p.hwRinv  = i.hwRinv; 
    p.hwInvPt = i.hwInvPt; 
    p.hwEta = i.hwEta; 
    p.hwPhi = i.hwPhi; 
    p.hwSinhEta = i.hwSinhEta; 
    p.hwPt  = i.hwPt;
    p.hwZ0  = i.hwZ0;
    p.hwQ   = i.hwQ;
    p.hwX2  = i.hwX2;
    p.VALID = i.VALID;
    p.BX = i.BX;
}
inline void clearProp(PropTkObj_tkmu & c) {
    clear(c);
    c.hwPropEta = 0; 
    c.hwPropPhi = 0; 
}


// reference and hardware functions
void tkmu_simple_ref( const TrackObj_tkmu& in, PropTrackObj_tkmu& out );   // software implementation
//void tkmu_simple_hw(  TkObj_tkmu& in, PropTkObj_tkmu& out );
PropTkObj_tkmu tkmu_simple_hw(  TkObj_tkmu& in );


//////////////
//// LUTS ////
//////////////

template<class data_T, int N_TABLE>
void init_deta_table(data_T table_out[N_TABLE]){
    /* deta_LUT  = track.hwZ0 * (1/550)*/
    for (int ii = 0; ii < N_TABLE; ii++) {
//        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (15)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table 
        data_T real_val = 0.002*in_val;  // convert to proper type
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
//    data_T tmp(Z0_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/15;   //tmp;

    if (index<0) res = deta_table[0];
    else if (index>TABLE_SIZE-1) res = deta_table[TABLE_SIZE-1];
    else res = deta_table[index];

    return;
}

template<class data_T, class res_T>
void deta_LUT(data_T &data, res_T &res) { 
    /* Gateway to deta_LUT (z0/550) : checks boundaries */
    res = 0;
//    deta_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 
    deta_LUT<data_T, res_T, 1024>(data, res); 

    return;
}


///////////////////////////
///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_minus_LUT(data_T table_out[N_TABLE]){
    /* delta_minus_LUT  z0 / (z0-850) */
    for (int ii = 0; ii < N_TABLE; ii++) {
//        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (15)*((N_TABLE-1)-ii)/float(N_TABLE);

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
//    data_T tmp(Z0_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/15;   // tmp;

    if (index<0) res = delta_minus_table[0];
    else if (index>TABLE_SIZE-1) res = delta_minus_table[TABLE_SIZE-1];
    else res = delta_minus_table[index];

    return;
}

template<class data_T, class res_T>
void delta_minus_LUT(data_T &data, res_T &res) { 
    /* Gateway to delta_minus_LUT (z0/(850-z0)) */
    res = 0;
//    delta_minus_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 
    delta_minus_LUT<data_T, res_T, 1024>(data, res); 

    return;
}



///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_plus_LUT(data_T table_out[N_TABLE]){
    /* delta_plus_LUT  z0 / (z0+850) */
    for (int ii = 0; ii < N_TABLE; ii++) {
//        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (15)*((N_TABLE-1)-ii)/float(N_TABLE);

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
//    data_T tmp(Z0_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/15;   //tmp;

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
    delta_plus_LUT<data_T, res_T, 1024>(data, res); 

    return;
}


///////////////////////////
template<class data_T, int N_TABLE>
void init_delta_LUT(data_T table_out[N_TABLE]){
    /* delta_LUT  = track.hwZ0 / 850 */
    for (int ii = 0; ii < N_TABLE; ii++) {
//        float in_val = (Z0_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (15)*((N_TABLE-1)-ii)/float(N_TABLE);

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
//    data_T tmp(Z0_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/15;   //tmp;

    if (index<0) res = delta_table[0];
    else if (index>TABLE_SIZE-1) res = delta_table[TABLE_SIZE-1];
    else res = delta_table[index];

    return;
}

template<class data_T, class res_T>
void delta_LUT(data_T &data, res_T &res) { 
    /* Gateway to delta_LUT (z0/850) : checks boundaries */
    res = 0;
//    delta_LUT<data_T, res_T, Z0_TABLE_SIZE>(data, res); 
    delta_LUT<data_T, res_T, 1024>(data, res); 

    return;
}


///////////////////////////
// -- TANH FUNCTION 
template<class data_T, int N_TABLE>
void init_tanh_table(data_T table_out[N_TABLE]) {
    /* Implement tanh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value (unsigned 4-bit, range 0 to +4)
//        float in_val = (ETA_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (3)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = tanh(in_val);
        if (DEBUG) std::cout << "Tanh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void tanh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T tanh_table[TABLE_SIZE];
    init_tanh_table<res_T, TABLE_SIZE>(tanh_table);

    #pragma HLS PIPELINE

    // convert input to index
    //data_T tmp(ETA_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/3;   //tmp;

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
//    tanh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 
    tanh<data_T, res_T, 8192>(data, res); 

    return;
}


///////////////////////////
// -- 1/cosh(x) LUT (follow tanh example)
template<class data_T, int N_TABLE>
void init_cosh_table(data_T table_out[N_TABLE]) {
    /* Implement cosh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value
//        float in_val = (COSH_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = (3)*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = 1./cosh(in_val);
        if (DEBUG) std::cout << "1/cosh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void invCosh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T cosh_table[TABLE_SIZE];
    init_cosh_table<res_T, TABLE_SIZE>(cosh_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
//    data_T tmp(COSH_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/3;   //tmp;

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
//    invCosh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 
    invCosh<data_T, res_T, 8192>(data, res); 

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
//        float in_val = (SINHETA_RANGE)*((N_TABLE-1)-ii)/float(N_TABLE);
        float in_val = 6*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        data_T real_val = log(in_val + sqrt(1+pow(in_val,2)));
        if (DEBUG) std::cout << "Arcsinh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}

template<class data_T, class res_T, int TABLE_SIZE>
void arcsinh(data_T &data, res_T &res) {
    // Initialize the lookup table
    res_T arcsinh_table[TABLE_SIZE];
    init_arcsinh_table<res_T, TABLE_SIZE>(arcsinh_table);

    #pragma HLS PIPELINE

    res = 0;

    // convert input to index
    // data_T tmp(SINHETA_RANGE);
    int index = TABLE_SIZE - data * TABLE_SIZE/6;  // TABLE_SIZE/tmp;

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
//    arcsinh<data_T, res_T, ETA_TABLE_SIZE>(data, res); 
    arcsinh<data_T, res_T, 8192>(data, res); 

    return;
}

#endif
