#ifndef SIMPLE_ALGO_TANH_H
#define SIMPLE_ALGO_TANH_H

#include "ap_fixed.h"
#include "ap_int.h"

// size of the LUT
#define N_TABLE_SIZE 1024
#define TANH_RANGE 4

// Type used for LUT (ap_fixed<X,Y>)
#define AP_FIXED_SIZE 12
#define AP_FIXED_DEC 3
//typedef ap_fixed<AP_FIXED_SIZE,AP_FIXED_DEC> val_t;
//typedef ap_fixed<AP_FIXED_SIZE,AP_FIXED_DEC> result_t;

typedef ap_int<10> val_t;
typedef ap_int<10> result_t;
typedef ap_int<10> etaphi_t;
typedef ap_int<10> pt_t;
typedef ap_uint<10> uetaphi_t;




// reference and hardware functions
void simple_algo_tanh_ref( const float& in, float& out );
void simple_algo_tanh_hw(val_t data, result_t& res);
//void simple_algo_tanh_hw(etaphi_t data, etaphi_t& res);


// *************************************************
//       TanH Activation
// Implemented following:
//  https://github.com/Xilinx/RFNoC-HLS-NeuralNet/blob/master/rfnoc/hls/nnet_lib/nnet_activation.h#L111-L153
//  -- remove references to NN layers
//  -- Make the range for tanh (0,4) [symmetric function]
//     -- Use +/- in the function call below
// *************************************************
template<class data_T, int N_TABLE>
void init_tanh_table(data_T table_out[N_TABLE]) {
    /* Implement tanh lookup */
    for (int ii = 0; ii < N_TABLE; ii++) {
        // Convert from table index to X-value (unsigned 4-bit, range 0 to +4)
        float in_val = TANH_RANGE*((N_TABLE-1)-ii)/float(N_TABLE);

        // Next, compute lookup table function
        float freal_val = tanh(in_val)*1E3;
        data_T real_val(freal_val); // = round( tanh(in_val)*1E3 );
        //std::cout << "Tanh:  Lookup table Index: " <<  ii<< " In Value: " << in_val << " Result: " << real_val << std::endl;
        table_out[ii] = real_val;
    }

    return;
}


template<class data_T, class res_T, int TABLE_SIZE/*=1024*/>
void tanh(data_T &data, res_T &res) {
    // Initialize the lookup table
    static res_T tanh_table[TABLE_SIZE];
    init_tanh_table<res_T, TABLE_SIZE>(tanh_table);

    // Index into the lookup table based on data
    int index;

    #pragma HLS PIPELINE

    index = (1.-data/(TANH_RANGE*1E3))*TABLE_SIZE;

    if (index < 0) index = 0;
    if (index > TABLE_SIZE-1) index = TABLE_SIZE-1;
    res = tanh_table[index];

    return;
}


// Default table size provided here:
template<class data_T, class res_T>
void tanh(data_T &data, res_T &res) { 
    /* Get the tanh value from the LUT */
    if (data < 0) {
        data = -1*data;
        tanh<data_T, res_T, N_TABLE_SIZE>(data, res); 
        res  = -1*res;
    }
    else{
        tanh<data_T, res_T, N_TABLE_SIZE>(data, res); 
    }

    return;
}

#endif
