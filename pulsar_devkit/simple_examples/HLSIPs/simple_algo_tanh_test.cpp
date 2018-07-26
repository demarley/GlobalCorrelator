/*
Example using tanh LUT
*/
#include <vector>
#include <cstdio>

#include "ap_fixed.h"
#include "src/simple_algo_tanh.h"


int main() {
    //etaphi_t in_hw;
    //etaphi_t out_hw;
    val_t in_hw;
    val_t out_hw;
//    result_t out_hw;

    float in,out;
    float values = 0.1;

    in  = 0.1;
    out = 0;

    // c++ implementation
    simple_algo_tanh_ref(in, out);
    std::cout << " REF : tanh(" << in << ") = " << out << std::endl;


    // hardware implementation

    in_hw  = values*1000;   //round(values*1E3);
    out_hw = 0;
    simple_algo_tanh_hw(in_hw, out_hw);

    std::cout << " HW  : tanh(" << in_hw << ") = " << out_hw << std::endl;

    return 0;
}
