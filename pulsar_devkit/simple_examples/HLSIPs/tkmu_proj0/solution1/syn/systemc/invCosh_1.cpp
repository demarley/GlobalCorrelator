// ==============================================================
// RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
// Version: 2016.4
// Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
// 
// ===========================================================

#include "invCosh_1.h"
#include "AESL_pkg.h"

using namespace std;

namespace ap_rtl {

const sc_logic invCosh_1::ap_const_logic_1 = sc_dt::Log_1;
const sc_logic invCosh_1::ap_const_logic_0 = sc_dt::Log_0;
const sc_lv<1> invCosh_1::ap_ST_fsm_pp0_stage0 = "1";
const sc_lv<32> invCosh_1::ap_const_lv32_0 = "00000000000000000000000000000000";
const sc_lv<1> invCosh_1::ap_const_lv1_1 = "1";
const sc_lv<13> invCosh_1::ap_const_lv13_0 = "0000000000000";
const sc_lv<49> invCosh_1::ap_const_lv49_1555556 = "1010101010101010101010110";
const sc_lv<32> invCosh_1::ap_const_lv32_1A = "11010";
const sc_lv<32> invCosh_1::ap_const_lv32_30 = "110000";
const sc_lv<25> invCosh_1::ap_const_lv25_800000 = "100000000000000000000000";
const sc_lv<32> invCosh_1::ap_const_lv32_A = "1010";
const sc_lv<32> invCosh_1::ap_const_lv32_18 = "11000";
const sc_lv<32> invCosh_1::ap_const_lv32_17 = "10111";
const sc_lv<2> invCosh_1::ap_const_lv2_0 = "00";
const sc_lv<11> invCosh_1::ap_const_lv11_400 = "10000000000";

invCosh_1::invCosh_1(sc_module_name name) : sc_module(name), mVcdFile(0) {
    cosh_table7_U = new invCosh_1_cosh_tajbC("cosh_table7_U");
    cosh_table7_U->clk(ap_clk);
    cosh_table7_U->reset(ap_rst);
    cosh_table7_U->address0(cosh_table7_address0);
    cosh_table7_U->ce0(cosh_table7_ce0);
    cosh_table7_U->q0(cosh_table7_q0);
    tkmu_simple_hw_mukbM_U19 = new tkmu_simple_hw_mukbM<1,5,26,24,49>("tkmu_simple_hw_mukbM_U19");
    tkmu_simple_hw_mukbM_U19->clk(ap_clk);
    tkmu_simple_hw_mukbM_U19->reset(ap_rst);
    tkmu_simple_hw_mukbM_U19->din0(grp_fu_80_p0);
    tkmu_simple_hw_mukbM_U19->din1(grp_fu_80_p1);
    tkmu_simple_hw_mukbM_U19->ce(ap_var_for_const0);
    tkmu_simple_hw_mukbM_U19->dout(grp_fu_80_p2);

    SC_METHOD(thread_ap_clk_no_reset_);
    dont_initialize();
    sensitive << ( ap_clk.pos() );

    SC_METHOD(thread_ap_CS_fsm_pp0_stage0);
    sensitive << ( ap_CS_fsm );

    SC_METHOD(thread_ap_done);
    sensitive << ( ap_start );
    sensitive << ( ap_CS_fsm_pp0_stage0 );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter6 );

    SC_METHOD(thread_ap_enable_reg_pp0_iter0);
    sensitive << ( ap_start );

    SC_METHOD(thread_ap_idle);
    sensitive << ( ap_start );
    sensitive << ( ap_CS_fsm_pp0_stage0 );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter1 );
    sensitive << ( ap_enable_reg_pp0_iter2 );
    sensitive << ( ap_enable_reg_pp0_iter3 );
    sensitive << ( ap_enable_reg_pp0_iter4 );
    sensitive << ( ap_enable_reg_pp0_iter5 );
    sensitive << ( ap_enable_reg_pp0_iter6 );

    SC_METHOD(thread_ap_pipeline_idle_pp0);
    sensitive << ( ap_start );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter1 );
    sensitive << ( ap_enable_reg_pp0_iter2 );
    sensitive << ( ap_enable_reg_pp0_iter3 );
    sensitive << ( ap_enable_reg_pp0_iter4 );
    sensitive << ( ap_enable_reg_pp0_iter5 );

    SC_METHOD(thread_ap_ready);
    sensitive << ( ap_start );
    sensitive << ( ap_CS_fsm_pp0_stage0 );
    sensitive << ( ap_enable_reg_pp0_iter0 );

    SC_METHOD(thread_ap_return);
    sensitive << ( ap_start );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter6 );
    sensitive << ( cosh_table7_q0 );
    sensitive << ( ap_pipeline_reg_pp0_iter5_icmp_reg_157 );

    SC_METHOD(thread_cosh_table7_address0);
    sensitive << ( ap_enable_reg_pp0_iter5 );
    sensitive << ( tmp_196_fu_135_p1 );

    SC_METHOD(thread_cosh_table7_ce0);
    sensitive << ( ap_start );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter5 );

    SC_METHOD(thread_grp_fu_80_p0);
    sensitive << ( ap_CS_fsm_pp0_stage0 );
    sensitive << ( ap_enable_reg_pp0_iter0 );

    SC_METHOD(thread_grp_fu_80_p1);
    sensitive << ( ap_CS_fsm_pp0_stage0 );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( grp_fu_80_p10 );

    SC_METHOD(thread_grp_fu_80_p10);
    sensitive << ( r_V_tr_fu_68_p3 );

    SC_METHOD(thread_icmp_fu_126_p2);
    sensitive << ( ap_start );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_enable_reg_pp0_iter4 );
    sensitive << ( tmp_233_fu_116_p4 );

    SC_METHOD(thread_r_V_fu_100_p2);
    sensitive << ( tmp_cast_cast_fu_96_p1 );

    SC_METHOD(thread_r_V_tr_fu_68_p3);
    sensitive << ( tmp_fu_64_p1 );

    SC_METHOD(thread_tmp_196_fu_135_p1);
    sensitive << ( tmp_198_fu_132_p1 );

    SC_METHOD(thread_tmp_198_fu_132_p1);
    sensitive << ( tmp_232_reg_152 );

    SC_METHOD(thread_tmp_231_fu_86_p4);
    sensitive << ( grp_fu_80_p2 );

    SC_METHOD(thread_tmp_233_fu_116_p4);
    sensitive << ( r_V_fu_100_p2 );

    SC_METHOD(thread_tmp_cast_cast_fu_96_p1);
    sensitive << ( tmp_231_fu_86_p4 );

    SC_METHOD(thread_tmp_fu_64_p1);
    sensitive << ( data_V_read );

    SC_METHOD(thread_ap_NS_fsm);
    sensitive << ( ap_start );
    sensitive << ( ap_CS_fsm );
    sensitive << ( ap_enable_reg_pp0_iter0 );
    sensitive << ( ap_pipeline_idle_pp0 );

    SC_THREAD(thread_ap_var_for_const0);

    ap_CS_fsm = "1";
    ap_enable_reg_pp0_iter1 = SC_LOGIC_0;
    ap_enable_reg_pp0_iter2 = SC_LOGIC_0;
    ap_enable_reg_pp0_iter3 = SC_LOGIC_0;
    ap_enable_reg_pp0_iter4 = SC_LOGIC_0;
    ap_enable_reg_pp0_iter5 = SC_LOGIC_0;
    ap_enable_reg_pp0_iter6 = SC_LOGIC_0;
    static int apTFileNum = 0;
    stringstream apTFilenSS;
    apTFilenSS << "invCosh_1_sc_trace_" << apTFileNum ++;
    string apTFn = apTFilenSS.str();
    mVcdFile = sc_create_vcd_trace_file(apTFn.c_str());
    mVcdFile->set_time_unit(1, SC_PS);
    if (1) {
#ifdef __HLS_TRACE_LEVEL_PORT_HIER__
    sc_trace(mVcdFile, ap_clk, "(port)ap_clk");
    sc_trace(mVcdFile, ap_rst, "(port)ap_rst");
    sc_trace(mVcdFile, ap_start, "(port)ap_start");
    sc_trace(mVcdFile, ap_done, "(port)ap_done");
    sc_trace(mVcdFile, ap_idle, "(port)ap_idle");
    sc_trace(mVcdFile, ap_ready, "(port)ap_ready");
    sc_trace(mVcdFile, data_V_read, "(port)data_V_read");
    sc_trace(mVcdFile, ap_return, "(port)ap_return");
#endif
#ifdef __HLS_TRACE_LEVEL_INT__
    sc_trace(mVcdFile, ap_CS_fsm, "ap_CS_fsm");
    sc_trace(mVcdFile, ap_CS_fsm_pp0_stage0, "ap_CS_fsm_pp0_stage0");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter0, "ap_enable_reg_pp0_iter0");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter1, "ap_enable_reg_pp0_iter1");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter2, "ap_enable_reg_pp0_iter2");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter3, "ap_enable_reg_pp0_iter3");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter4, "ap_enable_reg_pp0_iter4");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter5, "ap_enable_reg_pp0_iter5");
    sc_trace(mVcdFile, ap_enable_reg_pp0_iter6, "ap_enable_reg_pp0_iter6");
    sc_trace(mVcdFile, cosh_table7_address0, "cosh_table7_address0");
    sc_trace(mVcdFile, cosh_table7_ce0, "cosh_table7_ce0");
    sc_trace(mVcdFile, cosh_table7_q0, "cosh_table7_q0");
    sc_trace(mVcdFile, tmp_232_reg_152, "tmp_232_reg_152");
    sc_trace(mVcdFile, icmp_fu_126_p2, "icmp_fu_126_p2");
    sc_trace(mVcdFile, icmp_reg_157, "icmp_reg_157");
    sc_trace(mVcdFile, ap_pipeline_reg_pp0_iter5_icmp_reg_157, "ap_pipeline_reg_pp0_iter5_icmp_reg_157");
    sc_trace(mVcdFile, tmp_196_fu_135_p1, "tmp_196_fu_135_p1");
    sc_trace(mVcdFile, tmp_fu_64_p1, "tmp_fu_64_p1");
    sc_trace(mVcdFile, r_V_tr_fu_68_p3, "r_V_tr_fu_68_p3");
    sc_trace(mVcdFile, grp_fu_80_p0, "grp_fu_80_p0");
    sc_trace(mVcdFile, grp_fu_80_p1, "grp_fu_80_p1");
    sc_trace(mVcdFile, grp_fu_80_p2, "grp_fu_80_p2");
    sc_trace(mVcdFile, tmp_231_fu_86_p4, "tmp_231_fu_86_p4");
    sc_trace(mVcdFile, tmp_cast_cast_fu_96_p1, "tmp_cast_cast_fu_96_p1");
    sc_trace(mVcdFile, r_V_fu_100_p2, "r_V_fu_100_p2");
    sc_trace(mVcdFile, tmp_233_fu_116_p4, "tmp_233_fu_116_p4");
    sc_trace(mVcdFile, tmp_198_fu_132_p1, "tmp_198_fu_132_p1");
    sc_trace(mVcdFile, ap_NS_fsm, "ap_NS_fsm");
    sc_trace(mVcdFile, ap_pipeline_idle_pp0, "ap_pipeline_idle_pp0");
    sc_trace(mVcdFile, grp_fu_80_p10, "grp_fu_80_p10");
#endif

    }
}

invCosh_1::~invCosh_1() {
    if (mVcdFile) 
        sc_close_vcd_trace_file(mVcdFile);

    delete cosh_table7_U;
    delete tkmu_simple_hw_mukbM_U19;
}

void invCosh_1::thread_ap_var_for_const0() {
    ap_var_for_const0 = ap_const_logic_1;
}

void invCosh_1::thread_ap_clk_no_reset_() {
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_CS_fsm = ap_ST_fsm_pp0_stage0;
    } else {
        ap_CS_fsm = ap_NS_fsm.read();
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter1 = ap_const_logic_0;
    } else {
        if ((esl_seteq<1,1,1>(ap_CS_fsm_pp0_stage0.read(), ap_const_lv1_1) && 
             !(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0)))) {
            ap_enable_reg_pp0_iter1 = ap_start.read();
        }
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter2 = ap_const_logic_0;
    } else {
        if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
            ap_enable_reg_pp0_iter2 = ap_enable_reg_pp0_iter1.read();
        }
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter3 = ap_const_logic_0;
    } else {
        if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
            ap_enable_reg_pp0_iter3 = ap_enable_reg_pp0_iter2.read();
        }
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter4 = ap_const_logic_0;
    } else {
        if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
            ap_enable_reg_pp0_iter4 = ap_enable_reg_pp0_iter3.read();
        }
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter5 = ap_const_logic_0;
    } else {
        if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
            ap_enable_reg_pp0_iter5 = ap_enable_reg_pp0_iter4.read();
        }
    }
    if ( ap_rst.read() == ap_const_logic_1) {
        ap_enable_reg_pp0_iter6 = ap_const_logic_0;
    } else {
        if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
            ap_enable_reg_pp0_iter6 = ap_enable_reg_pp0_iter5.read();
        }
    }
    if (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0))) {
        ap_pipeline_reg_pp0_iter5_icmp_reg_157 = icmp_reg_157.read();
        icmp_reg_157 = icmp_fu_126_p2.read();
        tmp_232_reg_152 = r_V_fu_100_p2.read().range(24, 10);
    }
}

void invCosh_1::thread_ap_CS_fsm_pp0_stage0() {
    ap_CS_fsm_pp0_stage0 = ap_CS_fsm.read().range(0, 0);
}

void invCosh_1::thread_ap_done() {
    if (((esl_seteq<1,1,1>(ap_const_logic_0, ap_start.read()) && 
          esl_seteq<1,1,1>(ap_CS_fsm_pp0_stage0.read(), ap_const_lv1_1) && 
          esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read())) || 
         (!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0)) && 
          esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter6.read())))) {
        ap_done = ap_const_logic_1;
    } else {
        ap_done = ap_const_logic_0;
    }
}

void invCosh_1::thread_ap_enable_reg_pp0_iter0() {
    ap_enable_reg_pp0_iter0 = ap_start.read();
}

void invCosh_1::thread_ap_idle() {
    if ((esl_seteq<1,1,1>(ap_const_logic_0, ap_start.read()) && 
         esl_seteq<1,1,1>(ap_CS_fsm_pp0_stage0.read(), ap_const_lv1_1) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter0.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter1.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter2.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter3.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter4.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter5.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter6.read()))) {
        ap_idle = ap_const_logic_1;
    } else {
        ap_idle = ap_const_logic_0;
    }
}

void invCosh_1::thread_ap_pipeline_idle_pp0() {
    if ((esl_seteq<1,1,1>(ap_const_logic_0, ap_start.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter0.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter1.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter2.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter3.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter4.read()) && 
         esl_seteq<1,1,1>(ap_const_logic_0, ap_enable_reg_pp0_iter5.read()))) {
        ap_pipeline_idle_pp0 = ap_const_logic_1;
    } else {
        ap_pipeline_idle_pp0 = ap_const_logic_0;
    }
}

void invCosh_1::thread_ap_ready() {
    if ((esl_seteq<1,1,1>(ap_CS_fsm_pp0_stage0.read(), ap_const_lv1_1) && 
         esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && 
         !(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0)))) {
        ap_ready = ap_const_logic_1;
    } else {
        ap_ready = ap_const_logic_0;
    }
}

void invCosh_1::thread_ap_return() {
    ap_return = (!ap_pipeline_reg_pp0_iter5_icmp_reg_157.read()[0].is_01())? sc_lv<11>(): ((ap_pipeline_reg_pp0_iter5_icmp_reg_157.read()[0].to_bool())? ap_const_lv11_400: cosh_table7_q0.read());
}

void invCosh_1::thread_cosh_table7_address0() {
    cosh_table7_address0 =  (sc_lv<13>) (tmp_196_fu_135_p1.read());
}

void invCosh_1::thread_cosh_table7_ce0() {
    if ((!(esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter0.read()) && esl_seteq<1,1,1>(ap_start.read(), ap_const_logic_0)) && 
         esl_seteq<1,1,1>(ap_const_logic_1, ap_enable_reg_pp0_iter5.read()))) {
        cosh_table7_ce0 = ap_const_logic_1;
    } else {
        cosh_table7_ce0 = ap_const_logic_0;
    }
}

void invCosh_1::thread_grp_fu_80_p0() {
    grp_fu_80_p0 =  (sc_lv<26>) (ap_const_lv49_1555556);
}

void invCosh_1::thread_grp_fu_80_p1() {
    grp_fu_80_p1 =  (sc_lv<24>) (grp_fu_80_p10.read());
}

void invCosh_1::thread_grp_fu_80_p10() {
    grp_fu_80_p10 = esl_zext<49,24>(r_V_tr_fu_68_p3.read());
}

void invCosh_1::thread_icmp_fu_126_p2() {
    icmp_fu_126_p2 = (!tmp_233_fu_116_p4.read().is_01() || !ap_const_lv2_0.is_01())? sc_lv<1>(): sc_lv<1>(tmp_233_fu_116_p4.read() != ap_const_lv2_0);
}

void invCosh_1::thread_r_V_fu_100_p2() {
    r_V_fu_100_p2 = (!ap_const_lv25_800000.is_01() || !tmp_cast_cast_fu_96_p1.read().is_01())? sc_lv<25>(): (sc_biguint<25>(ap_const_lv25_800000) - sc_biguint<25>(tmp_cast_cast_fu_96_p1.read()));
}

void invCosh_1::thread_r_V_tr_fu_68_p3() {
    r_V_tr_fu_68_p3 = esl_concat<11,13>(tmp_fu_64_p1.read(), ap_const_lv13_0);
}

void invCosh_1::thread_tmp_196_fu_135_p1() {
    tmp_196_fu_135_p1 = esl_zext<64,16>(tmp_198_fu_132_p1.read());
}

void invCosh_1::thread_tmp_198_fu_132_p1() {
    tmp_198_fu_132_p1 = esl_sext<16,15>(tmp_232_reg_152.read());
}

void invCosh_1::thread_tmp_231_fu_86_p4() {
    tmp_231_fu_86_p4 = grp_fu_80_p2.read().range(48, 26);
}

void invCosh_1::thread_tmp_233_fu_116_p4() {
    tmp_233_fu_116_p4 = r_V_fu_100_p2.read().range(24, 23);
}

void invCosh_1::thread_tmp_cast_cast_fu_96_p1() {
    tmp_cast_cast_fu_96_p1 = esl_zext<25,23>(tmp_231_fu_86_p4.read());
}

void invCosh_1::thread_tmp_fu_64_p1() {
    tmp_fu_64_p1 = data_V_read.read().range(11-1, 0);
}

void invCosh_1::thread_ap_NS_fsm() {
    switch (ap_CS_fsm.read().to_uint64()) {
        case 1 : 
            ap_NS_fsm = ap_ST_fsm_pp0_stage0;
break;
        default : 
            ap_NS_fsm = "X";
            break;
    }
}

}
