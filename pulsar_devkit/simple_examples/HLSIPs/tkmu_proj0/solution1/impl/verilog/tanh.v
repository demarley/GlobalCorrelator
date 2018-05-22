// ==============================================================
// RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
// Version: 2016.4
// Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
// 
// ===========================================================

`timescale 1 ns / 1 ps 

module tanh (
        ap_clk,
        ap_rst,
        ap_start,
        ap_done,
        ap_idle,
        ap_ready,
        data_V_read,
        ap_return
);

parameter    ap_ST_fsm_pp0_stage0 = 1'b1;
parameter    ap_const_lv32_0 = 32'b00000000000000000000000000000000;
parameter    ap_const_lv13_0 = 13'b0000000000000;
parameter    ap_const_lv55_AAAAAAB = 55'b1010101010101010101010101011;
parameter    ap_const_lv32_C = 32'b1100;
parameter    ap_const_lv32_1D = 32'b11101;
parameter    ap_const_lv32_36 = 32'b110110;
parameter    ap_const_lv55_0 = 55'b0000000000000000000000000000000000000000000000000000000;
parameter    ap_const_lv26_0 = 26'b00000000000000000000000000;
parameter    ap_const_lv26_800000 = 26'b100000000000000000000000;
parameter    ap_const_lv32_A = 32'b1010;
parameter    ap_const_lv32_19 = 32'b11001;
parameter    ap_const_lv10_0 = 10'b0000000000;
parameter    ap_const_lv17_1 = 17'b1;
parameter    ap_const_lv32_1F = 32'b11111;
parameter    ap_const_lv32_D = 32'b1101;
parameter    ap_const_lv19_0 = 19'b0000000000000000000;
parameter    ap_const_lv10_3FA = 10'b1111111010;

input   ap_clk;
input   ap_rst;
input   ap_start;
output   ap_done;
output   ap_idle;
output   ap_ready;
input  [12:0] data_V_read;
output  [9:0] ap_return;

reg ap_done;
reg ap_idle;
reg ap_ready;

(* fsm_encoding = "none" *) reg   [0:0] ap_CS_fsm;
wire   [0:0] ap_CS_fsm_pp0_stage0;
wire    ap_enable_reg_pp0_iter0;
reg    ap_enable_reg_pp0_iter1;
reg    ap_enable_reg_pp0_iter2;
reg    ap_enable_reg_pp0_iter3;
reg    ap_enable_reg_pp0_iter4;
reg    ap_enable_reg_pp0_iter5;
reg    ap_enable_reg_pp0_iter6;
reg    ap_enable_reg_pp0_iter7;
reg    ap_enable_reg_pp0_iter8;
reg    ap_enable_reg_pp0_iter9;
reg    ap_enable_reg_pp0_iter10;
wire   [12:0] tanh_table2_address0;
reg    tanh_table2_ce0;
wire   [9:0] tanh_table2_q0;
reg   [0:0] tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter1_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter2_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter3_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter4_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter5_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter6_tmp_224_reg_283;
reg   [0:0] ap_pipeline_reg_pp0_iter7_tmp_224_reg_283;
wire   [54:0] grp_fu_96_p2;
reg   [54:0] mul_reg_289;
reg   [25:0] tmp_226_reg_294;
reg   [25:0] ap_pipeline_reg_pp0_iter7_tmp_226_reg_294;
wire   [25:0] neg_ti_fu_141_p2;
reg   [25:0] neg_ti_reg_300;
wire   [31:0] index_fu_212_p3;
reg   [31:0] index_reg_305;
reg   [0:0] tmp_229_reg_310;
reg   [0:0] ap_pipeline_reg_pp0_iter9_tmp_229_reg_310;
reg   [18:0] tmp_230_reg_316;
wire   [0:0] icmp_fu_238_p2;
reg   [0:0] icmp_reg_321;
wire   [63:0] tmp_192_fu_243_p1;
wire  signed [25:0] tmp_fu_84_p3;
wire   [28:0] grp_fu_96_p0;
wire   [54:0] neg_mul_fu_120_p2;
wire   [25:0] tmp_225_fu_125_p4;
wire   [25:0] p_v_v_fu_135_p3;
wire   [25:0] tmp_194_fu_147_p3;
wire   [25:0] r_V_fu_152_p2;
wire   [15:0] tmp_193_fu_158_p4;
wire   [9:0] tmp_228_fu_180_p1;
wire  signed [16:0] ret_V_cast_fu_168_p1;
wire   [16:0] ret_V_fu_190_p2;
wire   [0:0] tmp_s_fu_184_p2;
wire  signed [31:0] tmp_195_fu_196_p1;
wire  signed [31:0] tmp_196_fu_200_p1;
wire   [0:0] tmp_227_fu_172_p3;
wire   [31:0] tmp_198_fu_204_p3;
wire   [0:0] sel_tmp1_fu_247_p2;
wire   [0:0] sel_tmp2_fu_252_p2;
wire   [0:0] tmp_197_fu_265_p2;
wire   [9:0] sel_tmp_cast_fu_257_p3;
reg   [0:0] ap_NS_fsm;
reg    ap_pipeline_idle_pp0;

// power-on initialization
initial begin
#0 ap_CS_fsm = 1'b1;
#0 ap_enable_reg_pp0_iter1 = 1'b0;
#0 ap_enable_reg_pp0_iter2 = 1'b0;
#0 ap_enable_reg_pp0_iter3 = 1'b0;
#0 ap_enable_reg_pp0_iter4 = 1'b0;
#0 ap_enable_reg_pp0_iter5 = 1'b0;
#0 ap_enable_reg_pp0_iter6 = 1'b0;
#0 ap_enable_reg_pp0_iter7 = 1'b0;
#0 ap_enable_reg_pp0_iter8 = 1'b0;
#0 ap_enable_reg_pp0_iter9 = 1'b0;
#0 ap_enable_reg_pp0_iter10 = 1'b0;
end

tanh_tanh_table2 #(
    .DataWidth( 10 ),
    .AddressRange( 8192 ),
    .AddressWidth( 13 ))
tanh_table2_U(
    .clk(ap_clk),
    .reset(ap_rst),
    .address0(tanh_table2_address0),
    .ce0(tanh_table2_ce0),
    .q0(tanh_table2_q0)
);

tkmu_simple_hw_muhbi #(
    .ID( 1 ),
    .NUM_STAGE( 7 ),
    .din0_WIDTH( 29 ),
    .din1_WIDTH( 26 ),
    .dout_WIDTH( 55 ))
tkmu_simple_hw_muhbi_U13(
    .clk(ap_clk),
    .reset(ap_rst),
    .din0(grp_fu_96_p0),
    .din1(tmp_fu_84_p3),
    .ce(1'b1),
    .dout(grp_fu_96_p2)
);

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_CS_fsm <= ap_ST_fsm_pp0_stage0;
    end else begin
        ap_CS_fsm <= ap_NS_fsm;
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter1 <= 1'b0;
    end else begin
        if (((ap_CS_fsm_pp0_stage0 == 1'b1) & ~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)))) begin
            ap_enable_reg_pp0_iter1 <= ap_start;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter10 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter10 <= ap_enable_reg_pp0_iter9;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter2 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter2 <= ap_enable_reg_pp0_iter1;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter3 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter3 <= ap_enable_reg_pp0_iter2;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter4 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter4 <= ap_enable_reg_pp0_iter3;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter5 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter5 <= ap_enable_reg_pp0_iter4;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter6 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter6 <= ap_enable_reg_pp0_iter5;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter7 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter7 <= ap_enable_reg_pp0_iter6;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter8 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter8 <= ap_enable_reg_pp0_iter7;
        end
    end
end

always @ (posedge ap_clk) begin
    if (ap_rst == 1'b1) begin
        ap_enable_reg_pp0_iter9 <= 1'b0;
    end else begin
        if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
            ap_enable_reg_pp0_iter9 <= ap_enable_reg_pp0_iter8;
        end
    end
end

always @ (posedge ap_clk) begin
    if (((ap_CS_fsm_pp0_stage0 == 1'b1) & ~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)))) begin
        ap_pipeline_reg_pp0_iter1_tmp_224_reg_283 <= tmp_224_reg_283;
        tmp_224_reg_283 <= data_V_read[ap_const_lv32_C];
    end
end

always @ (posedge ap_clk) begin
    if (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0))) begin
        ap_pipeline_reg_pp0_iter2_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter1_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter3_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter2_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter4_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter3_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter5_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter4_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter6_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter5_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter7_tmp_224_reg_283 <= ap_pipeline_reg_pp0_iter6_tmp_224_reg_283;
        ap_pipeline_reg_pp0_iter7_tmp_226_reg_294 <= tmp_226_reg_294;
        ap_pipeline_reg_pp0_iter9_tmp_229_reg_310 <= tmp_229_reg_310;
        icmp_reg_321 <= icmp_fu_238_p2;
        index_reg_305 <= index_fu_212_p3;
        mul_reg_289 <= grp_fu_96_p2;
        tmp_226_reg_294 <= {{grp_fu_96_p2[ap_const_lv32_36 : ap_const_lv32_1D]}};
        tmp_229_reg_310 <= index_fu_212_p3[ap_const_lv32_1F];
        tmp_230_reg_316 <= {{index_fu_212_p3[ap_const_lv32_1F : ap_const_lv32_D]}};
    end
end

always @ (posedge ap_clk) begin
    if ((~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)) & ~(ap_pipeline_reg_pp0_iter6_tmp_224_reg_283 == 1'b0))) begin
        neg_ti_reg_300 <= neg_ti_fu_141_p2;
    end
end

always @ (*) begin
    if ((((1'b0 == ap_start) & (ap_CS_fsm_pp0_stage0 == 1'b1) & (1'b1 == ap_enable_reg_pp0_iter0)) | (~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)) & (1'b1 == ap_enable_reg_pp0_iter10)))) begin
        ap_done = 1'b1;
    end else begin
        ap_done = 1'b0;
    end
end

always @ (*) begin
    if (((1'b0 == ap_start) & (ap_CS_fsm_pp0_stage0 == 1'b1) & (1'b0 == ap_enable_reg_pp0_iter0) & (1'b0 == ap_enable_reg_pp0_iter1) & (1'b0 == ap_enable_reg_pp0_iter2) & (1'b0 == ap_enable_reg_pp0_iter3) & (1'b0 == ap_enable_reg_pp0_iter4) & (1'b0 == ap_enable_reg_pp0_iter5) & (1'b0 == ap_enable_reg_pp0_iter6) & (1'b0 == ap_enable_reg_pp0_iter7) & (1'b0 == ap_enable_reg_pp0_iter8) & (1'b0 == ap_enable_reg_pp0_iter9) & (1'b0 == ap_enable_reg_pp0_iter10))) begin
        ap_idle = 1'b1;
    end else begin
        ap_idle = 1'b0;
    end
end

always @ (*) begin
    if (((1'b0 == ap_start) & (1'b0 == ap_enable_reg_pp0_iter0) & (1'b0 == ap_enable_reg_pp0_iter1) & (1'b0 == ap_enable_reg_pp0_iter2) & (1'b0 == ap_enable_reg_pp0_iter3) & (1'b0 == ap_enable_reg_pp0_iter4) & (1'b0 == ap_enable_reg_pp0_iter5) & (1'b0 == ap_enable_reg_pp0_iter6) & (1'b0 == ap_enable_reg_pp0_iter7) & (1'b0 == ap_enable_reg_pp0_iter8) & (1'b0 == ap_enable_reg_pp0_iter9))) begin
        ap_pipeline_idle_pp0 = 1'b1;
    end else begin
        ap_pipeline_idle_pp0 = 1'b0;
    end
end

always @ (*) begin
    if (((ap_CS_fsm_pp0_stage0 == 1'b1) & (1'b1 == ap_enable_reg_pp0_iter0) & ~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)))) begin
        ap_ready = 1'b1;
    end else begin
        ap_ready = 1'b0;
    end
end

always @ (*) begin
    if ((~((1'b1 == ap_enable_reg_pp0_iter0) & (ap_start == 1'b0)) & (1'b1 == ap_enable_reg_pp0_iter9))) begin
        tanh_table2_ce0 = 1'b1;
    end else begin
        tanh_table2_ce0 = 1'b0;
    end
end

always @ (*) begin
    case (ap_CS_fsm)
        ap_ST_fsm_pp0_stage0 : begin
            ap_NS_fsm = ap_ST_fsm_pp0_stage0;
        end
        default : begin
            ap_NS_fsm = 'bx;
        end
    endcase
end

assign ap_CS_fsm_pp0_stage0 = ap_CS_fsm[ap_const_lv32_0];

assign ap_enable_reg_pp0_iter0 = ap_start;

assign ap_return = ((tmp_197_fu_265_p2[0:0] === 1'b1) ? sel_tmp_cast_fu_257_p3 : tanh_table2_q0);

assign grp_fu_96_p0 = ap_const_lv55_AAAAAAB;

assign icmp_fu_238_p2 = (($signed(tmp_230_reg_316) > $signed(19'b0000000000000000000)) ? 1'b1 : 1'b0);

assign index_fu_212_p3 = ((tmp_227_fu_172_p3[0:0] === 1'b1) ? tmp_198_fu_204_p3 : tmp_195_fu_196_p1);

assign neg_mul_fu_120_p2 = (ap_const_lv55_0 - mul_reg_289);

assign neg_ti_fu_141_p2 = (ap_const_lv26_0 - p_v_v_fu_135_p3);

assign p_v_v_fu_135_p3 = ((ap_pipeline_reg_pp0_iter6_tmp_224_reg_283[0:0] === 1'b1) ? tmp_225_fu_125_p4 : tmp_226_reg_294);

assign r_V_fu_152_p2 = (ap_const_lv26_800000 - tmp_194_fu_147_p3);

assign ret_V_cast_fu_168_p1 = $signed(tmp_193_fu_158_p4);

assign ret_V_fu_190_p2 = ($signed(ap_const_lv17_1) + $signed(ret_V_cast_fu_168_p1));

assign sel_tmp1_fu_247_p2 = (ap_pipeline_reg_pp0_iter9_tmp_229_reg_310 ^ 1'b1);

assign sel_tmp2_fu_252_p2 = (icmp_reg_321 & sel_tmp1_fu_247_p2);

assign sel_tmp_cast_fu_257_p3 = ((sel_tmp2_fu_252_p2[0:0] === 1'b1) ? ap_const_lv10_0 : ap_const_lv10_3FA);

assign tanh_table2_address0 = tmp_192_fu_243_p1;

assign tmp_192_fu_243_p1 = index_reg_305;

assign tmp_193_fu_158_p4 = {{r_V_fu_152_p2[ap_const_lv32_19 : ap_const_lv32_A]}};

assign tmp_194_fu_147_p3 = ((ap_pipeline_reg_pp0_iter7_tmp_224_reg_283[0:0] === 1'b1) ? neg_ti_reg_300 : ap_pipeline_reg_pp0_iter7_tmp_226_reg_294);

assign tmp_195_fu_196_p1 = $signed(tmp_193_fu_158_p4);

assign tmp_196_fu_200_p1 = $signed(ret_V_fu_190_p2);

assign tmp_197_fu_265_p2 = (sel_tmp2_fu_252_p2 | ap_pipeline_reg_pp0_iter9_tmp_229_reg_310);

assign tmp_198_fu_204_p3 = ((tmp_s_fu_184_p2[0:0] === 1'b1) ? tmp_195_fu_196_p1 : tmp_196_fu_200_p1);

assign tmp_225_fu_125_p4 = {{neg_mul_fu_120_p2[ap_const_lv32_36 : ap_const_lv32_1D]}};

assign tmp_227_fu_172_p3 = r_V_fu_152_p2[ap_const_lv32_19];

assign tmp_228_fu_180_p1 = r_V_fu_152_p2[9:0];

assign tmp_fu_84_p3 = {{data_V_read}, {ap_const_lv13_0}};

assign tmp_s_fu_184_p2 = ((tmp_228_fu_180_p1 == ap_const_lv10_0) ? 1'b1 : 1'b0);

endmodule //tanh