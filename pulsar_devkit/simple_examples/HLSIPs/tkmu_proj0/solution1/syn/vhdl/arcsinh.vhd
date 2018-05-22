-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
-- Version: 2016.4
-- Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity arcsinh is
port (
    ap_clk : IN STD_LOGIC;
    ap_rst : IN STD_LOGIC;
    ap_start : IN STD_LOGIC;
    ap_done : OUT STD_LOGIC;
    ap_idle : OUT STD_LOGIC;
    ap_ready : OUT STD_LOGIC;
    data_V_read : IN STD_LOGIC_VECTOR (13 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (11 downto 0) );
end;


architecture behav of arcsinh is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_ST_fsm_pp0_stage0 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv32_0 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000000000";
    constant ap_const_lv1_1 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv1_0 : STD_LOGIC_VECTOR (0 downto 0) := "0";
    constant ap_const_lv13_0 : STD_LOGIC_VECTOR (12 downto 0) := "0000000000000";
    constant ap_const_lv55_AAAAAAB : STD_LOGIC_VECTOR (54 downto 0) := "0000000000000000000000000001010101010101010101010101011";
    constant ap_const_lv32_D : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000001101";
    constant ap_const_lv55_0 : STD_LOGIC_VECTOR (54 downto 0) := "0000000000000000000000000000000000000000000000000000000";
    constant ap_const_lv32_1E : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011110";
    constant ap_const_lv32_36 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000110110";
    constant ap_const_lv26_0 : STD_LOGIC_VECTOR (25 downto 0) := "00000000000000000000000000";
    constant ap_const_lv26_800000 : STD_LOGIC_VECTOR (25 downto 0) := "00100000000000000000000000";
    constant ap_const_lv32_A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000001010";
    constant ap_const_lv32_19 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011001";
    constant ap_const_lv10_0 : STD_LOGIC_VECTOR (9 downto 0) := "0000000000";
    constant ap_const_lv17_1 : STD_LOGIC_VECTOR (16 downto 0) := "00000000000000001";
    constant ap_const_lv32_1F : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011111";
    constant ap_const_lv19_0 : STD_LOGIC_VECTOR (18 downto 0) := "0000000000000000000";
    constant ap_const_lv12_0 : STD_LOGIC_VECTOR (11 downto 0) := "000000000000";
    constant ap_const_lv12_9F7 : STD_LOGIC_VECTOR (11 downto 0) := "100111110111";

    signal ap_CS_fsm : STD_LOGIC_VECTOR (0 downto 0) := "1";
    attribute fsm_encoding : string;
    attribute fsm_encoding of ap_CS_fsm : signal is "none";
    signal ap_CS_fsm_pp0_stage0 : STD_LOGIC_VECTOR (0 downto 0);
    attribute fsm_encoding of ap_CS_fsm_pp0_stage0 : signal is "none";
    signal ap_enable_reg_pp0_iter0 : STD_LOGIC;
    signal ap_enable_reg_pp0_iter1 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter2 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter3 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter4 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter5 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter6 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter7 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter8 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter9 : STD_LOGIC := '0';
    signal ap_enable_reg_pp0_iter10 : STD_LOGIC := '0';
    signal arcsinh_table9_address0 : STD_LOGIC_VECTOR (12 downto 0);
    signal arcsinh_table9_ce0 : STD_LOGIC;
    signal arcsinh_table9_q0 : STD_LOGIC_VECTOR (11 downto 0);
    signal tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter1_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter2_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter3_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter4_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter5_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter6_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter7_tmp_reg_292 : STD_LOGIC_VECTOR (0 downto 0);
    signal p_v_v_fu_136_p3 : STD_LOGIC_VECTOR (24 downto 0);
    signal p_v_v_reg_298 : STD_LOGIC_VECTOR (24 downto 0);
    signal ap_pipeline_reg_pp0_iter7_p_v_v_reg_298 : STD_LOGIC_VECTOR (24 downto 0);
    signal neg_ti_fu_146_p2 : STD_LOGIC_VECTOR (25 downto 0);
    signal neg_ti_reg_304 : STD_LOGIC_VECTOR (25 downto 0);
    signal index_fu_221_p3 : STD_LOGIC_VECTOR (31 downto 0);
    signal index_reg_309 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_283_reg_314 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter9_tmp_283_reg_314 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_284_reg_320 : STD_LOGIC_VECTOR (18 downto 0);
    signal icmp_fu_247_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal icmp_reg_325 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_211_fu_252_p1 : STD_LOGIC_VECTOR (63 downto 0);
    signal r_V_fu_84_p3 : STD_LOGIC_VECTOR (26 downto 0);
    signal grp_fu_96_p0 : STD_LOGIC_VECTOR (28 downto 0);
    signal grp_fu_96_p2 : STD_LOGIC_VECTOR (54 downto 0);
    signal neg_mul_fu_110_p2 : STD_LOGIC_VECTOR (54 downto 0);
    signal tmp_279_fu_116_p4 : STD_LOGIC_VECTOR (24 downto 0);
    signal tmp_280_fu_126_p4 : STD_LOGIC_VECTOR (24 downto 0);
    signal trunc_fu_143_p1 : STD_LOGIC_VECTOR (25 downto 0);
    signal tmp_224_fu_152_p1 : STD_LOGIC_VECTOR (25 downto 0);
    signal tmp_225_fu_155_p3 : STD_LOGIC_VECTOR (25 downto 0);
    signal r_V_21_fu_161_p2 : STD_LOGIC_VECTOR (25 downto 0);
    signal tmp_213_fu_167_p4 : STD_LOGIC_VECTOR (15 downto 0);
    signal tmp_282_fu_189_p1 : STD_LOGIC_VECTOR (9 downto 0);
    signal ret_V_cast_fu_177_p1 : STD_LOGIC_VECTOR (16 downto 0);
    signal ret_V_fu_199_p2 : STD_LOGIC_VECTOR (16 downto 0);
    signal tmp_s_fu_193_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_226_fu_205_p1 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_227_fu_209_p1 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_281_fu_181_p3 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_228_fu_213_p3 : STD_LOGIC_VECTOR (31 downto 0);
    signal sel_tmp1_fu_256_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal sel_tmp2_fu_261_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_215_fu_274_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal sel_tmp_fu_266_p3 : STD_LOGIC_VECTOR (11 downto 0);
    signal ap_NS_fsm : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_idle_pp0 : STD_LOGIC;

    component tkmu_simple_hw_mucud IS
    generic (
        ID : INTEGER;
        NUM_STAGE : INTEGER;
        din0_WIDTH : INTEGER;
        din1_WIDTH : INTEGER;
        dout_WIDTH : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        din0 : IN STD_LOGIC_VECTOR (28 downto 0);
        din1 : IN STD_LOGIC_VECTOR (26 downto 0);
        ce : IN STD_LOGIC;
        dout : OUT STD_LOGIC_VECTOR (54 downto 0) );
    end component;


    component arcsinh_arcsinh_tbkb IS
    generic (
        DataWidth : INTEGER;
        AddressRange : INTEGER;
        AddressWidth : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        address0 : IN STD_LOGIC_VECTOR (12 downto 0);
        ce0 : IN STD_LOGIC;
        q0 : OUT STD_LOGIC_VECTOR (11 downto 0) );
    end component;



begin
    arcsinh_table9_U : component arcsinh_arcsinh_tbkb
    generic map (
        DataWidth => 12,
        AddressRange => 8192,
        AddressWidth => 13)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        address0 => arcsinh_table9_address0,
        ce0 => arcsinh_table9_ce0,
        q0 => arcsinh_table9_q0);

    tkmu_simple_hw_mucud_U1 : component tkmu_simple_hw_mucud
    generic map (
        ID => 1,
        NUM_STAGE => 7,
        din0_WIDTH => 29,
        din1_WIDTH => 27,
        dout_WIDTH => 55)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        din0 => grp_fu_96_p0,
        din1 => r_V_fu_84_p3,
        ce => ap_const_logic_1,
        dout => grp_fu_96_p2);





    ap_CS_fsm_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_CS_fsm <= ap_ST_fsm_pp0_stage0;
            else
                ap_CS_fsm <= ap_NS_fsm;
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter1_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter1 <= ap_const_logic_0;
            else
                if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))))) then 
                    ap_enable_reg_pp0_iter1 <= ap_start;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter10_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter10 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter10 <= ap_enable_reg_pp0_iter9;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter2_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter2 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter2 <= ap_enable_reg_pp0_iter1;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter3_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter3 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter3 <= ap_enable_reg_pp0_iter2;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter4_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter4 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter4 <= ap_enable_reg_pp0_iter3;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter5_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter5 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter5 <= ap_enable_reg_pp0_iter4;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter6_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter6 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter6 <= ap_enable_reg_pp0_iter5;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter7_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter7 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter7 <= ap_enable_reg_pp0_iter6;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter8_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter8 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter8 <= ap_enable_reg_pp0_iter7;
                end if; 
            end if;
        end if;
    end process;


    ap_enable_reg_pp0_iter9_assign_proc : process(ap_clk)
    begin
        if (ap_clk'event and ap_clk =  '1') then
            if (ap_rst = '1') then
                ap_enable_reg_pp0_iter9 <= ap_const_logic_0;
            else
                if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then 
                    ap_enable_reg_pp0_iter9 <= ap_enable_reg_pp0_iter8;
                end if; 
            end if;
        end if;
    end process;

    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))))) then
                ap_pipeline_reg_pp0_iter1_tmp_reg_292 <= tmp_reg_292;
                tmp_reg_292 <= data_V_read(13 downto 13);
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then
                ap_pipeline_reg_pp0_iter2_tmp_reg_292 <= ap_pipeline_reg_pp0_iter1_tmp_reg_292;
                ap_pipeline_reg_pp0_iter3_tmp_reg_292 <= ap_pipeline_reg_pp0_iter2_tmp_reg_292;
                ap_pipeline_reg_pp0_iter4_tmp_reg_292 <= ap_pipeline_reg_pp0_iter3_tmp_reg_292;
                ap_pipeline_reg_pp0_iter5_tmp_reg_292 <= ap_pipeline_reg_pp0_iter4_tmp_reg_292;
                ap_pipeline_reg_pp0_iter6_tmp_reg_292 <= ap_pipeline_reg_pp0_iter5_tmp_reg_292;
                ap_pipeline_reg_pp0_iter7_p_v_v_reg_298 <= p_v_v_reg_298;
                ap_pipeline_reg_pp0_iter7_tmp_reg_292 <= ap_pipeline_reg_pp0_iter6_tmp_reg_292;
                ap_pipeline_reg_pp0_iter9_tmp_283_reg_314 <= tmp_283_reg_314;
                icmp_reg_325 <= icmp_fu_247_p2;
                index_reg_309 <= index_fu_221_p3;
                p_v_v_reg_298 <= p_v_v_fu_136_p3;
                tmp_283_reg_314 <= index_fu_221_p3(31 downto 31);
                tmp_284_reg_320 <= index_fu_221_p3(31 downto 13);
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and not((ap_pipeline_reg_pp0_iter6_tmp_reg_292 = ap_const_lv1_0)))) then
                neg_ti_reg_304 <= neg_ti_fu_146_p2;
            end if;
        end if;
    end process;

    ap_NS_fsm_assign_proc : process (ap_start, ap_CS_fsm, ap_enable_reg_pp0_iter0, ap_pipeline_idle_pp0)
    begin
        case ap_CS_fsm is
            when ap_ST_fsm_pp0_stage0 => 
                ap_NS_fsm <= ap_ST_fsm_pp0_stage0;
            when others =>  
                ap_NS_fsm <= "X";
        end case;
    end process;
    ap_CS_fsm_pp0_stage0 <= ap_CS_fsm(0 downto 0);

    ap_done_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter10)
    begin
        if ((((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_1 = ap_enable_reg_pp0_iter0)) or (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_enable_reg_pp0_iter10)))) then 
            ap_done <= ap_const_logic_1;
        else 
            ap_done <= ap_const_logic_0;
        end if; 
    end process;

    ap_enable_reg_pp0_iter0 <= ap_start;

    ap_idle_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5, ap_enable_reg_pp0_iter6, ap_enable_reg_pp0_iter7, ap_enable_reg_pp0_iter8, ap_enable_reg_pp0_iter9, ap_enable_reg_pp0_iter10)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5) and (ap_const_logic_0 = ap_enable_reg_pp0_iter6) and (ap_const_logic_0 = ap_enable_reg_pp0_iter7) and (ap_const_logic_0 = ap_enable_reg_pp0_iter8) and (ap_const_logic_0 = ap_enable_reg_pp0_iter9) and (ap_const_logic_0 = ap_enable_reg_pp0_iter10))) then 
            ap_idle <= ap_const_logic_1;
        else 
            ap_idle <= ap_const_logic_0;
        end if; 
    end process;


    ap_pipeline_idle_pp0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5, ap_enable_reg_pp0_iter6, ap_enable_reg_pp0_iter7, ap_enable_reg_pp0_iter8, ap_enable_reg_pp0_iter9)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5) and (ap_const_logic_0 = ap_enable_reg_pp0_iter6) and (ap_const_logic_0 = ap_enable_reg_pp0_iter7) and (ap_const_logic_0 = ap_enable_reg_pp0_iter8) and (ap_const_logic_0 = ap_enable_reg_pp0_iter9))) then 
            ap_pipeline_idle_pp0 <= ap_const_logic_1;
        else 
            ap_pipeline_idle_pp0 <= ap_const_logic_0;
        end if; 
    end process;


    ap_ready_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0)
    begin
        if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_1 = ap_enable_reg_pp0_iter0) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))))) then 
            ap_ready <= ap_const_logic_1;
        else 
            ap_ready <= ap_const_logic_0;
        end if; 
    end process;

    ap_return <= 
        sel_tmp_fu_266_p3 when (tmp_215_fu_274_p2(0) = '1') else 
        arcsinh_table9_q0;
    arcsinh_table9_address0 <= tmp_211_fu_252_p1(13 - 1 downto 0);

    arcsinh_table9_ce0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter9)
    begin
        if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_enable_reg_pp0_iter9))) then 
            arcsinh_table9_ce0 <= ap_const_logic_1;
        else 
            arcsinh_table9_ce0 <= ap_const_logic_0;
        end if; 
    end process;

    grp_fu_96_p0 <= ap_const_lv55_AAAAAAB(29 - 1 downto 0);
    icmp_fu_247_p2 <= "1" when (signed(tmp_284_reg_320) > signed(ap_const_lv19_0)) else "0";
    index_fu_221_p3 <= 
        tmp_228_fu_213_p3 when (tmp_281_fu_181_p3(0) = '1') else 
        tmp_226_fu_205_p1;
    neg_mul_fu_110_p2 <= std_logic_vector(unsigned(ap_const_lv55_0) - unsigned(grp_fu_96_p2));
    neg_ti_fu_146_p2 <= std_logic_vector(unsigned(ap_const_lv26_0) - unsigned(trunc_fu_143_p1));
    p_v_v_fu_136_p3 <= 
        tmp_279_fu_116_p4 when (ap_pipeline_reg_pp0_iter5_tmp_reg_292(0) = '1') else 
        tmp_280_fu_126_p4;
    r_V_21_fu_161_p2 <= std_logic_vector(unsigned(ap_const_lv26_800000) - unsigned(tmp_225_fu_155_p3));
    r_V_fu_84_p3 <= (data_V_read & ap_const_lv13_0);
        ret_V_cast_fu_177_p1 <= std_logic_vector(resize(signed(tmp_213_fu_167_p4),17));

    ret_V_fu_199_p2 <= std_logic_vector(unsigned(ap_const_lv17_1) + unsigned(ret_V_cast_fu_177_p1));
    sel_tmp1_fu_256_p2 <= (ap_pipeline_reg_pp0_iter9_tmp_283_reg_314 xor ap_const_lv1_1);
    sel_tmp2_fu_261_p2 <= (icmp_reg_325 and sel_tmp1_fu_256_p2);
    sel_tmp_fu_266_p3 <= 
        ap_const_lv12_0 when (sel_tmp2_fu_261_p2(0) = '1') else 
        ap_const_lv12_9F7;
    tmp_211_fu_252_p1 <= std_logic_vector(resize(unsigned(index_reg_309),64));
    tmp_213_fu_167_p4 <= r_V_21_fu_161_p2(25 downto 10);
    tmp_215_fu_274_p2 <= (sel_tmp2_fu_261_p2 or ap_pipeline_reg_pp0_iter9_tmp_283_reg_314);
        tmp_224_fu_152_p1 <= std_logic_vector(resize(signed(ap_pipeline_reg_pp0_iter7_p_v_v_reg_298),26));

    tmp_225_fu_155_p3 <= 
        neg_ti_reg_304 when (ap_pipeline_reg_pp0_iter7_tmp_reg_292(0) = '1') else 
        tmp_224_fu_152_p1;
        tmp_226_fu_205_p1 <= std_logic_vector(resize(signed(tmp_213_fu_167_p4),32));

        tmp_227_fu_209_p1 <= std_logic_vector(resize(signed(ret_V_fu_199_p2),32));

    tmp_228_fu_213_p3 <= 
        tmp_226_fu_205_p1 when (tmp_s_fu_193_p2(0) = '1') else 
        tmp_227_fu_209_p1;
    tmp_279_fu_116_p4 <= neg_mul_fu_110_p2(54 downto 30);
    tmp_280_fu_126_p4 <= grp_fu_96_p2(54 downto 30);
    tmp_281_fu_181_p3 <= r_V_21_fu_161_p2(25 downto 25);
    tmp_282_fu_189_p1 <= r_V_21_fu_161_p2(10 - 1 downto 0);
    tmp_s_fu_193_p2 <= "1" when (tmp_282_fu_189_p1 = ap_const_lv10_0) else "0";
        trunc_fu_143_p1 <= std_logic_vector(resize(signed(p_v_v_reg_298),26));

end behav;
