-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
-- Version: 2016.4
-- Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity delta_plus_LUT is
port (
    ap_clk : IN STD_LOGIC;
    ap_rst : IN STD_LOGIC;
    ap_start : IN STD_LOGIC;
    ap_done : OUT STD_LOGIC;
    ap_idle : OUT STD_LOGIC;
    ap_ready : OUT STD_LOGIC;
    ap_ce : IN STD_LOGIC;
    data_V_read : IN STD_LOGIC_VECTOR (10 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (4 downto 0) );
end;


architecture behav of delta_plus_LUT is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_ST_fsm_pp0_stage0 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv32_0 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000000000";
    constant ap_const_lv1_1 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv1_0 : STD_LOGIC_VECTOR (0 downto 0) := "0";
    constant ap_const_lv10_0 : STD_LOGIC_VECTOR (9 downto 0) := "0000000000";
    constant ap_const_lv44_222223 : STD_LOGIC_VECTOR (43 downto 0) := "00000000000000000000001000100010001000100011";
    constant ap_const_lv32_A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000001010";
    constant ap_const_lv43_0 : STD_LOGIC_VECTOR (42 downto 0) := "0000000000000000000000000000000000000000000";
    constant ap_const_lv32_19 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011001";
    constant ap_const_lv32_2A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000101010";
    constant ap_const_lv32_2B : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000101011";
    constant ap_const_lv19_0 : STD_LOGIC_VECTOR (18 downto 0) := "0000000000000000000";
    constant ap_const_lv19_10000 : STD_LOGIC_VECTOR (18 downto 0) := "0010000000000000000";
    constant ap_const_lv32_6 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000000110";
    constant ap_const_lv32_12 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000010010";
    constant ap_const_lv6_0 : STD_LOGIC_VECTOR (5 downto 0) := "000000";
    constant ap_const_lv14_1 : STD_LOGIC_VECTOR (13 downto 0) := "00000000000001";
    constant ap_const_lv32_1F : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011111";
    constant ap_const_lv22_0 : STD_LOGIC_VECTOR (21 downto 0) := "0000000000000000000000";
    constant ap_const_lv5_0 : STD_LOGIC_VECTOR (4 downto 0) := "00000";
    constant ap_const_lv5_11 : STD_LOGIC_VECTOR (4 downto 0) := "10001";

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
    signal delta_plus_table4_address0 : STD_LOGIC_VECTOR (9 downto 0);
    signal delta_plus_table4_ce0 : STD_LOGIC;
    signal delta_plus_table4_q0 : STD_LOGIC_VECTOR (4 downto 0);
    signal tmp_250_reg_307 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter1_tmp_250_reg_307 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter2_tmp_250_reg_307 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter3_tmp_250_reg_307 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter4_tmp_250_reg_307 : STD_LOGIC_VECTOR (0 downto 0);
    signal p_v_fu_152_p3 : STD_LOGIC_VECTOR (21 downto 0);
    signal p_v_reg_313 : STD_LOGIC_VECTOR (21 downto 0);
    signal tmp_253_fu_159_p1 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_253_reg_318 : STD_LOGIC_VECTOR (18 downto 0);
    signal r_V_15_fu_178_p2 : STD_LOGIC_VECTOR (18 downto 0);
    signal r_V_15_reg_323 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_204_fu_184_p4 : STD_LOGIC_VECTOR (12 downto 0);
    signal tmp_204_reg_328 : STD_LOGIC_VECTOR (12 downto 0);
    signal tmp_s_fu_202_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_s_reg_333 : STD_LOGIC_VECTOR (0 downto 0);
    signal ret_V_fu_208_p2 : STD_LOGIC_VECTOR (13 downto 0);
    signal ret_V_reg_338 : STD_LOGIC_VECTOR (13 downto 0);
    signal tmp_257_reg_343 : STD_LOGIC_VECTOR (0 downto 0);
    signal icmp_fu_260_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal icmp_reg_349 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_202_fu_266_p1 : STD_LOGIC_VECTOR (63 downto 0);
    signal r_V_fu_88_p3 : STD_LOGIC_VECTOR (20 downto 0);
    signal grp_fu_100_p0 : STD_LOGIC_VECTOR (22 downto 0);
    signal grp_fu_100_p2 : STD_LOGIC_VECTOR (43 downto 0);
    signal tmp_249_fu_114_p1 : STD_LOGIC_VECTOR (42 downto 0);
    signal neg_mul_fu_118_p2 : STD_LOGIC_VECTOR (42 downto 0);
    signal tmp_251_fu_124_p4 : STD_LOGIC_VECTOR (17 downto 0);
    signal tmp_252_fu_138_p4 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_fu_134_p1 : STD_LOGIC_VECTOR (21 downto 0);
    signal tmp_209_fu_148_p1 : STD_LOGIC_VECTOR (21 downto 0);
    signal neg_ti_fu_163_p2 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_254_fu_168_p1 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_210_fu_171_p3 : STD_LOGIC_VECTOR (18 downto 0);
    signal tmp_256_fu_198_p1 : STD_LOGIC_VECTOR (5 downto 0);
    signal ret_V_cast_fu_194_p1 : STD_LOGIC_VECTOR (13 downto 0);
    signal tmp_211_fu_221_p1 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_212_fu_224_p1 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_255_fu_214_p3 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_213_fu_227_p3 : STD_LOGIC_VECTOR (31 downto 0);
    signal index_fu_234_p3 : STD_LOGIC_VECTOR (31 downto 0);
    signal tmp_258_fu_250_p4 : STD_LOGIC_VECTOR (21 downto 0);
    signal sel_tmp1_fu_271_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal sel_tmp2_fu_276_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_206_fu_289_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal sel_tmp_fu_281_p3 : STD_LOGIC_VECTOR (4 downto 0);
    signal grp_fu_100_ce : STD_LOGIC;
    signal ap_NS_fsm : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_idle_pp0 : STD_LOGIC;

    component tkmu_simple_hw_mueOg IS
    generic (
        ID : INTEGER;
        NUM_STAGE : INTEGER;
        din0_WIDTH : INTEGER;
        din1_WIDTH : INTEGER;
        dout_WIDTH : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        din0 : IN STD_LOGIC_VECTOR (22 downto 0);
        din1 : IN STD_LOGIC_VECTOR (20 downto 0);
        ce : IN STD_LOGIC;
        dout : OUT STD_LOGIC_VECTOR (43 downto 0) );
    end component;


    component delta_plus_LUT_defYi IS
    generic (
        DataWidth : INTEGER;
        AddressRange : INTEGER;
        AddressWidth : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        address0 : IN STD_LOGIC_VECTOR (9 downto 0);
        ce0 : IN STD_LOGIC;
        q0 : OUT STD_LOGIC_VECTOR (4 downto 0) );
    end component;



begin
    delta_plus_table4_U : component delta_plus_LUT_defYi
    generic map (
        DataWidth => 5,
        AddressRange => 1024,
        AddressWidth => 10)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        address0 => delta_plus_table4_address0,
        ce0 => delta_plus_table4_ce0,
        q0 => delta_plus_table4_q0);

    tkmu_simple_hw_mueOg_x_U7 : component tkmu_simple_hw_mueOg
    generic map (
        ID => 1,
        NUM_STAGE => 5,
        din0_WIDTH => 23,
        din1_WIDTH => 21,
        dout_WIDTH => 44)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        din0 => grp_fu_100_p0,
        din1 => r_V_fu_88_p3,
        ce => grp_fu_100_ce,
        dout => grp_fu_100_p2);





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
                if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce)))))) then 
                    ap_enable_reg_pp0_iter1 <= ap_start;
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
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
                if (not((((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)) or not((ap_const_logic_1 = ap_ce))))) then 
                    ap_enable_reg_pp0_iter7 <= ap_enable_reg_pp0_iter6;
                end if; 
            end if;
        end if;
    end process;

    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce))) then
                ap_pipeline_reg_pp0_iter1_tmp_250_reg_307 <= tmp_250_reg_307;
                tmp_250_reg_307 <= data_V_read(10 downto 10);
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce))) then
                ap_pipeline_reg_pp0_iter2_tmp_250_reg_307 <= ap_pipeline_reg_pp0_iter1_tmp_250_reg_307;
                ap_pipeline_reg_pp0_iter3_tmp_250_reg_307 <= ap_pipeline_reg_pp0_iter2_tmp_250_reg_307;
                ap_pipeline_reg_pp0_iter4_tmp_250_reg_307 <= ap_pipeline_reg_pp0_iter3_tmp_250_reg_307;
                icmp_reg_349 <= icmp_fu_260_p2;
                p_v_reg_313 <= p_v_fu_152_p3;
                r_V_15_reg_323 <= r_V_15_fu_178_p2;
                ret_V_reg_338 <= ret_V_fu_208_p2;
                tmp_204_reg_328 <= r_V_15_fu_178_p2(18 downto 6);
                tmp_257_reg_343 <= index_fu_234_p3(31 downto 31);
                tmp_s_reg_333 <= tmp_s_fu_202_p2;
            end if;
        end if;
    end process;
    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce) and not((ap_pipeline_reg_pp0_iter3_tmp_250_reg_307 = ap_const_lv1_0)))) then
                tmp_253_reg_318 <= tmp_253_fu_159_p1;
            end if;
        end if;
    end process;

    ap_NS_fsm_assign_proc : process (ap_start, ap_CS_fsm, ap_enable_reg_pp0_iter0, ap_ce, ap_pipeline_idle_pp0)
    begin
        case ap_CS_fsm is
            when ap_ST_fsm_pp0_stage0 => 
                ap_NS_fsm <= ap_ST_fsm_pp0_stage0;
            when others =>  
                ap_NS_fsm <= "X";
        end case;
    end process;
    ap_CS_fsm_pp0_stage0 <= ap_CS_fsm(0 downto 0);

    ap_done_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter7, ap_ce)
    begin
        if ((((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_1 = ap_enable_reg_pp0_iter0)) or (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce) and (ap_const_logic_1 = ap_enable_reg_pp0_iter7)))) then 
            ap_done <= ap_const_logic_1;
        else 
            ap_done <= ap_const_logic_0;
        end if; 
    end process;

    ap_enable_reg_pp0_iter0 <= ap_start;

    ap_idle_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5, ap_enable_reg_pp0_iter6, ap_enable_reg_pp0_iter7)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5) and (ap_const_logic_0 = ap_enable_reg_pp0_iter6) and (ap_const_logic_0 = ap_enable_reg_pp0_iter7))) then 
            ap_idle <= ap_const_logic_1;
        else 
            ap_idle <= ap_const_logic_0;
        end if; 
    end process;


    ap_pipeline_idle_pp0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5, ap_enable_reg_pp0_iter6)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5) and (ap_const_logic_0 = ap_enable_reg_pp0_iter6))) then 
            ap_pipeline_idle_pp0 <= ap_const_logic_1;
        else 
            ap_pipeline_idle_pp0 <= ap_const_logic_0;
        end if; 
    end process;


    ap_ready_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_ce)
    begin
        if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_1 = ap_enable_reg_pp0_iter0) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce))) then 
            ap_ready <= ap_const_logic_1;
        else 
            ap_ready <= ap_const_logic_0;
        end if; 
    end process;

    ap_return <= 
        sel_tmp_fu_281_p3 when (tmp_206_fu_289_p2(0) = '1') else 
        delta_plus_table4_q0;
    delta_plus_table4_address0 <= tmp_202_fu_266_p1(10 - 1 downto 0);

    delta_plus_table4_ce0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter6, ap_ce)
    begin
        if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce) and (ap_const_logic_1 = ap_enable_reg_pp0_iter6))) then 
            delta_plus_table4_ce0 <= ap_const_logic_1;
        else 
            delta_plus_table4_ce0 <= ap_const_logic_0;
        end if; 
    end process;


    grp_fu_100_ce_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_ce)
    begin
        if (((ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_ce))) then 
            grp_fu_100_ce <= ap_const_logic_1;
        else 
            grp_fu_100_ce <= ap_const_logic_0;
        end if; 
    end process;

    grp_fu_100_p0 <= ap_const_lv44_222223(23 - 1 downto 0);
    icmp_fu_260_p2 <= "1" when (signed(tmp_258_fu_250_p4) > signed(ap_const_lv22_0)) else "0";
    index_fu_234_p3 <= 
        tmp_213_fu_227_p3 when (tmp_255_fu_214_p3(0) = '1') else 
        tmp_211_fu_221_p1;
    neg_mul_fu_118_p2 <= std_logic_vector(unsigned(ap_const_lv43_0) - unsigned(tmp_249_fu_114_p1));
    neg_ti_fu_163_p2 <= std_logic_vector(unsigned(ap_const_lv19_0) - unsigned(tmp_253_reg_318));
    p_v_fu_152_p3 <= 
        tmp_fu_134_p1 when (ap_pipeline_reg_pp0_iter3_tmp_250_reg_307(0) = '1') else 
        tmp_209_fu_148_p1;
    r_V_15_fu_178_p2 <= std_logic_vector(unsigned(ap_const_lv19_10000) - unsigned(tmp_210_fu_171_p3));
    r_V_fu_88_p3 <= (data_V_read & ap_const_lv10_0);
        ret_V_cast_fu_194_p1 <= std_logic_vector(resize(signed(tmp_204_fu_184_p4),14));

    ret_V_fu_208_p2 <= std_logic_vector(unsigned(ap_const_lv14_1) + unsigned(ret_V_cast_fu_194_p1));
    sel_tmp1_fu_271_p2 <= (tmp_257_reg_343 xor ap_const_lv1_1);
    sel_tmp2_fu_276_p2 <= (icmp_reg_349 and sel_tmp1_fu_271_p2);
    sel_tmp_fu_281_p3 <= 
        ap_const_lv5_0 when (sel_tmp2_fu_276_p2(0) = '1') else 
        ap_const_lv5_11;
    tmp_202_fu_266_p1 <= std_logic_vector(resize(unsigned(index_fu_234_p3),64));
    tmp_204_fu_184_p4 <= r_V_15_fu_178_p2(18 downto 6);
    tmp_206_fu_289_p2 <= (sel_tmp2_fu_276_p2 or tmp_257_reg_343);
        tmp_209_fu_148_p1 <= std_logic_vector(resize(signed(tmp_252_fu_138_p4),22));

    tmp_210_fu_171_p3 <= 
        neg_ti_fu_163_p2 when (ap_pipeline_reg_pp0_iter4_tmp_250_reg_307(0) = '1') else 
        tmp_254_fu_168_p1;
        tmp_211_fu_221_p1 <= std_logic_vector(resize(signed(tmp_204_reg_328),32));

        tmp_212_fu_224_p1 <= std_logic_vector(resize(signed(ret_V_reg_338),32));

    tmp_213_fu_227_p3 <= 
        tmp_211_fu_221_p1 when (tmp_s_reg_333(0) = '1') else 
        tmp_212_fu_224_p1;
    tmp_249_fu_114_p1 <= grp_fu_100_p2(43 - 1 downto 0);
    tmp_251_fu_124_p4 <= neg_mul_fu_118_p2(42 downto 25);
    tmp_252_fu_138_p4 <= grp_fu_100_p2(43 downto 25);
    tmp_253_fu_159_p1 <= p_v_fu_152_p3(19 - 1 downto 0);
    tmp_254_fu_168_p1 <= p_v_reg_313(19 - 1 downto 0);
    tmp_255_fu_214_p3 <= r_V_15_reg_323(18 downto 18);
    tmp_256_fu_198_p1 <= r_V_15_fu_178_p2(6 - 1 downto 0);
    tmp_258_fu_250_p4 <= index_fu_234_p3(31 downto 10);
        tmp_fu_134_p1 <= std_logic_vector(resize(signed(tmp_251_fu_124_p4),22));

    tmp_s_fu_202_p2 <= "1" when (tmp_256_fu_198_p1 = ap_const_lv6_0) else "0";
end behav;
