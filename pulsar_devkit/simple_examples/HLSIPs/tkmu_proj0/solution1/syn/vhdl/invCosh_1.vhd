-- ==============================================================
-- RTL generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
-- Version: 2016.4
-- Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
-- 
-- ===========================================================

library IEEE;
use IEEE.std_logic_1164.all;
use IEEE.numeric_std.all;

entity invCosh_1 is
port (
    ap_clk : IN STD_LOGIC;
    ap_rst : IN STD_LOGIC;
    ap_start : IN STD_LOGIC;
    ap_done : OUT STD_LOGIC;
    ap_idle : OUT STD_LOGIC;
    ap_ready : OUT STD_LOGIC;
    data_V_read : IN STD_LOGIC_VECTOR (11 downto 0);
    ap_return : OUT STD_LOGIC_VECTOR (10 downto 0) );
end;


architecture behav of invCosh_1 is 
    constant ap_const_logic_1 : STD_LOGIC := '1';
    constant ap_const_logic_0 : STD_LOGIC := '0';
    constant ap_ST_fsm_pp0_stage0 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv32_0 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000000000";
    constant ap_const_lv1_1 : STD_LOGIC_VECTOR (0 downto 0) := "1";
    constant ap_const_lv13_0 : STD_LOGIC_VECTOR (12 downto 0) := "0000000000000";
    constant ap_const_lv49_1555556 : STD_LOGIC_VECTOR (48 downto 0) := "0000000000000000000000001010101010101010101010110";
    constant ap_const_lv32_1A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011010";
    constant ap_const_lv32_30 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000110000";
    constant ap_const_lv25_800000 : STD_LOGIC_VECTOR (24 downto 0) := "0100000000000000000000000";
    constant ap_const_lv32_A : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000001010";
    constant ap_const_lv32_18 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000011000";
    constant ap_const_lv32_17 : STD_LOGIC_VECTOR (31 downto 0) := "00000000000000000000000000010111";
    constant ap_const_lv2_0 : STD_LOGIC_VECTOR (1 downto 0) := "00";
    constant ap_const_lv11_400 : STD_LOGIC_VECTOR (10 downto 0) := "10000000000";

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
    signal cosh_table7_address0 : STD_LOGIC_VECTOR (12 downto 0);
    signal cosh_table7_ce0 : STD_LOGIC;
    signal cosh_table7_q0 : STD_LOGIC_VECTOR (10 downto 0);
    signal tmp_232_reg_152 : STD_LOGIC_VECTOR (14 downto 0);
    signal icmp_fu_126_p2 : STD_LOGIC_VECTOR (0 downto 0);
    signal icmp_reg_157 : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_reg_pp0_iter5_icmp_reg_157 : STD_LOGIC_VECTOR (0 downto 0);
    signal tmp_196_fu_135_p1 : STD_LOGIC_VECTOR (63 downto 0);
    signal tmp_fu_64_p1 : STD_LOGIC_VECTOR (10 downto 0);
    signal r_V_tr_fu_68_p3 : STD_LOGIC_VECTOR (23 downto 0);
    signal grp_fu_80_p0 : STD_LOGIC_VECTOR (25 downto 0);
    signal grp_fu_80_p1 : STD_LOGIC_VECTOR (23 downto 0);
    signal grp_fu_80_p2 : STD_LOGIC_VECTOR (48 downto 0);
    signal tmp_231_fu_86_p4 : STD_LOGIC_VECTOR (22 downto 0);
    signal tmp_cast_cast_fu_96_p1 : STD_LOGIC_VECTOR (24 downto 0);
    signal r_V_fu_100_p2 : STD_LOGIC_VECTOR (24 downto 0);
    signal tmp_233_fu_116_p4 : STD_LOGIC_VECTOR (1 downto 0);
    signal tmp_198_fu_132_p1 : STD_LOGIC_VECTOR (15 downto 0);
    signal ap_NS_fsm : STD_LOGIC_VECTOR (0 downto 0);
    signal ap_pipeline_idle_pp0 : STD_LOGIC;
    signal grp_fu_80_p10 : STD_LOGIC_VECTOR (48 downto 0);

    component tkmu_simple_hw_mukbM IS
    generic (
        ID : INTEGER;
        NUM_STAGE : INTEGER;
        din0_WIDTH : INTEGER;
        din1_WIDTH : INTEGER;
        dout_WIDTH : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        din0 : IN STD_LOGIC_VECTOR (25 downto 0);
        din1 : IN STD_LOGIC_VECTOR (23 downto 0);
        ce : IN STD_LOGIC;
        dout : OUT STD_LOGIC_VECTOR (48 downto 0) );
    end component;


    component invCosh_1_cosh_tajbC IS
    generic (
        DataWidth : INTEGER;
        AddressRange : INTEGER;
        AddressWidth : INTEGER );
    port (
        clk : IN STD_LOGIC;
        reset : IN STD_LOGIC;
        address0 : IN STD_LOGIC_VECTOR (12 downto 0);
        ce0 : IN STD_LOGIC;
        q0 : OUT STD_LOGIC_VECTOR (10 downto 0) );
    end component;



begin
    cosh_table7_U : component invCosh_1_cosh_tajbC
    generic map (
        DataWidth => 11,
        AddressRange => 8192,
        AddressWidth => 13)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        address0 => cosh_table7_address0,
        ce0 => cosh_table7_ce0,
        q0 => cosh_table7_q0);

    tkmu_simple_hw_mukbM_U19 : component tkmu_simple_hw_mukbM
    generic map (
        ID => 1,
        NUM_STAGE => 5,
        din0_WIDTH => 26,
        din1_WIDTH => 24,
        dout_WIDTH => 49)
    port map (
        clk => ap_clk,
        reset => ap_rst,
        din0 => grp_fu_80_p0,
        din1 => grp_fu_80_p1,
        ce => ap_const_logic_1,
        dout => grp_fu_80_p2);





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

    process (ap_clk)
    begin
        if (ap_clk'event and ap_clk = '1') then
            if (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0)))) then
                ap_pipeline_reg_pp0_iter5_icmp_reg_157 <= icmp_reg_157;
                icmp_reg_157 <= icmp_fu_126_p2;
                tmp_232_reg_152 <= r_V_fu_100_p2(24 downto 10);
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

    ap_done_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter6)
    begin
        if ((((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_1 = ap_enable_reg_pp0_iter0)) or (not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_enable_reg_pp0_iter6)))) then 
            ap_done <= ap_const_logic_1;
        else 
            ap_done <= ap_const_logic_0;
        end if; 
    end process;

    ap_enable_reg_pp0_iter0 <= ap_start;

    ap_idle_assign_proc : process(ap_start, ap_CS_fsm_pp0_stage0, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5, ap_enable_reg_pp0_iter6)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_CS_fsm_pp0_stage0 = ap_const_lv1_1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5) and (ap_const_logic_0 = ap_enable_reg_pp0_iter6))) then 
            ap_idle <= ap_const_logic_1;
        else 
            ap_idle <= ap_const_logic_0;
        end if; 
    end process;


    ap_pipeline_idle_pp0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter1, ap_enable_reg_pp0_iter2, ap_enable_reg_pp0_iter3, ap_enable_reg_pp0_iter4, ap_enable_reg_pp0_iter5)
    begin
        if (((ap_const_logic_0 = ap_start) and (ap_const_logic_0 = ap_enable_reg_pp0_iter0) and (ap_const_logic_0 = ap_enable_reg_pp0_iter1) and (ap_const_logic_0 = ap_enable_reg_pp0_iter2) and (ap_const_logic_0 = ap_enable_reg_pp0_iter3) and (ap_const_logic_0 = ap_enable_reg_pp0_iter4) and (ap_const_logic_0 = ap_enable_reg_pp0_iter5))) then 
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
        ap_const_lv11_400 when (ap_pipeline_reg_pp0_iter5_icmp_reg_157(0) = '1') else 
        cosh_table7_q0;
    cosh_table7_address0 <= tmp_196_fu_135_p1(13 - 1 downto 0);

    cosh_table7_ce0_assign_proc : process(ap_start, ap_enable_reg_pp0_iter0, ap_enable_reg_pp0_iter5)
    begin
        if ((not(((ap_const_logic_1 = ap_enable_reg_pp0_iter0) and (ap_start = ap_const_logic_0))) and (ap_const_logic_1 = ap_enable_reg_pp0_iter5))) then 
            cosh_table7_ce0 <= ap_const_logic_1;
        else 
            cosh_table7_ce0 <= ap_const_logic_0;
        end if; 
    end process;

    grp_fu_80_p0 <= ap_const_lv49_1555556(26 - 1 downto 0);
    grp_fu_80_p1 <= grp_fu_80_p10(24 - 1 downto 0);
    grp_fu_80_p10 <= std_logic_vector(resize(unsigned(r_V_tr_fu_68_p3),49));
    icmp_fu_126_p2 <= "0" when (tmp_233_fu_116_p4 = ap_const_lv2_0) else "1";
    r_V_fu_100_p2 <= std_logic_vector(unsigned(ap_const_lv25_800000) - unsigned(tmp_cast_cast_fu_96_p1));
    r_V_tr_fu_68_p3 <= (tmp_fu_64_p1 & ap_const_lv13_0);
    tmp_196_fu_135_p1 <= std_logic_vector(resize(unsigned(tmp_198_fu_132_p1),64));
        tmp_198_fu_132_p1 <= std_logic_vector(resize(signed(tmp_232_reg_152),16));

    tmp_231_fu_86_p4 <= grp_fu_80_p2(48 downto 26);
    tmp_233_fu_116_p4 <= r_V_fu_100_p2(24 downto 23);
    tmp_cast_cast_fu_96_p1 <= std_logic_vector(resize(unsigned(tmp_231_fu_86_p4),25));
    tmp_fu_64_p1 <= data_V_read(11 - 1 downto 0);
end behav;
