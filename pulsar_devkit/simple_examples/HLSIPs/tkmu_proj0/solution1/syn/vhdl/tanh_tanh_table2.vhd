-- ==============================================================
-- File generated by Vivado(TM) HLS - High-Level Synthesis from C, C++ and SystemC
-- Version: 2016.4
-- Copyright (C) 1986-2017 Xilinx, Inc. All Rights Reserved.
-- 
-- ==============================================================

library ieee; 
use ieee.std_logic_1164.all; 
use ieee.std_logic_unsigned.all;

entity tanh_tanh_table2_rom is 
    generic(
             dwidth     : integer := 10; 
             awidth     : integer := 13; 
             mem_size    : integer := 8192
    ); 
    port (
          addr0      : in std_logic_vector(awidth-1 downto 0); 
          ce0       : in std_logic; 
          q0         : out std_logic_vector(dwidth-1 downto 0);
          clk       : in std_logic
    ); 
end entity; 


architecture rtl of tanh_tanh_table2_rom is 

signal addr0_tmp : std_logic_vector(awidth-1 downto 0); 
type mem_array is array (0 to mem_size-1) of std_logic_vector (dwidth-1 downto 0); 
signal mem : mem_array := (
    0 to 231=> "1111111010", 232 to 442=> "1111111001", 443 to 625=> "1111111000", 626 to 786=> "1111110111", 
    787 to 931=> "1111110110", 932 to 1062=> "1111110101", 1063 to 1181=> "1111110100", 1182 to 1291=> "1111110011", 
    1292 to 1393=> "1111110010", 1394 to 1488=> "1111110001", 1489 to 1577=> "1111110000", 1578 to 1660=> "1111101111", 
    1661 to 1739=> "1111101110", 1740 to 1813=> "1111101101", 1814 to 1884=> "1111101100", 1885 to 1951=> "1111101011", 
    1952 to 2015=> "1111101010", 2016 to 2077=> "1111101001", 2078 to 2136=> "1111101000", 2137 to 2192=> "1111100111", 
    2193 to 2246=> "1111100110", 2247 to 2298=> "1111100101", 2299 to 2349=> "1111100100", 2350 to 2397=> "1111100011", 
    2398 to 2444=> "1111100010", 2445 to 2490=> "1111100001", 2491 to 2534=> "1111100000", 2535 to 2576=> "1111011111", 
    2577 to 2618=> "1111011110", 2619 to 2658=> "1111011101", 2659 to 2697=> "1111011100", 2698 to 2735=> "1111011011", 
    2736 to 2772=> "1111011010", 2773 to 2809=> "1111011001", 2810 to 2844=> "1111011000", 2845 to 2878=> "1111010111", 
    2879 to 2912=> "1111010110", 2913 to 2945=> "1111010101", 2946 to 2977=> "1111010100", 2978 to 3008=> "1111010011", 
    3009 to 3039=> "1111010010", 3040 to 3069=> "1111010001", 3070 to 3098=> "1111010000", 3099 to 3127=> "1111001111", 
    3128 to 3155=> "1111001110", 3156 to 3183=> "1111001101", 3184 to 3210=> "1111001100", 3211 to 3237=> "1111001011", 
    3238 to 3263=> "1111001010", 3264 to 3289=> "1111001001", 3290 to 3314=> "1111001000", 3315 to 3339=> "1111000111", 
    3340 to 3363=> "1111000110", 3364 to 3387=> "1111000101", 3388 to 3411=> "1111000100", 3412 to 3434=> "1111000011", 
    3435 to 3457=> "1111000010", 3458 to 3480=> "1111000001", 3481 to 3502=> "1111000000", 3503 to 3524=> "1110111111", 
    3525 to 3545=> "1110111110", 3546 to 3567=> "1110111101", 3568 to 3587=> "1110111100", 3588 to 3608=> "1110111011", 
    3609 to 3628=> "1110111010", 3629 to 3649=> "1110111001", 3650 to 3668=> "1110111000", 3669 to 3688=> "1110110111", 
    3689 to 3707=> "1110110110", 3708 to 3726=> "1110110101", 3727 to 3745=> "1110110100", 3746 to 3763=> "1110110011", 
    3764 to 3782=> "1110110010", 3783 to 3800=> "1110110001", 3801 to 3818=> "1110110000", 3819 to 3835=> "1110101111", 
    3836 to 3853=> "1110101110", 3854 to 3870=> "1110101101", 3871 to 3887=> "1110101100", 3888 to 3904=> "1110101011", 
    3905 to 3921=> "1110101010", 3922 to 3937=> "1110101001", 3938 to 3953=> "1110101000", 3954 to 3969=> "1110100111", 
    3970 to 3985=> "1110100110", 3986 to 4001=> "1110100101", 4002 to 4017=> "1110100100", 4018 to 4032=> "1110100011", 
    4033 to 4048=> "1110100010", 4049 to 4063=> "1110100001", 4064 to 4078=> "1110100000", 4079 to 4093=> "1110011111", 
    4094 to 4107=> "1110011110", 4108 to 4122=> "1110011101", 4123 to 4136=> "1110011100", 4137 to 4151=> "1110011011", 
    4152 to 4165=> "1110011010", 4166 to 4179=> "1110011001", 4180 to 4193=> "1110011000", 4194 to 4206=> "1110010111", 
    4207 to 4220=> "1110010110", 4221 to 4234=> "1110010101", 4235 to 4247=> "1110010100", 4248 to 4260=> "1110010011", 
    4261 to 4273=> "1110010010", 4274 to 4287=> "1110010001", 4288 to 4299=> "1110010000", 4300 to 4312=> "1110001111", 
    4313 to 4325=> "1110001110", 4326 to 4338=> "1110001101", 4339 to 4350=> "1110001100", 4351 to 4363=> "1110001011", 
    4364 to 4375=> "1110001010", 4376 to 4387=> "1110001001", 4388 to 4399=> "1110001000", 4400 to 4411=> "1110000111", 
    4412 to 4423=> "1110000110", 4424 to 4435=> "1110000101", 4436 to 4447=> "1110000100", 4448 to 4459=> "1110000011", 
    4460 to 4470=> "1110000010", 4471 to 4482=> "1110000001", 4483 to 4493=> "1110000000", 4494 to 4504=> "1101111111", 
    4505 to 4516=> "1101111110", 4517 to 4527=> "1101111101", 4528 to 4538=> "1101111100", 4539 to 4549=> "1101111011", 
    4550 to 4560=> "1101111010", 4561 to 4571=> "1101111001", 4572 to 4582=> "1101111000", 4583 to 4592=> "1101110111", 
    4593 to 4603=> "1101110110", 4604 to 4614=> "1101110101", 4615 to 4624=> "1101110100", 4625 to 4634=> "1101110011", 
    4635 to 4645=> "1101110010", 4646 to 4655=> "1101110001", 4656 to 4665=> "1101110000", 4666 to 4676=> "1101101111", 
    4677 to 4686=> "1101101110", 4687 to 4696=> "1101101101", 4697 to 4706=> "1101101100", 4707 to 4716=> "1101101011", 
    4717 to 4725=> "1101101010", 4726 to 4735=> "1101101001", 4736 to 4745=> "1101101000", 4746 to 4755=> "1101100111", 
    4756 to 4764=> "1101100110", 4765 to 4774=> "1101100101", 4775 to 4783=> "1101100100", 4784 to 4793=> "1101100011", 
    4794 to 4802=> "1101100010", 4803 to 4811=> "1101100001", 4812 to 4821=> "1101100000", 4822 to 4830=> "1101011111", 
    4831 to 4839=> "1101011110", 4840 to 4848=> "1101011101", 4849 to 4857=> "1101011100", 4858 to 4866=> "1101011011", 
    4867 to 4875=> "1101011010", 4876 to 4884=> "1101011001", 4885 to 4893=> "1101011000", 4894 to 4902=> "1101010111", 
    4903 to 4911=> "1101010110", 4912 to 4919=> "1101010101", 4920 to 4928=> "1101010100", 4929 to 4937=> "1101010011", 
    4938 to 4945=> "1101010010", 4946 to 4954=> "1101010001", 4955 to 4962=> "1101010000", 4963 to 4971=> "1101001111", 
    4972 to 4979=> "1101001110", 4980 to 4988=> "1101001101", 4989 to 4996=> "1101001100", 4997 to 5004=> "1101001011", 
    5005 to 5013=> "1101001010", 5014 to 5021=> "1101001001", 5022 to 5029=> "1101001000", 5030 to 5037=> "1101000111", 
    5038 to 5045=> "1101000110", 5046 to 5053=> "1101000101", 5054 to 5061=> "1101000100", 5062 to 5069=> "1101000011", 
    5070 to 5077=> "1101000010", 5078 to 5085=> "1101000001", 5086 to 5093=> "1101000000", 5094 to 5101=> "1100111111", 
    5102 to 5109=> "1100111110", 5110 to 5116=> "1100111101", 5117 to 5124=> "1100111100", 5125 to 5132=> "1100111011", 
    5133 to 5139=> "1100111010", 5140 to 5147=> "1100111001", 5148 to 5155=> "1100111000", 5156 to 5162=> "1100110111", 
    5163 to 5170=> "1100110110", 5171 to 5177=> "1100110101", 5178 to 5185=> "1100110100", 5186 to 5192=> "1100110011", 
    5193 to 5199=> "1100110010", 5200 to 5207=> "1100110001", 5208 to 5214=> "1100110000", 5215 to 5221=> "1100101111", 
    5222 to 5229=> "1100101110", 5230 to 5236=> "1100101101", 5237 to 5243=> "1100101100", 5244 to 5250=> "1100101011", 
    5251 to 5257=> "1100101010", 5258 to 5264=> "1100101001", 5265 to 5272=> "1100101000", 5273 to 5279=> "1100100111", 
    5280 to 5286=> "1100100110", 5287 to 5293=> "1100100101", 5294 to 5300=> "1100100100", 5301 to 5307=> "1100100011", 
    5308 to 5313=> "1100100010", 5314 to 5320=> "1100100001", 5321 to 5327=> "1100100000", 5328 to 5334=> "1100011111", 
    5335 to 5341=> "1100011110", 5342 to 5348=> "1100011101", 5349 to 5354=> "1100011100", 5355 to 5361=> "1100011011", 
    5362 to 5368=> "1100011010", 5369 to 5374=> "1100011001", 5375 to 5381=> "1100011000", 5382 to 5388=> "1100010111", 
    5389 to 5394=> "1100010110", 5395 to 5401=> "1100010101", 5402 to 5407=> "1100010100", 5408 to 5414=> "1100010011", 
    5415 to 5420=> "1100010010", 5421 to 5427=> "1100010001", 5428 to 5433=> "1100010000", 5434 to 5440=> "1100001111", 
    5441 to 5446=> "1100001110", 5447 to 5453=> "1100001101", 5454 to 5459=> "1100001100", 5460 to 5465=> "1100001011", 
    5466 to 5472=> "1100001010", 5473 to 5478=> "1100001001", 5479 to 5484=> "1100001000", 5485 to 5491=> "1100000111", 
    5492 to 5497=> "1100000110", 5498 to 5503=> "1100000101", 5504 to 5509=> "1100000100", 5510 to 5515=> "1100000011", 
    5516 to 5521=> "1100000010", 5522 to 5528=> "1100000001", 5529 to 5534=> "1100000000", 5535 to 5540=> "1011111111", 
    5541 to 5546=> "1011111110", 5547 to 5552=> "1011111101", 5553 to 5558=> "1011111100", 5559 to 5564=> "1011111011", 
    5565 to 5570=> "1011111010", 5571 to 5576=> "1011111001", 5577 to 5582=> "1011111000", 5583 to 5588=> "1011110111", 
    5589 to 5594=> "1011110110", 5595 to 5600=> "1011110101", 5601 to 5605=> "1011110100", 5606 to 5611=> "1011110011", 
    5612 to 5617=> "1011110010", 5618 to 5623=> "1011110001", 5624 to 5629=> "1011110000", 5630 to 5634=> "1011101111", 
    5635 to 5640=> "1011101110", 5641 to 5646=> "1011101101", 5647 to 5652=> "1011101100", 5653 to 5657=> "1011101011", 
    5658 to 5663=> "1011101010", 5664 to 5669=> "1011101001", 5670 to 5674=> "1011101000", 5675 to 5680=> "1011100111", 
    5681 to 5686=> "1011100110", 5687 to 5691=> "1011100101", 5692 to 5697=> "1011100100", 5698 to 5702=> "1011100011", 
    5703 to 5708=> "1011100010", 5709 to 5714=> "1011100001", 5715 to 5719=> "1011100000", 5720 to 5725=> "1011011111", 
    5726 to 5730=> "1011011110", 5731 to 5736=> "1011011101", 5737 to 5741=> "1011011100", 5742 to 5746=> "1011011011", 
    5747 to 5752=> "1011011010", 5753 to 5757=> "1011011001", 5758 to 5763=> "1011011000", 5764 to 5768=> "1011010111", 
    5769 to 5773=> "1011010110", 5774 to 5779=> "1011010101", 5780 to 5784=> "1011010100", 5785 to 5789=> "1011010011", 
    5790 to 5795=> "1011010010", 5796 to 5800=> "1011010001", 5801 to 5805=> "1011010000", 5806 to 5811=> "1011001111", 
    5812 to 5816=> "1011001110", 5817 to 5821=> "1011001101", 5822 to 5826=> "1011001100", 5827 to 5832=> "1011001011", 
    5833 to 5837=> "1011001010", 5838 to 5842=> "1011001001", 5843 to 5847=> "1011001000", 5848 to 5852=> "1011000111", 
    5853 to 5857=> "1011000110", 5858 to 5863=> "1011000101", 5864 to 5868=> "1011000100", 5869 to 5873=> "1011000011", 
    5874 to 5878=> "1011000010", 5879 to 5883=> "1011000001", 5884 to 5888=> "1011000000", 5889 to 5893=> "1010111111", 
    5894 to 5898=> "1010111110", 5899 to 5903=> "1010111101", 5904 to 5908=> "1010111100", 5909 to 5913=> "1010111011", 
    5914 to 5918=> "1010111010", 5919 to 5923=> "1010111001", 5924 to 5928=> "1010111000", 5929 to 5933=> "1010110111", 
    5934 to 5938=> "1010110110", 5939 to 5943=> "1010110101", 5944 to 5948=> "1010110100", 5949 to 5953=> "1010110011", 
    5954 to 5958=> "1010110010", 5959 to 5962=> "1010110001", 5963 to 5967=> "1010110000", 5968 to 5972=> "1010101111", 
    5973 to 5977=> "1010101110", 5978 to 5982=> "1010101101", 5983 to 5987=> "1010101100", 5988 to 5991=> "1010101011", 
    5992 to 5996=> "1010101010", 5997 to 6001=> "1010101001", 6002 to 6006=> "1010101000", 6007 to 6011=> "1010100111", 
    6012 to 6015=> "1010100110", 6016 to 6020=> "1010100101", 6021 to 6025=> "1010100100", 6026 to 6030=> "1010100011", 
    6031 to 6034=> "1010100010", 6035 to 6039=> "1010100001", 6040 to 6044=> "1010100000", 6045 to 6048=> "1010011111", 
    6049 to 6053=> "1010011110", 6054 to 6058=> "1010011101", 6059 to 6062=> "1010011100", 6063 to 6067=> "1010011011", 
    6068 to 6072=> "1010011010", 6073 to 6076=> "1010011001", 6077 to 6081=> "1010011000", 6082 to 6085=> "1010010111", 
    6086 to 6090=> "1010010110", 6091 to 6095=> "1010010101", 6096 to 6099=> "1010010100", 6100 to 6104=> "1010010011", 
    6105 to 6108=> "1010010010", 6109 to 6113=> "1010010001", 6114 to 6117=> "1010010000", 6118 to 6122=> "1010001111", 
    6123 to 6126=> "1010001110", 6127 to 6131=> "1010001101", 6132 to 6135=> "1010001100", 6136 to 6140=> "1010001011", 
    6141 to 6144=> "1010001010", 6145 to 6149=> "1010001001", 6150 to 6153=> "1010001000", 6154 to 6158=> "1010000111", 
    6159 to 6162=> "1010000110", 6163 to 6166=> "1010000101", 6167 to 6171=> "1010000100", 6172 to 6175=> "1010000011", 
    6176 to 6180=> "1010000010", 6181 to 6184=> "1010000001", 6185 to 6188=> "1010000000", 6189 to 6193=> "1001111111", 
    6194 to 6197=> "1001111110", 6198 to 6202=> "1001111101", 6203 to 6206=> "1001111100", 6207 to 6210=> "1001111011", 
    6211 to 6215=> "1001111010", 6216 to 6219=> "1001111001", 6220 to 6223=> "1001111000", 6224 to 6227=> "1001110111", 
    6228 to 6232=> "1001110110", 6233 to 6236=> "1001110101", 6237 to 6240=> "1001110100", 6241 to 6245=> "1001110011", 
    6246 to 6249=> "1001110010", 6250 to 6253=> "1001110001", 6254 to 6257=> "1001110000", 6258 to 6262=> "1001101111", 
    6263 to 6266=> "1001101110", 6267 to 6270=> "1001101101", 6271 to 6274=> "1001101100", 6275 to 6278=> "1001101011", 
    6279 to 6283=> "1001101010", 6284 to 6287=> "1001101001", 6288 to 6291=> "1001101000", 6292 to 6295=> "1001100111", 
    6296 to 6299=> "1001100110", 6300 to 6304=> "1001100101", 6305 to 6308=> "1001100100", 6309 to 6312=> "1001100011", 
    6313 to 6316=> "1001100010", 6317 to 6320=> "1001100001", 6321 to 6324=> "1001100000", 6325 to 6328=> "1001011111", 
    6329 to 6332=> "1001011110", 6333 to 6337=> "1001011101", 6338 to 6341=> "1001011100", 6342 to 6345=> "1001011011", 
    6346 to 6349=> "1001011010", 6350 to 6353=> "1001011001", 6354 to 6357=> "1001011000", 6358 to 6361=> "1001010111", 
    6362 to 6365=> "1001010110", 6366 to 6369=> "1001010101", 6370 to 6373=> "1001010100", 6374 to 6377=> "1001010011", 
    6378 to 6381=> "1001010010", 6382 to 6385=> "1001010001", 6386 to 6389=> "1001010000", 6390 to 6393=> "1001001111", 
    6394 to 6397=> "1001001110", 6398 to 6401=> "1001001101", 6402 to 6405=> "1001001100", 6406 to 6409=> "1001001011", 
    6410 to 6413=> "1001001010", 6414 to 6417=> "1001001001", 6418 to 6421=> "1001001000", 6422 to 6425=> "1001000111", 
    6426 to 6429=> "1001000110", 6430 to 6433=> "1001000101", 6434 to 6437=> "1001000100", 6438 to 6441=> "1001000011", 
    6442 to 6445=> "1001000010", 6446 to 6449=> "1001000001", 6450 to 6452=> "1001000000", 6453 to 6456=> "1000111111", 
    6457 to 6460=> "1000111110", 6461 to 6464=> "1000111101", 6465 to 6468=> "1000111100", 6469 to 6472=> "1000111011", 
    6473 to 6476=> "1000111010", 6477 to 6480=> "1000111001", 6481 to 6483=> "1000111000", 6484 to 6487=> "1000110111", 
    6488 to 6491=> "1000110110", 6492 to 6495=> "1000110101", 6496 to 6499=> "1000110100", 6500 to 6503=> "1000110011", 
    6504 to 6506=> "1000110010", 6507 to 6510=> "1000110001", 6511 to 6514=> "1000110000", 6515 to 6518=> "1000101111", 
    6519 to 6522=> "1000101110", 6523 to 6525=> "1000101101", 6526 to 6529=> "1000101100", 6530 to 6533=> "1000101011", 
    6534 to 6537=> "1000101010", 6538 to 6541=> "1000101001", 6542 to 6544=> "1000101000", 6545 to 6548=> "1000100111", 
    6549 to 6552=> "1000100110", 6553 to 6556=> "1000100101", 6557 to 6559=> "1000100100", 6560 to 6563=> "1000100011", 
    6564 to 6567=> "1000100010", 6568 to 6571=> "1000100001", 6572 to 6574=> "1000100000", 6575 to 6578=> "1000011111", 
    6579 to 6582=> "1000011110", 6583 to 6585=> "1000011101", 6586 to 6589=> "1000011100", 6590 to 6593=> "1000011011", 
    6594 to 6596=> "1000011010", 6597 to 6600=> "1000011001", 6601 to 6604=> "1000011000", 6605 to 6607=> "1000010111", 
    6608 to 6611=> "1000010110", 6612 to 6615=> "1000010101", 6616 to 6618=> "1000010100", 6619 to 6622=> "1000010011", 
    6623 to 6626=> "1000010010", 6627 to 6629=> "1000010001", 6630 to 6633=> "1000010000", 6634 to 6637=> "1000001111", 
    6638 to 6640=> "1000001110", 6641 to 6644=> "1000001101", 6645 to 6648=> "1000001100", 6649 to 6651=> "1000001011", 
    6652 to 6655=> "1000001010", 6656 to 6658=> "1000001001", 6659 to 6662=> "1000001000", 6663 to 6666=> "1000000111", 
    6667 to 6669=> "1000000110", 6670 to 6673=> "1000000101", 6674 to 6676=> "1000000100", 6677 to 6680=> "1000000011", 
    6681 to 6683=> "1000000010", 6684 to 6687=> "1000000001", 6688 to 6691=> "1000000000", 6692 to 6694=> "0111111111", 
    6695 to 6698=> "0111111110", 6699 to 6701=> "0111111101", 6702 to 6705=> "0111111100", 6706 to 6708=> "0111111011", 
    6709 to 6712=> "0111111010", 6713 to 6715=> "0111111001", 6716 to 6719=> "0111111000", 6720 to 6722=> "0111110111", 
    6723 to 6726=> "0111110110", 6727 to 6729=> "0111110101", 6730 to 6733=> "0111110100", 6734 to 6736=> "0111110011", 
    6737 to 6740=> "0111110010", 6741 to 6743=> "0111110001", 6744 to 6747=> "0111110000", 6748 to 6750=> "0111101111", 
    6751 to 6754=> "0111101110", 6755 to 6757=> "0111101101", 6758 to 6761=> "0111101100", 6762 to 6764=> "0111101011", 
    6765 to 6768=> "0111101010", 6769 to 6771=> "0111101001", 6772 to 6775=> "0111101000", 6776 to 6778=> "0111100111", 
    6779 to 6781=> "0111100110", 6782 to 6785=> "0111100101", 6786 to 6788=> "0111100100", 6789 to 6792=> "0111100011", 
    6793 to 6795=> "0111100010", 6796 to 6799=> "0111100001", 6800 to 6802=> "0111100000", 6803 to 6805=> "0111011111", 
    6806 to 6809=> "0111011110", 6810 to 6812=> "0111011101", 6813 to 6816=> "0111011100", 6817 to 6819=> "0111011011", 
    6820 to 6822=> "0111011010", 6823 to 6826=> "0111011001", 6827 to 6829=> "0111011000", 6830 to 6833=> "0111010111", 
    6834 to 6836=> "0111010110", 6837 to 6839=> "0111010101", 6840 to 6843=> "0111010100", 6844 to 6846=> "0111010011", 
    6847 to 6850=> "0111010010", 6851 to 6853=> "0111010001", 6854 to 6856=> "0111010000", 6857 to 6860=> "0111001111", 
    6861 to 6863=> "0111001110", 6864 to 6866=> "0111001101", 6867 to 6870=> "0111001100", 6871 to 6873=> "0111001011", 
    6874 to 6876=> "0111001010", 6877 to 6880=> "0111001001", 6881 to 6883=> "0111001000", 6884 to 6886=> "0111000111", 
    6887 to 6890=> "0111000110", 6891 to 6893=> "0111000101", 6894 to 6896=> "0111000100", 6897 to 6900=> "0111000011", 
    6901 to 6903=> "0111000010", 6904 to 6906=> "0111000001", 6907 to 6909=> "0111000000", 6910 to 6913=> "0110111111", 
    6914 to 6916=> "0110111110", 6917 to 6919=> "0110111101", 6920 to 6923=> "0110111100", 6924 to 6926=> "0110111011", 
    6927 to 6929=> "0110111010", 6930 to 6932=> "0110111001", 6933 to 6936=> "0110111000", 6937 to 6939=> "0110110111", 
    6940 to 6942=> "0110110110", 6943 to 6946=> "0110110101", 6947 to 6949=> "0110110100", 6950 to 6952=> "0110110011", 
    6953 to 6955=> "0110110010", 6956 to 6959=> "0110110001", 6960 to 6962=> "0110110000", 6963 to 6965=> "0110101111", 
    6966 to 6968=> "0110101110", 6969 to 6971=> "0110101101", 6972 to 6975=> "0110101100", 6976 to 6978=> "0110101011", 
    6979 to 6981=> "0110101010", 6982 to 6984=> "0110101001", 6985 to 6988=> "0110101000", 6989 to 6991=> "0110100111", 
    6992 to 6994=> "0110100110", 6995 to 6997=> "0110100101", 6998 to 7000=> "0110100100", 7001 to 7004=> "0110100011", 
    7005 to 7007=> "0110100010", 7008 to 7010=> "0110100001", 7011 to 7013=> "0110100000", 7014 to 7016=> "0110011111", 
    7017 to 7020=> "0110011110", 7021 to 7023=> "0110011101", 7024 to 7026=> "0110011100", 7027 to 7029=> "0110011011", 
    7030 to 7032=> "0110011010", 7033 to 7036=> "0110011001", 7037 to 7039=> "0110011000", 7040 to 7042=> "0110010111", 
    7043 to 7045=> "0110010110", 7046 to 7048=> "0110010101", 7049 to 7051=> "0110010100", 7052 to 7055=> "0110010011", 
    7056 to 7058=> "0110010010", 7059 to 7061=> "0110010001", 7062 to 7064=> "0110010000", 7065 to 7067=> "0110001111", 
    7068 to 7070=> "0110001110", 7071 to 7073=> "0110001101", 7074 to 7077=> "0110001100", 7078 to 7080=> "0110001011", 
    7081 to 7083=> "0110001010", 7084 to 7086=> "0110001001", 7087 to 7089=> "0110001000", 7090 to 7092=> "0110000111", 
    7093 to 7095=> "0110000110", 7096 to 7098=> "0110000101", 7099 to 7102=> "0110000100", 7103 to 7105=> "0110000011", 
    7106 to 7108=> "0110000010", 7109 to 7111=> "0110000001", 7112 to 7114=> "0110000000", 7115 to 7117=> "0101111111", 
    7118 to 7120=> "0101111110", 7121 to 7123=> "0101111101", 7124 to 7126=> "0101111100", 7127 to 7129=> "0101111011", 
    7130 to 7133=> "0101111010", 7134 to 7136=> "0101111001", 7137 to 7139=> "0101111000", 7140 to 7142=> "0101110111", 
    7143 to 7145=> "0101110110", 7146 to 7148=> "0101110101", 7149 to 7151=> "0101110100", 7152 to 7154=> "0101110011", 
    7155 to 7157=> "0101110010", 7158 to 7160=> "0101110001", 7161 to 7163=> "0101110000", 7164 to 7166=> "0101101111", 
    7167 to 7169=> "0101101110", 7170 to 7172=> "0101101101", 7173 to 7176=> "0101101100", 7177 to 7179=> "0101101011", 
    7180 to 7182=> "0101101010", 7183 to 7185=> "0101101001", 7186 to 7188=> "0101101000", 7189 to 7191=> "0101100111", 
    7192 to 7194=> "0101100110", 7195 to 7197=> "0101100101", 7198 to 7200=> "0101100100", 7201 to 7203=> "0101100011", 
    7204 to 7206=> "0101100010", 7207 to 7209=> "0101100001", 7210 to 7212=> "0101100000", 7213 to 7215=> "0101011111", 
    7216 to 7218=> "0101011110", 7219 to 7221=> "0101011101", 7222 to 7224=> "0101011100", 7225 to 7227=> "0101011011", 
    7228 to 7230=> "0101011010", 7231 to 7233=> "0101011001", 7234 to 7236=> "0101011000", 7237 to 7239=> "0101010111", 
    7240 to 7242=> "0101010110", 7243 to 7245=> "0101010101", 7246 to 7248=> "0101010100", 7249 to 7251=> "0101010011", 
    7252 to 7254=> "0101010010", 7255 to 7257=> "0101010001", 7258 to 7260=> "0101010000", 7261 to 7263=> "0101001111", 
    7264 to 7266=> "0101001110", 7267 to 7269=> "0101001101", 7270 to 7272=> "0101001100", 7273 to 7275=> "0101001011", 
    7276 to 7278=> "0101001010", 7279 to 7281=> "0101001001", 7282 to 7284=> "0101001000", 7285 to 7287=> "0101000111", 
    7288 to 7290=> "0101000110", 7291 to 7293=> "0101000101", 7294 to 7296=> "0101000100", 7297 to 7299=> "0101000011", 
    7300 to 7302=> "0101000010", 7303 to 7305=> "0101000001", 7306 to 7308=> "0101000000", 7309 to 7311=> "0100111111", 
    7312 to 7314=> "0100111110", 7315 to 7316=> "0100111101", 7317 to 7319=> "0100111100", 7320 to 7322=> "0100111011", 
    7323 to 7325=> "0100111010", 7326 to 7328=> "0100111001", 7329 to 7331=> "0100111000", 7332 to 7334=> "0100110111", 
    7335 to 7337=> "0100110110", 7338 to 7340=> "0100110101", 7341 to 7343=> "0100110100", 7344 to 7346=> "0100110011", 
    7347 to 7349=> "0100110010", 7350 to 7352=> "0100110001", 7353 to 7355=> "0100110000", 7356 to 7358=> "0100101111", 
    7359 to 7361=> "0100101110", 7362 to 7363=> "0100101101", 7364 to 7366=> "0100101100", 7367 to 7369=> "0100101011", 
    7370 to 7372=> "0100101010", 7373 to 7375=> "0100101001", 7376 to 7378=> "0100101000", 7379 to 7381=> "0100100111", 
    7382 to 7384=> "0100100110", 7385 to 7387=> "0100100101", 7388 to 7390=> "0100100100", 7391 to 7393=> "0100100011", 
    7394 to 7395=> "0100100010", 7396 to 7398=> "0100100001", 7399 to 7401=> "0100100000", 7402 to 7404=> "0100011111", 
    7405 to 7407=> "0100011110", 7408 to 7410=> "0100011101", 7411 to 7413=> "0100011100", 7414 to 7416=> "0100011011", 
    7417 to 7419=> "0100011010", 7420 to 7421=> "0100011001", 7422 to 7424=> "0100011000", 7425 to 7427=> "0100010111", 
    7428 to 7430=> "0100010110", 7431 to 7433=> "0100010101", 7434 to 7436=> "0100010100", 7437 to 7439=> "0100010011", 
    7440 to 7442=> "0100010010", 7443 to 7444=> "0100010001", 7445 to 7447=> "0100010000", 7448 to 7450=> "0100001111", 
    7451 to 7453=> "0100001110", 7454 to 7456=> "0100001101", 7457 to 7459=> "0100001100", 7460 to 7462=> "0100001011", 
    7463 to 7465=> "0100001010", 7466 to 7467=> "0100001001", 7468 to 7470=> "0100001000", 7471 to 7473=> "0100000111", 
    7474 to 7476=> "0100000110", 7477 to 7479=> "0100000101", 7480 to 7482=> "0100000100", 7483 to 7485=> "0100000011", 
    7486 to 7487=> "0100000010", 7488 to 7490=> "0100000001", 7491 to 7493=> "0100000000", 7494 to 7496=> "0011111111", 
    7497 to 7499=> "0011111110", 7500 to 7502=> "0011111101", 7503 to 7504=> "0011111100", 7505 to 7507=> "0011111011", 
    7508 to 7510=> "0011111010", 7511 to 7513=> "0011111001", 7514 to 7516=> "0011111000", 7517 to 7519=> "0011110111", 
    7520 to 7521=> "0011110110", 7522 to 7524=> "0011110101", 7525 to 7527=> "0011110100", 7528 to 7530=> "0011110011", 
    7531 to 7533=> "0011110010", 7534 to 7536=> "0011110001", 7537 to 7538=> "0011110000", 7539 to 7541=> "0011101111", 
    7542 to 7544=> "0011101110", 7545 to 7547=> "0011101101", 7548 to 7550=> "0011101100", 7551 to 7552=> "0011101011", 
    7553 to 7555=> "0011101010", 7556 to 7558=> "0011101001", 7559 to 7561=> "0011101000", 7562 to 7564=> "0011100111", 
    7565 to 7567=> "0011100110", 7568 to 7569=> "0011100101", 7570 to 7572=> "0011100100", 7573 to 7575=> "0011100011", 
    7576 to 7578=> "0011100010", 7579 to 7581=> "0011100001", 7582 to 7583=> "0011100000", 7584 to 7586=> "0011011111", 
    7587 to 7589=> "0011011110", 7590 to 7592=> "0011011101", 7593 to 7595=> "0011011100", 7596 to 7597=> "0011011011", 
    7598 to 7600=> "0011011010", 7601 to 7603=> "0011011001", 7604 to 7606=> "0011011000", 7607 to 7609=> "0011010111", 
    7610 to 7611=> "0011010110", 7612 to 7614=> "0011010101", 7615 to 7617=> "0011010100", 7618 to 7620=> "0011010011", 
    7621 to 7622=> "0011010010", 7623 to 7625=> "0011010001", 7626 to 7628=> "0011010000", 7629 to 7631=> "0011001111", 
    7632 to 7634=> "0011001110", 7635 to 7636=> "0011001101", 7637 to 7639=> "0011001100", 7640 to 7642=> "0011001011", 
    7643 to 7645=> "0011001010", 7646 to 7647=> "0011001001", 7648 to 7650=> "0011001000", 7651 to 7653=> "0011000111", 
    7654 to 7656=> "0011000110", 7657 to 7659=> "0011000101", 7660 to 7661=> "0011000100", 7662 to 7664=> "0011000011", 
    7665 to 7667=> "0011000010", 7668 to 7670=> "0011000001", 7671 to 7672=> "0011000000", 7673 to 7675=> "0010111111", 
    7676 to 7678=> "0010111110", 7679 to 7681=> "0010111101", 7682 to 7683=> "0010111100", 7684 to 7686=> "0010111011", 
    7687 to 7689=> "0010111010", 7690 to 7692=> "0010111001", 7693 to 7694=> "0010111000", 7695 to 7697=> "0010110111", 
    7698 to 7700=> "0010110110", 7701 to 7703=> "0010110101", 7704 to 7705=> "0010110100", 7706 to 7708=> "0010110011", 
    7709 to 7711=> "0010110010", 7712 to 7714=> "0010110001", 7715 to 7716=> "0010110000", 7717 to 7719=> "0010101111", 
    7720 to 7722=> "0010101110", 7723 to 7725=> "0010101101", 7726 to 7727=> "0010101100", 7728 to 7730=> "0010101011", 
    7731 to 7733=> "0010101010", 7734 to 7736=> "0010101001", 7737 to 7738=> "0010101000", 7739 to 7741=> "0010100111", 
    7742 to 7744=> "0010100110", 7745 to 7747=> "0010100101", 7748 to 7749=> "0010100100", 7750 to 7752=> "0010100011", 
    7753 to 7755=> "0010100010", 7756 to 7758=> "0010100001", 7759 to 7760=> "0010100000", 7761 to 7763=> "0010011111", 
    7764 to 7766=> "0010011110", 7767 to 7769=> "0010011101", 7770 to 7771=> "0010011100", 7772 to 7774=> "0010011011", 
    7775 to 7777=> "0010011010", 7778 to 7779=> "0010011001", 7780 to 7782=> "0010011000", 7783 to 7785=> "0010010111", 
    7786 to 7788=> "0010010110", 7789 to 7790=> "0010010101", 7791 to 7793=> "0010010100", 7794 to 7796=> "0010010011", 
    7797 to 7798=> "0010010010", 7799 to 7801=> "0010010001", 7802 to 7804=> "0010010000", 7805 to 7807=> "0010001111", 
    7808 to 7809=> "0010001110", 7810 to 7812=> "0010001101", 7813 to 7815=> "0010001100", 7816 to 7818=> "0010001011", 
    7819 to 7820=> "0010001010", 7821 to 7823=> "0010001001", 7824 to 7826=> "0010001000", 7827 to 7828=> "0010000111", 
    7829 to 7831=> "0010000110", 7832 to 7834=> "0010000101", 7835 to 7837=> "0010000100", 7838 to 7839=> "0010000011", 
    7840 to 7842=> "0010000010", 7843 to 7845=> "0010000001", 7846 to 7847=> "0010000000", 7848 to 7850=> "0001111111", 
    7851 to 7853=> "0001111110", 7854 to 7855=> "0001111101", 7856 to 7858=> "0001111100", 7859 to 7861=> "0001111011", 
    7862 to 7864=> "0001111010", 7865 to 7866=> "0001111001", 7867 to 7869=> "0001111000", 7870 to 7872=> "0001110111", 
    7873 to 7874=> "0001110110", 7875 to 7877=> "0001110101", 7878 to 7880=> "0001110100", 7881 to 7883=> "0001110011", 
    7884 to 7885=> "0001110010", 7886 to 7888=> "0001110001", 7889 to 7891=> "0001110000", 7892 to 7893=> "0001101111", 
    7894 to 7896=> "0001101110", 7897 to 7899=> "0001101101", 7900 to 7901=> "0001101100", 7902 to 7904=> "0001101011", 
    7905 to 7907=> "0001101010", 7908 to 7910=> "0001101001", 7911 to 7912=> "0001101000", 7913 to 7915=> "0001100111", 
    7916 to 7918=> "0001100110", 7919 to 7920=> "0001100101", 7921 to 7923=> "0001100100", 7924 to 7926=> "0001100011", 
    7927 to 7928=> "0001100010", 7929 to 7931=> "0001100001", 7932 to 7934=> "0001100000", 7935 to 7936=> "0001011111", 
    7937 to 7939=> "0001011110", 7940 to 7942=> "0001011101", 7943 to 7945=> "0001011100", 7946 to 7947=> "0001011011", 
    7948 to 7950=> "0001011010", 7951 to 7953=> "0001011001", 7954 to 7955=> "0001011000", 7956 to 7958=> "0001010111", 
    7959 to 7961=> "0001010110", 7962 to 7963=> "0001010101", 7964 to 7966=> "0001010100", 7967 to 7969=> "0001010011", 
    7970 to 7971=> "0001010010", 7972 to 7974=> "0001010001", 7975 to 7977=> "0001010000", 7978 to 7979=> "0001001111", 
    7980 to 7982=> "0001001110", 7983 to 7985=> "0001001101", 7986 to 7987=> "0001001100", 7988 to 7990=> "0001001011", 
    7991 to 7993=> "0001001010", 7994 to 7996=> "0001001001", 7997 to 7998=> "0001001000", 7999 to 8001=> "0001000111", 
    8002 to 8004=> "0001000110", 8005 to 8006=> "0001000101", 8007 to 8009=> "0001000100", 8010 to 8012=> "0001000011", 
    8013 to 8014=> "0001000010", 8015 to 8017=> "0001000001", 8018 to 8020=> "0001000000", 8021 to 8022=> "0000111111", 
    8023 to 8025=> "0000111110", 8026 to 8028=> "0000111101", 8029 to 8030=> "0000111100", 8031 to 8033=> "0000111011", 
    8034 to 8036=> "0000111010", 8037 to 8038=> "0000111001", 8039 to 8041=> "0000111000", 8042 to 8044=> "0000110111", 
    8045 to 8046=> "0000110110", 8047 to 8049=> "0000110101", 8050 to 8052=> "0000110100", 8053 to 8054=> "0000110011", 
    8055 to 8057=> "0000110010", 8058 to 8060=> "0000110001", 8061 to 8062=> "0000110000", 8063 to 8065=> "0000101111", 
    8066 to 8068=> "0000101110", 8069 to 8070=> "0000101101", 8071 to 8073=> "0000101100", 8074 to 8076=> "0000101011", 
    8077 to 8078=> "0000101010", 8079 to 8081=> "0000101001", 8082 to 8084=> "0000101000", 8085 to 8086=> "0000100111", 
    8087 to 8089=> "0000100110", 8090 to 8092=> "0000100101", 8093 to 8094=> "0000100100", 8095 to 8097=> "0000100011", 
    8098 to 8100=> "0000100010", 8101 to 8102=> "0000100001", 8103 to 8105=> "0000100000", 8106 to 8108=> "0000011111", 
    8109 to 8110=> "0000011110", 8111 to 8113=> "0000011101", 8114 to 8116=> "0000011100", 8117 to 8118=> "0000011011", 
    8119 to 8121=> "0000011010", 8122 to 8124=> "0000011001", 8125 to 8126=> "0000011000", 8127 to 8129=> "0000010111", 
    8130 to 8132=> "0000010110", 8133 to 8134=> "0000010101", 8135 to 8137=> "0000010100", 8138 to 8140=> "0000010011", 
    8141 to 8142=> "0000010010", 8143 to 8145=> "0000010001", 8146 to 8148=> "0000010000", 8149 to 8150=> "0000001111", 
    8151 to 8153=> "0000001110", 8154 to 8156=> "0000001101", 8157 to 8158=> "0000001100", 8159 to 8161=> "0000001011", 
    8162 to 8164=> "0000001010", 8165 to 8166=> "0000001001", 8167 to 8169=> "0000001000", 8170 to 8172=> "0000000111", 
    8173 to 8174=> "0000000110", 8175 to 8177=> "0000000101", 8178 to 8180=> "0000000100", 8181 to 8182=> "0000000011", 
    8183 to 8185=> "0000000010", 8186 to 8188=> "0000000001", 8189 to 8191=> "0000000000" );


attribute EQUIVALENT_REGISTER_REMOVAL : string;
begin 


memory_access_guard_0: process (addr0) 
begin
      addr0_tmp <= addr0;
--synthesis translate_off
      if (CONV_INTEGER(addr0) > mem_size-1) then
           addr0_tmp <= (others => '0');
      else 
           addr0_tmp <= addr0;
      end if;
--synthesis translate_on
end process;

p_rom_access: process (clk)  
begin 
    if (clk'event and clk = '1') then
        if (ce0 = '1') then 
            q0 <= mem(CONV_INTEGER(addr0_tmp)); 
        end if;
    end if;
end process;

end rtl;


Library IEEE;
use IEEE.std_logic_1164.all;

entity tanh_tanh_table2 is
    generic (
        DataWidth : INTEGER := 10;
        AddressRange : INTEGER := 8192;
        AddressWidth : INTEGER := 13);
    port (
        reset : IN STD_LOGIC;
        clk : IN STD_LOGIC;
        address0 : IN STD_LOGIC_VECTOR(AddressWidth - 1 DOWNTO 0);
        ce0 : IN STD_LOGIC;
        q0 : OUT STD_LOGIC_VECTOR(DataWidth - 1 DOWNTO 0));
end entity;

architecture arch of tanh_tanh_table2 is
    component tanh_tanh_table2_rom is
        port (
            clk : IN STD_LOGIC;
            addr0 : IN STD_LOGIC_VECTOR;
            ce0 : IN STD_LOGIC;
            q0 : OUT STD_LOGIC_VECTOR);
    end component;



begin
    tanh_tanh_table2_rom_U :  component tanh_tanh_table2_rom
    port map (
        clk => clk,
        addr0 => address0,
        ce0 => ce0,
        q0 => q0);

end architecture;

