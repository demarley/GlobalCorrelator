/*
Propagate track to muon station

main script for running tests comparing
HLS with c++

    // :: FIRMWARE :: //
    // [rinv (15 bit) | phi0 (19 bit)     | t (14 bit)  | z0 (11 bit)]
    // 111100101101011|0000110001100000011|11111100000101|11110110000
*/
#include <cstdio>
#include <vector> 
#include <bitset>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "ap_fixed.h"
#include "src/tk-mu_simple.h"


void read_file(const std::string &file_name, std::vector<std::string> &values, const std::string &comment);
void split(const std::string &s, char delim, std::vector<std::string> &elems);
float rinv2pt(const float& rinv);
float sinhEta2eta(const float& sinhEta);
void isNegative( std::string& bit_value, bool& isNeg);

int main() {
    /* Run the program -- 
        Reading in data and predicting new eta/phi values in software and firmware 
        Data from text file provided by track trigger group (see GitLab group)
    */
    // software data I/O
    std::vector<std::string> real_data;
//    read_file("../../../../config/CleanTracks_muminus_2-10_real.dat",real_data,"BX");
    read_file("../../../../new_small_data_files/sw_track_data_0.dat",real_data,"BX");


    std::ofstream software_output;
    software_output.open("../../../../software_prop_7Jun2018.txt");

    // firmware data I/O
    std::vector<std::string> real_hw_data;
//    read_file("../../../../config/CleanTracks_muminus_2-10.dat",real_hw_data,"BX");
    read_file("../../../../new_small_data_files/hw_track_data_0.dat",real_hw_data,"BX");

    std::ofstream firmware_output;
    firmware_output.open("../../../../firmware_prop_7Jun2018.txt");


    // Conversions between binary and floating point (using example file to derive)




    for (unsigned int i=0,size=real_data.size();i<size;i++){
//    for (unsigned int i=1;i<2;i++){

        // :: SOFTWARE :: //
        std::string data_sw = real_data[i];
        std::vector<std::string> values_sw;
        split(data_sw,' ',values_sw);

        TrackObj_tkmu in_ref;

        float rinv    = std::atof(values_sw.at(1).c_str());
        float sinhEta = std::atof(values_sw.at(3).c_str());

        in_ref.pt  = rinv2pt(rinv);           // 1.360636778;
        in_ref.eta = sinhEta2eta(sinhEta);    //-2.265;
        in_ref.phi = std::atof(values_sw.at(2).c_str()); // 1.26735;
        in_ref.z0  = std::atof(values_sw.at(4).c_str()); //-4.72;
        in_ref.q   = (rinv>0) ? 0 : 1;       //-1;

        PropTrackObj_tkmu out_ref;

        tkmu_simple_ref(in_ref, out_ref);
        //std::cout << " REF : eta = " << in_ref.eta << " => " << out_ref.propEta << std::endl;
        //std::cout << "     : phi = " << in_ref.phi << " => " << out_ref.propPhi << std::endl;
        //std::cout << "     : in z0 = " << (in_ref.z0) << std::endl;
        //std::cout << "     : in invpT = " << (1/in_ref.pt) << std::endl;

        software_output << in_ref.eta << "," << out_ref.propEta << "," 
                        << in_ref.phi << "," << out_ref.propPhi << ","
                        << rinv       << "," << in_ref.z0       << "\n";


        // :: FIRMWARE :: //
        std::string data_fw = real_hw_data[i];
        std::vector<std::string> values_fw;
        split(data_fw,' ',values_fw);

        TkObj_tkmu in_hw;
        PropTkObj_tkmu out_hw;

        // setup Rinv
        bool isNegativeCharge(false);
        std::string rinv_str = values_fw.at(1);
        isNegative(rinv_str, isNegativeCharge);
        in_hw.hwRinv = std::bitset<15>(rinv_str).to_ulong();   //std::stoi(values_fw.at(1).c_str(),NULL,2);
        if (isNegativeCharge) in_hw.hwRinv *= -1;
        in_hw.hwQ = (isNegativeCharge) ? 1 : 0;

        if (DEBUG) std::cout << " Looping over data - hwRinv = " << in_hw.hwRinv << std::endl;

        // setup sinhEta
        bool isNegativeEta(false);
        std::string eta_str = values_fw.at(3);
        isNegative(eta_str, isNegativeEta);
        in_hw.hwSinhEta = std::bitset<14>(eta_str).to_ulong();   //std::stoi(values_fw.at(3).c_str(),NULL,2);
        if (isNegativeEta) in_hw.hwSinhEta *= -1;

        if (DEBUG) std::cout << " Looping over data - hwSinhEta = " << in_hw.hwSinhEta << std::endl;

        // setup phi0
        bool isNegativePhi(false);
        std::string phi_str = values_fw.at(2);
        isNegative(phi_str, isNegativePhi);
        in_hw.hwPhi = std::bitset<19>(phi_str).to_ulong();   //std::stoi(values_fw.at(2).c_str(),NULL,2);
        if (isNegativePhi) in_hw.hwPhi *= -1;

        if (DEBUG) std::cout << " Looping over data - hwPhi = " << in_hw.hwPhi << std::endl;

        // setup z0
        bool isNegativeZ0(false);
        std::string z0_str = values_fw.at(4);
        isNegative(z0_str, isNegativeZ0);
        in_hw.hwZ0 = std::bitset<12>(z0_str).to_ulong();   //std::stoi(values_fw.at(4).c_str(),NULL,2);
        if (isNegativeZ0) in_hw.hwZ0 *= -1;

        // get sector value
        unsigned int nEntries( values_fw.size() );
        in_hw.hwSector = std::atoi( values_fw.at(nEntries-1).c_str() );

        if (DEBUG) std::cout << " Looping over data - hwZ0 = " << in_hw.hwZ0 << std::endl;

        //in_hw.hwInvPt: convert in algorithm with '1/pt = 87.719298*in_hw.hwRinv;' // pt = (0.3*3.8*0.01)/rinv
        //in_hw.hwEta:   convert in algorithm with LUT

        out_hw = tkmu_simple_hw(in_hw);

        //std::cout << " HW  : eta   = " << float(in_hw.hwEta)/(ETA_CONVERSION) << " => " << float(out_hw.hwPropEta)/(ETA_CONVERSION) << std::endl;
        //std::cout << "     : phi   = " << float(in_hw.hwPhi)/(PHI_CONVERSION) << " => " << float(out_hw.hwPropPhi)/(PHI_CONVERSION) << std::endl;

        // save results to file
        fphi_t inhwPhi; 
        fphi_t outhwPhi;
        if (isNegativePhi){
            inhwPhi = -1*(in_hw.hwPhi+262144)*INVPHI_CONVERSION; 
            outhwPhi = -1*(out_hw.hwPropPhi+262144)*INVPHI_CONVERSION;
        }
        else{
            inhwPhi  = in_hw.hwPhi*INVPHI_CONVERSION;
            outhwPhi = out_hw.hwPropPhi*INVPHI_CONVERSION;
        }

        // sector-dependent conversion
        fphi_t m_phi_add_conv(-0.0387851);  // conversion in phi, addition (sector-dependent)
        fphi_t m_phi_mult_conv(0.232711);   // conversion in phi, multiplication (sector-dependent)
        inhwPhi = inhwPhi - m_phi_add_conv + (in_hw.hwSector-1)*m_phi_mult_conv;

         
        outhwPhi = outhwPhi - m_phi_add_conv + (in_hw.hwSector-1)*m_phi_mult_conv;

        finvpt_t inhwRinv;
        if (isNegativeCharge)
            inhwRinv = -1*(in_hw.hwRinv+16384)*INVRINV_CONVERSION;
        else
            inhwRinv = in_hw.hwRinv*INVRINV_CONVERSION;

        feta_t inhwZ0;
        if (isNegativeZ0)
            inhwZ0 = -1*(in_hw.hwZ0+1024)*INVZ_CONVERSION;
        else
            inhwZ0 = in_hw.hwZ0*INVZ_CONVERSION;

        firmware_output << in_hw.hwEta*INVETA_CONVERSION   << "," << out_hw.hwPropEta*INVETA_CONVERSION << ","
                        << inhwPhi   << "," << outhwPhi << ","
                        << inhwRinv << "," << inhwZ0 << "\n";
    }

    std::cout << " Finished " << std::endl;

    software_output.close();
    firmware_output.close();

    return 0;
}




void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    /* Split a string into pieces based on delimiter */
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return;
}


void read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment ) {
    /* Read in a generic file and put it into a vector of strings */
    std::ifstream ifile(file_name.c_str());

    if (!ifile) {
        std::cout << "File does not exist: " << file_name << std::endl;
        std::cout << "Exiting. " << std::endl;
        exit(0);
    }

    // open the file and put the data into a vector
    std::string line("");
    if (ifile.is_open()){

        std::string sector("");
        while (std::getline(ifile, line)) {
            std::string newstring(line);

            // remove all white spaces at the end of the string
            std::size_t space_pos = newstring.rfind(" ");
            while ( space_pos != std::string::npos && space_pos == newstring.size()-1 ) {
                newstring = newstring.substr(0, newstring.rfind(" "));
                space_pos = newstring.rfind(" ");
            }

            // BX= 000 event= 1 seed= CT_L1L2 phisector= 13
            // get the sector number
            std::size_t lineComment = newstring.find(comment);
            if (lineComment == 0){
                std::size_t total_size = newstring.size();
                sector = (newstring.at(total_size-2)!=' ') ? newstring.substr(total_size-2,2) : newstring.substr(total_size-1,1);
                std::cout << " SECTOR " << newstring << std::endl;
                continue;
            }

            // ignore empty lines
            if (newstring.length()==0) continue;
            if (sector.size()>0) newstring += " "+sector;

            values.push_back(newstring); // put values into vector
        }

        ifile.close();
    }

    return;
}


float rinv2pt(const float& rinv){
    /* Convert rinv to pT 
       pt = (0.3*3.8/100.0)/irinv
    */
    float pt(-999.);
    float r = (rinv>0) ? rinv : -rinv;
    pt = 0.0114 / r;
    return pt;
}


float sinhEta2eta(const float& sinhEta){
    /* Convert sinhEta to eta */
    float eta(-999.);
    eta = log( sinhEta + sqrt(1+pow(sinhEta,2)) );
    return eta;
}


void isNegative(std::string& bit_value, bool& isNeg){
    /* Convert string of signed bit value to unsigned */
    isNeg = false;
    const char bit_value_0 = bit_value[0];
    if (bit_value_0 == '1'){
        bit_value = bit_value.replace(0,1,"0");
        isNeg = true;
    }

    return;
}


// THE END
