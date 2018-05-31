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

// functions to decode data
void decode_sw_track_data(const std::string &word, SwTrack&);
void decode_hw_track_data(const std::string &word, HwTrack&);
void decode_sw_muon_data(const std::string &word, SwMuon&);
void decode_hw_muon_data(const std::string &word, HwMuon&);

// helpers
void split(const std::string &s, 
	   char delim, std::vector<std::string> &elems);
float rinv2pt(const float& rinv);
float sinhEta2eta(const float& sinhEta);
void isNegative( std::string& bit_value, bool& isNeg);
float normalizePhi(float phi);

// dump output
void writeOutput(std::ofstream&,		   
		 const Event&);

void eventReader(std::vector<Event>& events, 
		 const std::string& simTrackFile,
		 const std::string& swTrack,
		 const std::string& swMuon,
		 const std::string& hwTrack,
		 const std::string& hwMuon); 

int main() 
{
    /* Run the program -- 
        Reading in data and predicting new eta/phi values 
	in software and firmware 
        Data from text file provided by track trigger group 
	(see GitLab group)
    */
  
    std::vector<std::string> muonComment;
    muonComment.push_back("Empty");
    muonComment.push_back("MuonB");
    muonComment.push_back("MuonO");
    muonComment.push_back("MuonE");

    // truth information I/O (CMSSW SimTrack)
    std::vector<Event> events;
    eventReader(events, 
		"../../../../config/sim_track_data.dat",
		"../../../../config/sw_track_data.dat",
		"../../../../config/sw_muon_data.dat",
		"../../../../config/hw_track_data.dat",
		"../../../../config/hw_muon_data.dat");

    for (unsigned int iEvent = 0; iEvent < events.size(); ++iEvent){
      // propagate the tracks
      for (unsigned int iHwTrack = 0; iHwTrack < events[iEvent].hwTracks.size(); ++iHwTrack){
	HwTrack hwTrack = events[iEvent].hwTracks[iHwTrack];
	HwPropTrack prop_track_hw = tkmu_simple_hw(hwTrack);
	events[iEvent].hwPropTracks.push_back(prop_track_hw);
      } 
      // propagate the tracks
      for (unsigned int iSwTrack = 0; iSwTrack < events[iEvent].swTracks.size(); ++iSwTrack){
	SwTrack swTrack = events[iEvent].swTracks[iSwTrack];
	SwPropTrack prop_track_sw = tkmu_simple_ref(swTrack);
	events[iEvent].swPropTracks.push_back(prop_track_sw);
      } 

      for (unsigned int iHwMuon = 0; iHwMuon < events[iEvent].hwMuons.size(); ++iHwMuon){
	HwMuon hwMuon = events[iEvent].hwMuons[iHwMuon];
	// check all possible matches!
	for (unsigned int iHwPropTrack = 0; 
	     iHwPropTrack < events[iEvent].hwPropTracks.size(); ++iHwPropTrack){
	  HwPropTrack hwPropTrack = events[iEvent].hwPropTracks[iHwPropTrack];
	  HwTrackMuon hwTrackMuon = match_hw(hwPropTrack, hwMuon);

	  if (hwTrackMuon.hwValid)
	    events[iEvent].hwTrackMuons.push_back(hwTrackMuon);
	}
      } 
      for (unsigned int iSwMuon = 0; iSwMuon < events[iEvent].swMuons.size(); ++iSwMuon){
	SwMuon swMuon = events[iEvent].swMuons[iSwMuon];
	// check all possible matches!
	for (unsigned int iSwPropTrack = 0; 
	     iSwPropTrack < events[iEvent].swPropTracks.size(); ++iSwPropTrack){
	  SwPropTrack swPropTrack = events[iEvent].swPropTracks[iSwPropTrack];
	  SwTrackMuon swTrackMuon = match_sw(swPropTrack, swMuon);

	  if (swTrackMuon.valid)
	    events[iEvent].swTrackMuons.push_back(swTrackMuon);
	} 
      }
    }

    
    std::ofstream data_output;
    data_output.open("../../../../output_data.txt");
    for (unsigned int i = 0; i < events.size(); ++i){
      std::cout << events[i] << std::endl;
      writeOutput(data_output, events[i]);
    }
    data_output.close();
    
    std::cout << " Finished " << std::endl;
    
    return 0;
}

std::string keepOnlyDigits(std::string input)
{
  std::string target = "";
  for(std::string::size_type i = 0; i < input.size(); ++i) {
    if (std::isdigit(input[i])) target += input[i];
  }
  return target;
}

void eventReader(std::vector<Event>& events, 
		 const std::string& simTrackFile,
		 const std::string& swTrack,
		 const std::string& swMuon,
		 const std::string& hwTrack,
		 const std::string& hwMuon)
{
  events.clear();

  std::ifstream ifile(simTrackFile.c_str());
  
  if (!ifile) {
    std::cout << "File does not exist: " << simTrackFile << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(0);
  }
  
  // step 1: generate as many events as there are in the truth file
  // open the file and put the truth data into a vector
  std::string line("");
  if (ifile.is_open()){
    
    while (std::getline(ifile, line)) {
      std::string newString(line);

      // create a new event
      if (newString.substr(0,5) == "Event"){
	Event newEvent;
	std::vector<std::string> values;
	split(newString,' ',values);
	newEvent.eventNumber = std::atoi( values.at(1).c_str() );
	newEvent.BX = std::atoi( values.at(3).c_str() );	
	// if (newEvent.eventNumber > 1) break;
	events.push_back(newEvent);
      }
      
      if (newString.substr(0,8) == "SimTrack"){
	SimTrack newSimTrack;
	std::vector<std::string> values;
	split(newString,' ',values);
	newSimTrack.index =  std::atof( values.at(1).c_str() );
	newSimTrack.pt =  std::atof( values.at(2).c_str() );
	newSimTrack.eta = std::atof( values.at(3).c_str() );
	newSimTrack.phi = std::atof( values.at(4).c_str() );
	newSimTrack.q  =  std::atoi( values.at(5).c_str() );
	events.back().simTracks.push_back(newSimTrack);
      }
    }
  }

  ifile.close();

  // step 2: read the Sw tracks
  std::ifstream swTrackFile(swTrack.c_str());
  
  if (!swTrackFile) {
    std::cout << "File does not exist: " << swTrack << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(0);
  }
  
  // open the file and put the truth data into a vector
  line = "";
  if (swTrackFile.is_open()){
    
    std::string phisector;
    int eventNumber;
    int BX;
    while (std::getline(swTrackFile, line)) {
      std::string newString(line);
      // create a new event
      if (newString.substr(0,2) == "BX"){
	std::vector<std::string> values;
	split(newString,' ',values);
	eventNumber = std::atoi( values.at(3).c_str() );
	if (eventNumber > events.size()) break;
	phisector = values.at(7);
      }
      // add attributes
      else{
       	SwTrack newTrack;
	newTrack.BX = 0;
	// std::cout << "eventNumber "  << eventNumber <<  " " << newString + " " + phisector << std::endl;
	decode_sw_track_data(newString + " " + phisector, newTrack);
	if (newTrack.valid)
	  events[eventNumber-1].swTracks.push_back(newTrack);
      }
    }
  }
  swTrackFile.close();


  // step 3: read the Hw tracks
  std::ifstream hwTrackFile(hwTrack.c_str());
  
  if (!hwTrackFile) {
    std::cout << "File does not exist: " << hwTrack << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(0);
  }
  
  // open the file and put the truth data into a vector
  line = "";
  if (hwTrackFile.is_open()){
    
    std::string phisector;
    int eventNumber;
    int BX;
    while (std::getline(hwTrackFile, line)) {
      std::string newString(line);

      // create a new event
      if (newString.substr(0,2) == "BX"){
	std::vector<std::string> values;
	split(newString,' ',values);
	eventNumber = std::atoi( values.at(3).c_str() );
	if (eventNumber > events.size()) break;
	phisector = values.at(7);
      }
      // add attributes
      else{
	HwTrack newTrack;
	newTrack.hwBX = std::bitset<3>(eventNumber).to_ulong();
	decode_hw_track_data(newString + " " + phisector, newTrack);
	//if (newTrack.hwValid)
	  events[eventNumber-1].hwTracks.push_back(newTrack);
      }
    }
  }
  hwTrackFile.close();


  // step 4: read the Sw muons
  std::ifstream swMuonFile(swMuon.c_str());
  
  if (!swMuonFile) {
    std::cout << "File does not exist: " << swMuon << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(0);
  }
  
  // open the file and put the truth data into a vector
  line = "";
  if (swMuonFile.is_open()){
    
    int eventNumber;
    int BX;
    while (std::getline(swMuonFile, line)) {
      std::string newString(line);
      // create a new event
      if (newString.substr(0,5) == "Begin"){
	std::vector<std::string> values;
	split(newString,' ',values);
	eventNumber = std::atoi(keepOnlyDigits( values.at(3)).c_str() );
	if (eventNumber > events.size()) break;
      }
      // add attributes
      else{
       	SwMuon newMuon;
	newMuon.BX = eventNumber;
	decode_sw_muon_data(newString, newMuon);
	if (newMuon.valid)
	  events[eventNumber-1].swMuons.push_back(newMuon);
      }
    }
  }
  swMuonFile.close();


  // step 5: read the Hw muons
  std::ifstream hwMuonFile(hwMuon.c_str());
  
  if (!hwMuonFile) {
    std::cout << "File does not exist: " << hwMuon << std::endl;
    std::cout << "Exiting. " << std::endl;
    exit(0);
  }
  
  // open the file and put the truth data into a vector
  line = "";
  if (hwMuonFile.is_open()){
    
    int eventNumber;
    int BX;
    while (std::getline(hwMuonFile, line)) {
      std::string newString(line);
      // create a new event
      if (newString.substr(0,5) == "Begin"){
	std::vector<std::string> values;
	split(newString,' ',values);
	eventNumber = std::atoi(keepOnlyDigits( values.at(3)).c_str() );
	if (eventNumber > events.size()) break;
      }
      // add attributes
      else{
       	HwMuon newMuon;
	newMuon.hwBX = std::bitset<3>(eventNumber).to_ulong();
	decode_hw_muon_data(newString, newMuon);
	
	//if (newMuon.hwValid)
	events[eventNumber-1].hwMuons.push_back(newMuon);
      }
    }
  }
  hwMuonFile.close();
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


float rinv2pt(const float& rinv)
{
  /* Convert rinv to pT 
     pt = (0.3*3.8/100.0)/irinv
  */
  float pt(-999.);
  float r = (rinv>0) ? rinv : -rinv;
  pt = 0.0114 / r;
  return pt;
}


float sinhEta2eta(const float& sinhEta)
{
  /* Convert sinhEta to eta */
  float eta(-999.);
  eta = log( sinhEta + sqrt(1+pow(sinhEta,2)) );
  return eta;
}


void isNegative(std::string& bit_value, bool& isNeg)
{
  /* Convert string of signed bit value to unsigned */
  isNeg = false;
  const char bit_value_0 = bit_value[0];
  if (bit_value_0 == '1'){
    bit_value = bit_value.replace(0,1,"0");
    isNeg = true;
  }
  return;
}

void decode_sw_track_data(const std::string &data_sw, SwTrack& in_track_sw)
{
  std::vector<std::string> values_sw;
  split(data_sw,' ',values_sw);
  
  // setup Rinv
  float rinv    = std::atof(values_sw.at(1).c_str());
  float sinhEta = std::atof(values_sw.at(3).c_str());

  in_track_sw.rinv = rinv;
  in_track_sw.pt  = rinv2pt(rinv);           // 1.360636778;
  in_track_sw.sinheta = sinhEta;
  in_track_sw.eta = sinhEta2eta(sinhEta);    //-2.265;
  in_track_sw.phi = std::atof(values_sw.at(2).c_str()); // 1.26735;
  in_track_sw.z0  = std::atof(values_sw.at(4).c_str()); //-4.72;
  in_track_sw.q   = (rinv>0) ? 1 : -1;       //-1;
  in_track_sw.valid = in_track_sw.pt > 0;
  in_track_sw.sector = std::atoi(values_sw.at(11).c_str());
}

void decode_hw_track_data(const std::string &data_fw, HwTrack& in_track_hw)
{
  std::vector<std::string> values_fw;
  split(data_fw,' ',values_fw);

  // setup Rinv
  bool isNegativeCharge(false);
  std::string rinv_str = values_fw.at(1);
  isNegative(rinv_str, isNegativeCharge);
  in_track_hw.hwRinv = std::bitset<15>(rinv_str).to_ulong();   
  //std::stoi(values_fw.at(1).c_str(),NULL,2);
  if (isNegativeCharge) in_track_hw.hwRinv *= -1;
  in_track_hw.hwQ = (isNegativeCharge) ? -1 : 1;

  // setup sinhEta
  bool isNegativeEta(false);
  std::string eta_str = values_fw.at(3);
  isNegative(eta_str, isNegativeEta);
  in_track_hw.hwSinhEta = std::bitset<14>(eta_str).to_ulong();   //std::stoi(values_fw.at(3).c_str(),NULL,2);
  if (isNegativeEta) in_track_hw.hwSinhEta *= -1;
  
  if (DEBUG) std::cout << " Looping over data - hwSinhEta = " << in_track_hw.hwSinhEta << std::endl;
  
  // setup phi0
  bool isNegativePhi(false);
  std::string phi_str = values_fw.at(2);
  isNegative(phi_str, isNegativePhi);
  in_track_hw.hwPhi = std::bitset<19>(phi_str).to_ulong();   //std::stoi(values_fw.at(2).c_str(),NULL,2);
  if (isNegativePhi) in_track_hw.hwPhi *= -1;
  
  if (DEBUG) std::cout << " Looping over data - hwPhi = " << in_track_hw.hwPhi << std::endl;
  
  // setup z0
  bool isNegativeZ0(false);
  std::string z0_str = values_fw.at(4);
  isNegative(z0_str, isNegativeZ0);
  in_track_hw.hwZ0 = std::bitset<12>(z0_str).to_ulong();   //std::stoi(values_fw.at(4).c_str(),NULL,2);
  if (isNegativeZ0) in_track_hw.hwZ0 *= -1;
  
  if (DEBUG) std::cout << " Looping over data - hwZ0 = " << in_track_hw.hwZ0 << std::endl;

  in_track_hw.hwValid = std::bitset<1>(rinv_str != "000000000000000").to_ulong();
}

void decode_sw_muon_data(const std::string &data_sw, SwMuon& in_muon_sw)
{
  std::vector<std::string> values_sw;
  split(data_sw,' ',values_sw);

  // note that even empty muons must have some kind of data in them, i.e. 
  // all members set to 0!
  in_muon_sw.pt  = std::atof(values_sw.at(1).c_str());
  in_muon_sw.eta = std::atof(values_sw.at(2).c_str());
  // keep in mind that muon phi is from 0 to 2pi
  in_muon_sw.phi = normalizePhi(std::atof(values_sw.at(3).c_str()));
  in_muon_sw.q   = std::atoi(values_sw.at(4).c_str());
  in_muon_sw.valid = int(in_muon_sw.pt > 0);
}

void decode_hw_muon_data(const std::string &data_fw, HwMuon& in_muon_hw)
{
  std::vector<std::string> values_fw;
  split(data_fw,' ',values_fw);

  // muon 1 or muon 2?
  std::string muonword = values_fw.at(0);
  int iMuon = std::atoi( &muonword.at(5) );

  std::string dataword = values_fw.at(1);

  // split up in 3 frames
  std::string frame1 = dataword.substr(64,32);
  std::string frame2 = dataword.substr(32,32);
  std::string frame3 = dataword.substr(0,32);

  // note that here, you must be careful with how the information is encoded
  /* now insert the bits into the datawords

   96-bit words look like

   first frame:
   BX0 phi H/F eta qual pT
   31  30  22  21  12   8  0
   
   second frame:
   SE track addresses VCH CH
   31 30               1  0

   third frame:
   B0  Empty
   31   0                                                                                                                                                                                                   

   fourth frame:                                                                                                                                              
   B1 phi H/F eta qual pT                                                                                                                                     
   31  30  22  21  12   8  0                                                                                                                                                                                
   
   fifth frame:                                                                                                                                             
   B2 track addresses VCH CH                                                                                                                          
   31 30               1  0                                                                                                                                                                                
   
   sixth frame:                                                                                                                               
   Empty                                                                                                                                                    
   31   0                    
  */
  

  // it is slightly difference for the first muon than for the second muon
  // int B0, B1, B2, BX0, SE;
  std::string pt_str, eta_str, phi_str, q_str;
  pt_str = frame2.substr(23,9);
  eta_str = frame2.substr(10,10);
  phi_str = frame2.substr(1,9);
  q_str = frame3.substr(0,1);

  in_muon_hw.hwPt = std::bitset<9>(pt_str).to_ulong();
  in_muon_hw.hwEta = std::bitset<10>(eta_str).to_ulong();
  in_muon_hw.hwPhi = std::bitset<9>(phi_str).to_ulong();
  in_muon_hw.hwQ = std::bitset<1>(q_str).to_ulong();
  //in_muon_hw.hwValid = in_muon_hw.hwEta != 0 and in_muon_hw.hwPhi != 0 and in_muon_hw.hwQ != 0;

  if (false) {
    std::cout << "data_fw " << data_fw << std::endl;
    std::cout << "dataword " << dataword << std::endl;
    std::cout << "frame 1: " << frame1 << std::endl;
    std::cout << "frame 2: " << frame2 << std::endl;
    std::cout << "frame 3: " << frame3 << std::endl;
    std::cout << "check muon data " 
	      << " pt " << pt_str << " " 
	      << " eta " << eta_str << " "
	      << " phi " << phi_str << " "
	      << " q " << q_str << std::endl;
  }
}

float normalizePhi(float outPhi)
{
  float returnValue = outPhi;
  if (returnValue <= -M_PI)
    returnValue += 2*M_PI;
  if (returnValue > M_PI)
    returnValue -= 2*M_PI;
  return returnValue;
}

template <class T> 
T matchingTrack(const SimTrack& simTrack, const std::vector<T> collection)
{
  T matchingT;
  // find the best matching T
  for (unsigned iT = 0; iT < collection.size(); iT++){
    // do eta based matching for now! 
    // relatively loose cut
    // phi and or charge may be misassigned 
    if (fabs(collection[iT].eta - simTrack.eta) < 0.2){
      matchingT = collection[iT];
      break;
    }
  }
  return matchingT;
}

void writeOutput(std::ofstream& output,		   
		   const Event& event)
{
  for (unsigned iSimTrack = 0; iSimTrack < event.simTracks.size(); iSimTrack++){
    output << "Event: " << event.eventNumber << ","
	   << " BX: " << event.BX << ","
	   << " pT: " << event.simTracks[iSimTrack].pt << "," 
	   << " eta: " << event.simTracks[iSimTrack].eta << ","
	   << " phi: " << event.simTracks[iSimTrack].phi << "," 
	   << " Q: " << event.simTracks[iSimTrack].q << ",";
    const SwTrack& matchingSwTrack =         matchingTrack(event.simTracks[iSimTrack], event.swTracks);
    const SwPropTrack& matchingSwPropTrack = matchingTrack(event.simTracks[iSimTrack], event.swPropTracks);
    const SwMuon& matchingSwMuon =           matchingTrack(event.simTracks[iSimTrack], event.swMuons);
    const SwTrackMuon& matchingSwTrackMuon = matchingTrack(event.simTracks[iSimTrack], event.swTrackMuons);

    //    output << matchingSwTrack << ",";
    output << matchingSwPropTrack << ",";
    output << matchingSwMuon << ",";
    output << matchingSwTrackMuon << ",";
    output << std::endl;
  }
}

// THE END
