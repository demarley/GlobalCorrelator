#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#include "ap_fixed.h"
#include "ap_int.h"
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <bitset>

// general stuff
const int N_BITS_GENERAL_Q = 1;
const int N_BITS_GENERAL_BX = 3;
const int N_BITS_GENERAL_B0 = 1;
const int N_BITS_GENERAL_SE = 1;
const int N_BITS_GENERAL_RES = 1;

// track specific
const int N_BITS_TRACK_INVPT = 15;
const int N_BITS_TRACK_PT = 14;
const int N_BITS_TRACK_ETA = 14;
const int N_BITS_TRACK_PHI = 19;
const int N_BITS_TRACK_PHIGLOBAL = 23;
const int N_BITS_TRACK_CHI2 = 10;
const int N_BITS_TRACK_Z0 = 11;
const int N_BITS_TRACK_SECTOR = 5;

// muon specific
const int N_BITS_MUON_PT = 9;
const int N_BITS_MUON_ETA = 9;
const int N_BITS_MUON_PHI = 10;
const int N_BITS_MUON_QUALITY = 4;

typedef ap_int<N_BITS_TRACK_INVPT> invpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_uint<N_BITS_TRACK_PT> pt_t;     // convert from RINV
typedef ap_int<N_BITS_TRACK_ETA> eta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_int<N_BITS_TRACK_PHI> phi_t;    // phi (50 micro-rad)
typedef ap_int<N_BITS_TRACK_PHIGLOBAL> phiglobal_t;    // phi (50 micro-rad)
typedef ap_int<N_BITS_TRACK_CHI2> chisq_t;  // chi^2 (0 - 100; 0.1 steps)
typedef ap_uint<1> q_t;       // charge
typedef ap_int<N_BITS_TRACK_Z0> z0_t;     // z0  (1 mm over +/-14.9 cm)
typedef ap_int<3> bx_t;     // z0  (1 mm over +/-14.9 cm)
typedef ap_uint<N_BITS_TRACK_SECTOR> sector_t;

// before the decimal point, after the decimal point
typedef ap_fixed<24,2> finvpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_fixed<14,9> fpt_t;     // 1/Rinv
typedef ap_fixed<14,4> feta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<19,3> fphi_t;    // phi (50 micro-rad)
typedef ap_fixed<23,3> fphiglobal_t;    // global phi (50 micro-rad)
typedef ap_fixed<10,7> fchisq_t;  // chi^2 (0 - 100; 0.1 steps) 
typedef ap_fixed<11,5> fz0_t;     // z0  (1 mm over +/-14.9 cm) 

typedef ap_fixed<10,4> feta_m;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<9,3> fphi_m;    // phi (50 micro-rad)

// muon data
typedef ap_uint<9> pt_m;
typedef ap_int<9> eta_m; // muon eta goes from -2.4 to 2.4 
typedef ap_uint<10> phi_m; // muon phi goes from 0 to 2pi
typedef ap_uint<4> quality_m;

// propetaphi
typedef std::pair<eta_t, phiglobal_t> etaphiglobal_t;
typedef std::pair<eta_m, phi_m> etaphiglobal_m;

// size of the LUTs
#define PT_TABLE_SIZE 16384  // 14 unsigned bits
#define RINV_TABLE_SIZE 16384 // 32768 // 16384  // 14 unsigned bits
#define ETA_TABLE_SIZE 8192  // 13 unsigned bits
#define Z0_TABLE_SIZE 1024   // 10 unsigned bits

// range for LUTs
#define PT_RANGE 175
#define INV_PT_RANGE 1/PT_RANGE

#define RINV_RANGE 57E-4  // (0.3*3.8*0.01)/2
#define INV_RINV_RANGE 175

#define SINHETA_RANGE 6
#define INV_SINHETA_RANGE 1/SINHETA_RANGE

#define ETA_RANGE 3
#define INV_ETA_RANGE 1/ETA_RANGE

#define COSH_RANGE 3
#define INV_COSH_RANGE 1/COSH_RANGE

#define Z0_RANGE 15
#define INV_Z0_RANGE 1/Z0_RANGE

// Conversions between binary and floating point (using example file to derive)
#define RINV_CONVERSION 792057              //1314229             // 1/0.000000760902077
#define INVRINV_CONVERSION 1262535E-12   //0.000001262538462  //original: 760902077E-15   
#define PT_CONVERSION 877192982456E-10  // 1/(0.01*0.3*3.8); 87719298E-6
#define INVPT_CONVERSION 114E-4

#define ETA_CONVERSION 512
#define INVETA_CONVERSION 19531261E-10

#define PHI_CONVERSION 211216
#define INVPHI_CONVERSION 47345E-10

#define Z_CONVERSION 17
#define INVZ_CONVERSION 5859375E-8

// Muon conversions
#define MUONPT_CONVERSION 0.5
#define MUONETA_CONVERSION 0.010875
#define MUONPHI_CONVERSION 0.010908

template<int N>
int from_twos_complement(int x) {
  return (x & (1<<(N-1))) ? (~(x-1) & ((1<<N)-1)) : x;
}

namespace {

template<class T>
T normalizePhi(T outPhi)
{
  T returnValue = outPhi;
  if (returnValue <= -T(M_PI))
    returnValue += T(2*M_PI);
  if (returnValue > T(M_PI))
    returnValue -= T(2*M_PI);
  return returnValue;
}

fphi_t phiOffSetValues[27] = {
  -0.0387851, // 1
  0.193925, // 2
  0.426636, // 3
  0.659347, // 4
  0.892057, // 5
  1.124768, // 6
  1.357478, // 7
  1.590189, // 8
  1.822899, // 9
  2.055610, // 10
  2.288321, // 11
  2.521031, // 12
  2.753742, // 13
  2.986452, // 14
  -3.064022, // 15
  -2.831312, // 16
  -2.598601, // 17
  -2.365891, // 18
  -2.133180, // 19
  -1.900470, // 20
  -1.667759, // 21
  -1.435049, // 22
  -1.202338, // 23
  -0.969627, // 24
  -0.736917, // 25
  -0.504206, // 26
  -0.271496  // 27
};

}
// -- Define structs for physics objects in hardware
struct HwTrack 
{
  invpt_t hwRinv;
  pt_t hwPt;
  eta_t hwSinhEta;
  eta_t hwEta;
  phi_t hwPhi;
  phiglobal_t hwPhiGlobal;
  z0_t hwZ0;  // same precision at eta_t (Accoring to Dan: z0 precision needed to get precision on eta, which derives from z0)
  sector_t hwSector;
  // test precision of z0 11 to 14 
  q_t hwQ;
  chisq_t hwX2;
  q_t hwValid;   // valid bit
  bx_t hwBX;    // bunch crossing 3-bit counter
  // constructor
  HwTrack() : 
    hwRinv(0),
    hwPt(0),
    hwSinhEta(0),  
    hwEta(0),
    hwPhi(0),
    hwPhiGlobal(0),
    hwZ0(0),
    hwSector(0),
    hwQ(0),
    hwX2(0),
    hwValid(0),
      hwBX(0)
  {
  }
};

struct HwPropTrack : public HwTrack 
{
  eta_t hwPropEta;
  phiglobal_t hwPropPhi;
  // constructor
 HwPropTrack() : 
  HwTrack(),
    hwPropEta(0),
    hwPropPhi(0)
      {
      }
  // copy constructor
  HwPropTrack(const HwTrack& ref)
    {
      hwRinv = ref.hwRinv;
      hwPt = ref.hwPt;
      hwSinhEta = ref.hwSinhEta;
      hwEta = ref.hwEta;
      hwPhi = ref.hwPhi;
      hwPhiGlobal = ref.hwPhiGlobal;
      hwZ0 = ref.hwZ0;
      hwSector = ref.hwSector;
      hwQ = ref.hwQ;
      hwX2 = ref.hwX2;
      hwValid = ref.hwValid;
      hwBX = ref.hwBX;
      hwPropEta = 0;
      hwPropPhi = 0;
    }
};

struct HwMuon {
    pt_m hwPt;
    eta_m hwEta;
    phi_m hwPhi;
    q_t hwQ;
    q_t hwValid;   // valid bit
    bx_t hwBX;    // bunch crossing
    quality_m hwQuality;
  // constructor
  HwMuon() : 
    hwPt(0),
    hwEta(0),
    hwPhi(0),
    hwQ(0),
    hwValid(0),
    hwBX(0),
    hwQuality(0)
  {
  }
};

// need a 64-bit trackmuon at least
struct HwTrackMuon {
  pt_t hwPt; //14
  eta_m hwEta; //9
  phi_m hwPhi; //9
  q_t hwQ; //1
  q_t hwValid;   //1 valid bit
  bx_t hwBX;    // 3 bunch crossing
  quality_m hwQuality;//4
  // constructor
  HwTrackMuon() : 
    hwPt(0),
    hwEta(0),
    hwPhi(0),
    hwQ(0),
    hwValid(0),
    hwBX(0),
    hwQuality(0)
  {
  }
};

#endif
