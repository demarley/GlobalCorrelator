#ifndef DATAFORMATS_H
#define DATAFORMATS_H

#include "ap_fixed.h"
#include "ap_int.h"
#include <stdlib.h>


typedef ap_int<15> invpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_int<12> pt_t;     // convert from RINV
typedef ap_int<14> eta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_int<19> phi_t;    // phi (50 micro-rad)
typedef ap_int<10> chisq_t;  // chi^2 (0 - 100; 0.1 steps)
typedef ap_int<1> q_t;       // charge
typedef ap_int<11> z0_t;     // z0  (1 mm over +/-14.9 cm)
typedef ap_int<3> bx_t;     // z0  (1 mm over +/-14.9 cm)
typedef ap_int<5> sector_t;

// before the decimal point, after the decimal point
typedef ap_fixed<15,2> finvpt_t;  // inverse pt [1% at 100 GeV]
typedef ap_fixed<12,9> fpt_t;     // 1/Rinv
typedef ap_fixed<14,4> feta_t;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<19,3> fphi_t;    // phi (50 micro-rad)
typedef ap_fixed<10,7> fchisq_t;  // chi^2 (0 - 100; 0.1 steps) 
typedef ap_fixed<11,5> fz0_t;     // z0  (1 mm over +/-14.9 cm) 

typedef ap_fixed<10,4> feta_m;    // eta [sinh(eta) measure to 0.005]
typedef ap_fixed<9,3> fphi_m;    // phi (50 micro-rad)

// muon data
typedef ap_int<9> pt_m;
typedef ap_int<10> eta_m;
typedef ap_int<9> phi_m;
typedef ap_int<4> quality_m;

// size of the LUTs
#define ETA_TABLE_SIZE 8192  // 13 unsigned bits
#define Z0_TABLE_SIZE 1024   // 10 unsigned bits

// range for LUTs
#define SINHETA_RANGE 6
#define ETA_RANGE 3
#define COSH_RANGE 3
#define Z0_RANGE 15
#define INV_SINHETA_RANGE 1/SINHETA_RANGE
#define INV_ETA_RANGE 1/ETA_RANGE
#define INV_COSH_RANGE 1/COSH_RANGE
#define INV_Z0_RANGE 1/Z0_RANGE

// Conversions between binary and floating point (using example file to derive)
#define RINV_CONVERSION 1314233             // 1/(76090E-11)
#define PT_CONVERSION 87719298E-6           // 1/(0.01*0.3*3.8); 87719298E-6
#define ETA_CONVERSION 855                  // 1/0.0011698 = 854.84698
#define PHI_CONVERSION 219037
#define Z_CONVERSION 18                     // 1/0.05615 = 17.81
#define INVRINV_CONVERSION 76090E-11
#define INVETA_CONVERSION 11698E-7
#define INVPHI_CONVERSION 456544E-11        // compare floating point and bit values from example file
#define INVZ_CONVERSION 5615E-5             // 0.05615 -- shift of 1024 for negative values!


// -- Define structs for physics objects in software
struct TrackObj_tkmu {
  float rinv;
  float pt;
  float sinheta;
  float eta;
  float phi;
  float z0;
  int q;
  int VALID;
  int BX;
  // constructor
  TrackObj_tkmu() :
    rinv(0),
    pt(0),
    sinheta(0),
    eta(0),
    phi(0),
    z0(0),
    q(0),
    VALID(0),
    BX(0)
  {
  }
};

struct MuonObj_tkmu {
  float pt;
  float eta;
  float phi;
  int q;
  int VALID;
  int BX;
  // constructor
  MuonObj_tkmu() : 
    pt(0),
    eta(0),
    phi(0),
    q(0),
    VALID(0),
    BX(0)
  {
  }
};

struct TrackMuonObj_tkmu 
{
  float pt;
  float eta;
  float phi;
  int q;
  int VALID;   // VALID bit
  int BX;    // bunch crossing
  // constructor
  TrackMuonObj_tkmu() : 
    pt(0),
    eta(0),
    phi(0),
    q(0),
    VALID(0),
    BX(0)
  {
  }
};


struct PropTrackObj_tkmu : public TrackObj_tkmu {
  float propEta;
  float propPhi;
  // constructor
  PropTrackObj_tkmu() : 
    TrackObj_tkmu(),
    propEta(0),
    propPhi(0)
  {
  }
};

// -- Define structs for physics objects in hardware
struct TkObj_tkmu 
{
  invpt_t hwRinv;
  pt_t hwPt;
  eta_t hwSinhEta;
  eta_t hwEta;
  phi_t hwPhi;
  z0_t hwZ0;  // same precision at eta_t (Accoring to Dan: z0 precision needed to get precision on eta, which derives from z0)
  // test precision of z0 11 to 14 
  q_t hwQ;
  chisq_t hwX2;
  q_t VALID;   // VALID bit
  bx_t hwBX;    // bunch crossing 3-bit counter
  /* sector_t hwSector; */
  // constructor
  TkObj_tkmu() : 
    hwRinv(0),
    hwPt(0),
    hwSinhEta(0),  
    hwEta(0),
    hwPhi(0),
    hwZ0(0),
    hwQ(0),
    hwX2(0),
    VALID(0),
    hwBX(0)
    /* hwSector(0) */
  {
  }
};

struct PropTkObj_tkmu : public TkObj_tkmu 
{
  eta_t hwPropEta;
  phi_t hwPropPhi;
  // constructor
 PropTkObj_tkmu() : 
  TkObj_tkmu(),
    hwPropEta(0),
    hwPropPhi(0)
      {
      }
  // copy constructor
  PropTkObj_tkmu(const TkObj_tkmu& ref)
    {
      hwRinv = ref.hwRinv;
      hwPt = ref.hwPt;
      hwSinhEta = ref.hwSinhEta;
      hwEta = ref.hwEta;
      hwPhi = ref.hwPhi;
      hwZ0 = ref.hwZ0;
      hwQ = ref.hwQ;
      hwX2 = ref.hwX2;
      VALID = ref.VALID;
      hwBX = ref.hwBX;
    }
};

struct MuObj_tkmu {
    pt_t hwPt;
    eta_t hwEta;
    phi_t hwPhi;
    q_t hwQ;
    q_t VALID;   // VALID bit
    bx_t hwBX;    // bunch crossing
  // constructor
  MuObj_tkmu() : 
    hwPt(0),
    hwEta(0),
    hwPhi(0),
    hwQ(0),
    VALID(0),
    hwBX(0)
  {
  }
};

struct TkMuObj_tkmu {
    pt_t hwPt;
    eta_m hwEta;
    phi_m hwPhi;
    q_t hwQ;
    q_t VALID;   // VALID bit
    bx_t hwBX;    // bunch crossing
  // constructor
  TkMuObj_tkmu() : 
    hwPt(0),
    hwEta(0),
    hwPhi(0),
    hwQ(0),
    VALID(0),
    hwBX(0)
  {
  }
};



namespace {

std::ostream& operator << (std::ostream& os, const TkObj_tkmu& rhs)
{
    os << "Rinv " << rhs.hwRinv << " " 
       << "pT " << rhs.hwPt << " " 
       << "eta " << rhs.hwEta << " "
       << "phi " << rhs.hwPhi << " " 
       << "Z0 " << rhs.hwZ0 << " " 
       << "Q " << rhs.hwQ << " "
       << "X2 " << rhs.hwX2 << " " 
       << "Valid " << rhs.VALID << " "
       << "BX " << rhs.hwBX << " ";
   return os;
}

std::ostream& operator << (std::ostream& os, const PropTkObj_tkmu& rhs)
{
    os << "Rinv " << rhs.hwRinv << " " 
       << "pT " << rhs.hwPt << " " 
       << "eta " << rhs.hwEta << " "
       << "phi " << rhs.hwPhi << " " 
       << "eta_prop " << rhs.hwPropEta << " " 
       << "phi_prop " << rhs.hwPropPhi << " "
       << "Z0 " << rhs.hwZ0 << " " 
       << "Q " << rhs.hwQ << " "
       << "X2 " << rhs.hwX2 << " " 
       << "Valid " << rhs.VALID << " "
       << "BX " << rhs.hwBX << " ";
   return os;
}

std::ostream& operator << (std::ostream& os, const PropTrackObj_tkmu& rhs)
{
    os << rhs.pt << " " 
       << rhs.eta << " "
       << rhs.phi << " " 
       << rhs.propEta << " " 
       << rhs.propPhi << " "
       << rhs.q << " ";
    return os;
}

std::ostream& operator << (std::ostream& os, const MuonObj_tkmu& rhs)
{
    os << rhs.pt << " " 
       << rhs.eta << " "
       << rhs.phi << " " 
       << rhs.q << " ";
    return os;
}

std::ostream& operator << (std::ostream& os, const TrackMuonObj_tkmu& rhs)
{
    os << rhs.pt << " " 
       << rhs.eta << " "
       << rhs.phi << " " 
       << rhs.q << " ";
    return os;
}

}

#endif
