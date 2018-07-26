#ifndef DATAFORMATS_EXTENDED_H
#define DATAFORMATS_EXTENDED_H

#include "dataformats.h"

namespace {

float getMuonPtFloat(int hwPt)
{
  return hwPt * MUONPT_CONVERSION;
}

float getMuonEtaFloat(int hwEta)
{
  return (1- 2*std::bitset<9>(hwEta)[8]) * from_twos_complement<9>(hwEta) * MUONETA_CONVERSION;
}

float getMuonPhiFloat(int hwPhi)
{
  return normalizePhi(hwPhi * MUONPHI_CONVERSION);
}
 
float getTrackPtFloat(int hwPt)
{
  return hwPt * INVPT_CONVERSION;
}

float getTrackEtaFloat(int hwEta)
{
  return hwEta * INVETA_CONVERSION;
}

float getTrackPhiFloat(int hwPhi)
{
  return normalizePhi(hwPhi * INVPHI_CONVERSION);
}

float getTrackEtaFloatFromSinhEta(eta_t hwSinhEta)
{
  feta_t absSinhEta;
  if (hwSinhEta<0){
    absSinhEta = (hwSinhEta+8192)*INVETA_CONVERSION;
  }
  else{
    absSinhEta = hwSinhEta*INVETA_CONVERSION;
  }
  feta_t inhwEta = asinh(absSinhEta);
  if (hwSinhEta<0) inhwEta*=-1;

  return inhwEta;
}

float getTrackPhiFloatFromLocalPhi(phi_t hwPhi, sector_t hwSector)
{
  fphi_t inhwPhi;
  if (hwPhi<0)
    inhwPhi = -1*(hwPhi+262144)*INVPHI_CONVERSION;
  else
    inhwPhi = hwPhi*INVPHI_CONVERSION;
  
  // convert to a global phi value
  fphi_t hwPhiGlobal = inhwPhi + phiOffSetValues[hwSector-1];
  return hwPhiGlobal;
}

}
// -- Define structs for physics objects in software
struct SimTrack {
  int index;
  float pt;
  float eta;
  float phi;
  int q;
  // constructor
  SimTrack() :
    index(0),
    pt(0),
    eta(0),
    phi(0),
    q(0)
  {
  }
};

struct SwTrack {
  int index;
  float rinv;
  float pt;
  float sinhEta;
  float eta;
  float phi;
  float z0;
  int q;
  int valid;
  int BX;
  int sector;
  // constructor
  SwTrack() :
    index(0),
    rinv(0),
    pt(0),
    sinhEta(0),
    eta(0),
    phi(0),
    z0(0),
    q(0),
    valid(0),
    BX(0),
    sector(0)
  {
  }
};

struct SwPropTrack : public SwTrack {
  float propEta;
  float propPhi;
  // constructor
  SwPropTrack() : 
    SwTrack(),
    propEta(0),
    propPhi(0)
  {
  }
  // copy constructor
  SwPropTrack(const SwTrack& ref)
    {
      index = ref.index;
      rinv = ref.rinv;
      pt = ref.pt;
      sinhEta = ref.sinhEta;
      eta = ref.eta;
      phi = ref.phi;
      z0 = ref.z0;
      q = ref.q;
      valid = ref.valid;
      BX = ref.BX;
      propEta = 0;
      propPhi = 0;
      sector = ref.sector;
    }
};

struct SwMuon {
  int index;
  float pt;
  float eta;
  float phi;
  int q;
  int valid;
  int BX;
  // constructor
  SwMuon() : 
    index(0),
    pt(0),
    eta(0),
    phi(0),
    q(0),
    valid(0),
    BX(0)
  {
  }
};

struct SwTrackMuon 
{
  int index;
  float pt;
  float eta;
  float phi;
  int q;
  int valid;   // valid bit
  int BX;    // bunch crossing
  int quality;
  // constructor
  SwTrackMuon() : 
    index(0),
    pt(0),
    eta(0),
    phi(0),
    q(0),
    valid(0),
    BX(0),
    quality(0)
  {
  }
};

bool sameSwMuon(const SwMuon& lhs,const SwMuon& rhs)
{
  float eta1 = lhs.eta;
  float eta2 = rhs.eta;
  float phi1 = lhs.phi;
  float phi2 = rhs.phi;
  float dR2 = (eta1 - eta2)*(eta1 - eta2) + normalizePhi(phi1 - phi2)*normalizePhi(phi1 - phi2);
  return dR2 < 0.3*0.3;
}

bool sameSwTrack(const SwTrack& lhs,const SwTrack& rhs)
{
  float eta1 = lhs.eta;
  float eta2 = rhs.eta;
  float phi1 = lhs.phi;
  float phi2 = rhs.phi;
  float dR2 = (eta1 - eta2)*(eta1 - eta2) + normalizePhi(phi1 - phi2)*normalizePhi(phi1 - phi2);
  return dR2 < 0.3*0.3;
}

bool sameHwMuon(const HwMuon& lhs,const HwMuon& rhs)
{
  float eta1 = getMuonEtaFloat(lhs.hwEta);
  float eta2 = getMuonEtaFloat(rhs.hwEta);
  float phi1 = getMuonPhiFloat(lhs.hwPhi);
  float phi2 = getMuonPhiFloat(rhs.hwPhi);
  float dR2 = (eta1 - eta2)*(eta1 - eta2) + normalizePhi(phi1 - phi2)*normalizePhi(phi1 - phi2);
  return dR2 < 0.3*0.3;
}

bool sameHwTrack(const HwTrack& lhs,const HwTrack& rhs)
{
  float eta1 = getTrackEtaFloatFromSinhEta(lhs.hwSinhEta);
  float eta2 = getTrackEtaFloatFromSinhEta(rhs.hwSinhEta);
  float phi1 = getTrackPhiFloatFromLocalPhi(lhs.hwPhi, lhs.hwSector);
  float phi2 = getTrackPhiFloatFromLocalPhi(rhs.hwPhi, rhs.hwSector);
  float dR2 = (eta1 - eta2)*(eta1 - eta2) + normalizePhi(phi1 - phi2)*normalizePhi(phi1 - phi2);
  return dR2 < 0.3*0.3;
}

struct Event 
{
  int eventNumber;
  int BX;
  std::vector<SimTrack> simTracks;

  std::vector<SwTrack> swTracks;
  std::vector<SwPropTrack> swPropTracks;
  std::vector<SwMuon> swMuons;
  // non-propagated track and muons
  std::vector<SwTrackMuon> swTrackMuons;
  // propagated track and muons
  std::vector<SwTrackMuon> swPropTrackMuons;

  std::vector<HwTrack> hwTracks;
  std::vector<HwPropTrack> hwPropTracks;
  std::vector<HwMuon> hwMuons;
  // non-propagated track and muons
  std::vector<HwTrackMuon> hwTrackMuons;
  // propagated track and muons
  std::vector<HwTrackMuon> hwPropTrackMuons;

  Event(){
    simTracks.clear();

    swTracks.clear();
    swPropTracks.clear();
    swMuons.clear();
    swTrackMuons.clear();
    swPropTrackMuons.clear();

    hwTracks.clear();
    hwPropTracks.clear();
    hwMuons.clear();
    hwTrackMuons.clear();
    hwPropTrackMuons.clear();
  }
};


namespace {

std::ostream& operator << (std::ostream& os, const SimTrack& rhs)
{
  os << "pT: " << rhs.pt << " " 
     << "eta: " << rhs.eta << " "
     << "phi: " << rhs.phi << " " 
     << "Q: " << rhs.q;
  return os;
}

std::ostream& operator << (std::ostream& os, const SwTrack& rhs)
{
  os << "pT: " << rhs.pt << " " 
     /* << "rinv: " << rhs.rinv << " " */
     /* << "sinhEta: " << rhs.sinhEta << " " */
     << "eta: " << rhs.eta << " "
     << "phi: " << rhs.phi << " " 
     << "sector: " << rhs.sector << " "
     << "Q: " << rhs.q;
  return os;
}

std::ostream& operator << (std::ostream& os, const SwPropTrack& rhs)
{
  os << "pT: " << rhs.pt << " " 
     /* << "sinhEta: " << rhs.sinhEta << " " */
     /* << "eta: " << rhs.eta << " " */
     /* << "phi: " << rhs.phi << " "  */
     << "eta_prop: " << rhs.propEta << " " 
     << "phi_prop: " << rhs.propPhi << " "
     << "sector: " << rhs.sector << " "
     << "Q: " << rhs.q;
  return os;
}

std::ostream& operator << (std::ostream& os, const SwMuon& rhs)
{
  os << "pT: " << rhs.pt << " " 
     << "eta: " << rhs.eta << " "
     << "phi: " << rhs.phi << " " 
     << "Q: " << rhs.q;
  return os;
}

std::ostream& operator << (std::ostream& os, const SwTrackMuon& rhs)
{
  os << "pT: " << rhs.pt << " " 
     << "eta: " << rhs.eta << " "
     << "phi: " << rhs.phi << " " 
     << "Q: " << rhs.q;
  return os;
}

std::ostream& operator << (std::ostream& os, const HwTrack& rhs)
{
  os << "pT: " << getTrackPtFloat(rhs.hwPt) << " "
     << "eta: " << getTrackEtaFloat(rhs.hwEta) << " "
     << "phi: " << getTrackPhiFloat(rhs.hwPhiGlobal) << " "
     << "Z0: " << rhs.hwZ0*INVZ_CONVERSION << " " 
     << "Q: " << 1- 2*rhs.hwQ << " "
     << "HwValid: " << std::bitset<1>(rhs.hwValid) << " "
     << "BX: " << std::bitset<3>(rhs.hwBX);
   return os;
}

std::ostream& operator << (std::ostream& os, const HwPropTrack& rhs)
{
  os << "pT: " << getTrackPtFloat(rhs.hwPt) << " "
     << "eta: " << getTrackEtaFloat(rhs.hwPropEta) << " "
     << "phi: " << getTrackPhiFloat(rhs.hwPropPhi) << " "
     << "Z0: " << rhs.hwZ0*INVZ_CONVERSION << " " 
     << "Q: " << 1- 2*rhs.hwQ << " "
     << "HwValid: " << std::bitset<1>(rhs.hwValid) << " "
     << "BX: " << std::bitset<3>(rhs.hwBX);
  return os;
}

std::ostream& operator << (std::ostream& os, const HwMuon& rhs)
{
  os << "pT: " << getMuonPtFloat(rhs.hwPt) << " " 
     << "eta: " << getMuonEtaFloat(rhs.hwEta) << " "
     << "phi: " << getMuonPhiFloat(rhs.hwPhi) << " "
     << "Q: " << 1- 2*rhs.hwQ << " "
     << "HwValid: " << std::bitset<1>(rhs.hwValid) << " "
     << "BX: " << std::bitset<3>(rhs.hwBX);
   return os;
}

std::ostream& operator << (std::ostream& os, const HwTrackMuon& rhs)
{
  os << "pT: " << getTrackPtFloat(rhs.hwPt) << " "
     << "eta: " << getMuonEtaFloat(rhs.hwEta) << " "
     << "phi: " << getMuonPhiFloat(rhs.hwPhi) << " "
     << "Q: " << 1- 2*rhs.hwQ << " "
     << "HwValid: " << std::bitset<1>(rhs.hwValid) << " "
     << "BX: " << std::bitset<3>(rhs.hwBX);
   return os;
}

std::ostream& operator << (std::ostream& os, const Event& rhs)
{
  os << "Event: " << rhs.eventNumber << " " 
     << "BX: " << rhs.BX << " " << std::endl
     << "SimTracks: " << rhs.simTracks.size() << std::endl;

  for (unsigned i=0; i<rhs.simTracks.size(); ++i){
    os << "  " << rhs.simTracks[i] << std::endl;
  }
  os << "SwTracks: " << rhs.swTracks.size() << std::endl;
  for (unsigned i=0; i<rhs.swTracks.size(); ++i){
    os << "  " << rhs.swTracks[i] << std::endl;
  }
  os << "HwTracks: " << rhs.hwTracks.size() << std::endl;
  for (unsigned i=0; i<rhs.hwTracks.size(); ++i) {
    os << "  " << rhs.hwTracks[i] << std::endl;
  }
  os << "SwPropTracks: " << rhs.swPropTracks.size() << std::endl;
  for (unsigned i=0; i<rhs.swPropTracks.size(); ++i) {
    os << "  " << rhs.swPropTracks[i] << std::endl;
  }
  os << "HwPropTracks: " << rhs.hwPropTracks.size() << std::endl;
  for (unsigned i=0; i<rhs.hwPropTracks.size(); ++i) {
    os << "  " << rhs.hwPropTracks[i] << std::endl;
  }
  os << "SwMuons: " << rhs.swMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.swMuons.size(); ++i) { 
    os << "  " << rhs.swMuons[i] << std::endl;
  }
  os << "HwMuons: " << rhs.hwMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.hwMuons.size(); ++i) { 
    os << "  " << rhs.hwMuons[i] << std::endl;
  }
  os << "SwTrackMuons: " << rhs.swTrackMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.swTrackMuons.size(); ++i) { 
    os << "  " << rhs.swTrackMuons[i] << std::endl;
  }
  os << "HwTrackMuons: " << rhs.hwTrackMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.hwTrackMuons.size(); ++i) { 
    os << "  " << rhs.hwTrackMuons[i] << std::endl;
  } 
  os << "SwPropTrackMuons: " << rhs.swPropTrackMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.swPropTrackMuons.size(); ++i) { 
    os << "  " << rhs.swPropTrackMuons[i] << std::endl;
  }
  os << "HwPropTrackMuons: " << rhs.hwPropTrackMuons.size() << std::endl; 
  for (unsigned i=0; i<rhs.hwPropTrackMuons.size(); ++i) { 
    os << "  " << rhs.hwPropTrackMuons[i] << std::endl;
  } 
  return os;  
}

}

#endif
