/*
Propagate track to muon system

> Reference algorithm
  https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-inegration-CMSSW_9_3_2/L1Trigger/L1TTrackMatch/plugins/L1TkMuonProducer.cc#L271-L298
*/
#include <math.h>
#include <cmath>
#include <algorithm>
#include "src/tk-mu_simple.h"


SwPropTrack tkmu_simple_ref( const SwTrack& in ){

  SwPropTrack out(in);

   /* C++ implementation of the track propagation algorithm */
    float tk_pt   = in.pt;
    float tk_eta  = in.eta;
    float tk_aeta = std::abs(tk_eta);
    float tk_phi  = in.phi;
    int tk_q      = in.q;
    float tk_z    = in.z0;

    float deta = 0;
    float dzCorrPhi = 1.;
    float etaProp = tk_aeta;

    if (DEBUG) std::cout << " SOFTWARE : Eta calculation " << std::endl;
    if (tk_aeta < 1.1){
        if (DEBUG) std::cout << " SOFTWARE : -- barrel " << std::endl;
        etaProp = 1.1;
        deta    = tk_z/550./cosh(tk_aeta);
    }
    else {
        if (DEBUG) std::cout << " SOFTWARE : -- endcap " << std::endl;
        float delta = tk_z/850.; //roughly scales as distance to 2nd station
        if (tk_eta > 0) delta *=-1;
        dzCorrPhi = 1. + delta;

        float zOzs = tk_z/850.;
        if (tk_eta > 0) deta = zOzs/(1. - zOzs);
        else deta = zOzs/(1.+zOzs);
        deta = deta*tanh(tk_eta);
        if (DEBUG) std::cout << " SOFTWARE :       delta = " << delta << std::endl;
    }
    if (DEBUG) std::cout << " SOFTWARE :       deta = " << deta << std::endl;
    if (DEBUG) std::cout << " SOFTWARE :       dzcorrphi  = " << dzCorrPhi << std::endl;

    if (DEBUG) std::cout << " SOFTWARE : -- Phi calculation " << std::endl;
    if (DEBUG) std::cout << " SOFTWARE :       invCoshEta  = " << 1/cosh(etaProp) << std::endl;
    if (DEBUG) std::cout << " SOFTWARE :       1.464/invPt = " << 1.464/tk_pt << std::endl;
    if (DEBUG) std::cout << " SOFTWARE :       cosh(1.7)/cosh(etaProp) = " << 2.828/cosh(etaProp) << std::endl;
    // cosh(1.7)*1.464 = 4.14065
    float resPhi = tk_phi - tk_q*4.1407/cosh(etaProp)/tk_pt*dzCorrPhi - M_PI/144.;

    while (resPhi > M_PI)  resPhi -= 2.*M_PI;
    while (resPhi < -M_PI) resPhi += 2.*M_PI;

    // out.pt = in.pt;
    // out.eta = in.eta;
    // out.phi = in.phi;
    out.propEta = tk_eta + deta;
    out.propPhi = resPhi;

    return out;
}

SwTrackMuon match_sw(const SwTrack& inTrack, const SwMuon& inMuon)
{
  bool debug = false;

  SwTrackMuon outTrack;
  float tketa = inTrack.eta;
  float tkphi = inTrack.phi;
  
  float mueta = inMuon.eta;
  float muphi = inMuon.phi;

  // dR calculation
  float dR2_tk_mu = deltaR (tketa, tkphi, mueta, muphi);
  
  if (debug){
    std::cout 
      << "CheckMatch: tketa " << tketa
      << " tkphi " << tkphi
      << " mueta " << mueta
      << " muphi " << muphi
      << " dR2_tk_mu " << dR2_tk_mu
      << std::endl;
  }
  
  if (dR2_tk_mu < 0.2) {
    if (debug){
      std::cout << ">>>> MATCHED! <<<<" << std::endl;
    }
    outTrack.pt = inTrack.pt;
    outTrack.eta = inMuon.eta;
    outTrack.phi = inMuon.phi;
    outTrack.q = inTrack.q;
    outTrack.valid = 1; 
    outTrack.BX = inTrack.BX;
  } else {
    outTrack.valid = 0;     
    if (debug) std::cout << ">>>> NOT MATCHED! <<<<" << std::endl;
  }
  return outTrack;
}

SwTrackMuon match_prop_sw(const SwPropTrack& inTrack, const SwMuon& inMuon)
{
  bool debug = false;

  SwTrackMuon outTrack;
  float tketa = inTrack.propEta;
  float tkphi = inTrack.propPhi;
  
  float mueta = inMuon.eta;
  float muphi = inMuon.phi;

  // dR calculation
  float dR2_tk_mu = deltaR (tketa, tkphi, mueta, muphi);
  
  if (debug){
    std::cout 
      << "CheckMatch: tketa " << tketa
      << " tkphi " << tkphi
      << " mueta " << mueta
      << " muphi " << muphi
      << " dR2_tk_mu " << dR2_tk_mu
      << std::endl;
  }
  
  if (dR2_tk_mu < 0.2) {
    if (debug){
      std::cout << ">>>> MATCHED! <<<<" << std::endl;
    }
    outTrack.pt = inTrack.pt;
    outTrack.eta = inMuon.eta;
    outTrack.phi = inMuon.phi;
    outTrack.q = inTrack.q;
    outTrack.valid = 1; 
    outTrack.BX = inTrack.BX;
  } else {
    outTrack.valid = 0;     
    if (debug) std::cout << ">>>> NOT MATCHED! <<<<" << std::endl;
  }
  return outTrack;
}

float deltaR(float eta1, float phi1, float eta2, float phi2)
{
  float dEta = eta1 - eta1;
  float dPhi = normalizePhi(normalizePhi(phi1) - normalizePhi(phi2));
  return sqrt(dEta*dEta + dPhi*dPhi);
}
