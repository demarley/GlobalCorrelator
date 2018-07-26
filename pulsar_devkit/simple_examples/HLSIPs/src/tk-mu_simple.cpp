/*
Created: 
Last Updated:   17 February 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-1;95;0c----

HLS implementation of TK-MU Linking
  > https://github.com/cms-l1t-offline/cmssw/blob/phase2-l1t-inegration-CMSSW_9_3_2/L1Trigger/L1TTrackMatch/plugins/L1TkMuonProducer.cc

**
LUTs are simplified to only include positive calculations,
as long as negative calculations are simple transformations.

 pt = (0.3*3.8*0.01)/rinv

Negative values in binary are generated assuming "One's complement"
*/
#include "tk-mu_simple.h"
#include <cmath>
#include <cassert>
#ifndef __SYNTHESIS__
#include <cstdio>
#endif


etaphiglobal_t prop_hw(HwTrack& in)
{
  bool debug(false);

  /* Hardware implementation of the track propagation */
  feta_t boundary(1.1);           // barrel/endcap boundary
  feta_t unity(1.0);
  fphi_t M_PI_144(0.0218);        // used in phi propagation
  fphi_t tmp_A(1.464);            // constant in phi extrapolation
  fphi_t tmp_B(2.82832);          // cosh(1.7): constant in phi extrapolation

  // declare some variables
  fphi_t dzCorrPhi(1.0);
  fphi_t delta(0.0);

  feta_t inhwEta;
  feta_t deta(0.0);
  feta_t etaProp(1.1);

  fphi_t invCoshEta_Phi(0.0);
  feta_t invCoshEta_EtaBarrel(0.0);

  // ** convert inputs (ap_int<>) to ap_fixed<> for internal use ** //

  // sinhEta -> eta (8192 = 2^13; number of unsigned bits)
  feta_t sinhEta;
  feta_t absSinhEta;
  if (in.hwSinhEta<0){
    absSinhEta = (in.hwSinhEta+8192)*INVETA_CONVERSION;     // use ap_fixed<>
    //    sinhEta    = -1*(in.hwSinhEta+8192)*INVETA_CONVERSION;
  }
  else{
    absSinhEta = in.hwSinhEta*INVETA_CONVERSION;
    //sinhEta    = in.hwSinhEta*INVETA_CONVERSION;
  }
  arcsinh(absSinhEta, inhwEta);
  feta_t abshwEta = inhwEta;
  if (in.hwSinhEta<0) inhwEta*=-1;
  in.hwEta = inhwEta*ETA_CONVERSION;             // set input eta ap_int

  if (debug) {
    std::cout << " FIRMWARE : SinhEta calculation " << in.hwSinhEta << std::endl
	      << " -- |sinheta| = " << absSinhEta << std::endl
      // << " -- sinheta   = " << sinhEta << std::endl
	      << " -- eta       = " << inhwEta << std::endl
	      << " -- in.hwEta  = " << in.hwEta << std::endl;
  }
  // Z0 (1024 = 2^10; number of unsigned bits)
  fz0_t inhwZ0;
  fz0_t absInhwZ0;
  if (in.hwZ0<0) {
    inhwZ0    = -1*(in.hwZ0+1024)*INVZ_CONVERSION;
    absInhwZ0 = (in.hwZ0+1024)*INVZ_CONVERSION;
  }
  else{
    inhwZ0    = in.hwZ0*INVZ_CONVERSION;
    absInhwZ0 = in.hwZ0*INVZ_CONVERSION;
  }

  if (debug) {
     std::cout << " FIRMWARE : z0 calculation " << std::endl
	       << " -- inz0 = " << in.hwZ0 << std::endl
	       << " -- z0   = " << inhwZ0 << std::endl
	       << " -- |z0| = " << absInhwZ0 << std::endl;
  }

  // Phi0 (262144 = 2^18; number of unsigned bits)
  fphi_t inhwPhi;
  if (in.hwPhi<0)
    inhwPhi = -1*(in.hwPhi+262144)*INVPHI_CONVERSION;
  else
    inhwPhi = in.hwPhi*INVPHI_CONVERSION;
  
  // convert to a global phi value -- 1+22 bits (for 360 degrees)
  in.hwPhiGlobal = (inhwPhi + phiOffSetValues[in.hwSector-1]) * PHI_CONVERSION;
  if (debug) {
     std::cout << " FIRMWARE : phi calculation " << std::endl
	       << " -- local                     = " << inhwPhi << std::endl
	       << " -- sector                    = " << in.hwSector << std::endl
	       << " -- offset                    = " << phiOffSetValues[in.hwSector-1] << std::endl
	       << " -- global                    = " << inhwPhi + phiOffSetValues[in.hwSector-1] << std::endl
	       << " -- global (directly printed) = " << (inhwPhi + phiOffSetValues[in.hwSector-1]) * PHI_CONVERSION << std::endl
	       << " -- Phi Converted             = " << in.hwPhiGlobal << std::endl
	       << " -- InvPhi-Phi Converted      = " << in.hwPhiGlobal*INVPHI_CONVERSION << std::endl;
  }
  if (debug) {
    std::cout << " FIRMWARE : Charge calculation " << std::endl 
	      << " -- in.hwQ = " << in.hwQ << std::endl
	      << " -- Q      = " << 1- 2*in.hwQ << std::endl;
  }

  // Rinv (16384 = 2^14; number of unsigned bits)
  invpt_t absInvRinv;
  if (in.hwRinv < 0) 
    absInvRinv = (in.hwRinv+16384);
  else               
    absInvRinv = in.hwRinv;

  finvpt_t inhwRinv;
  inhwRinv = absInvRinv * INVRINV_CONVERSION;
  if (debug) {
    std::cout << " FIRMWARE : Rinv calculation " << std::endl
	      << " -- in.hwRinv  =       " << in.hwRinv << std::endl
	      << " -- absInvRinv =       " << absInvRinv << std::endl
	      << " -- inhwRinv (float) = " << inhwRinv << std::endl;
  }

  // calculate 1/pT with a LUT
  finvpt_t inhwInvPt;
  rinvToInvPt(inhwRinv, inhwInvPt);

  //finvpt_t absInhwRinv = inhwRinv;
  // hwQ==0 maps to +1; hwQ==1 maps to -1  
  if (in.hwQ > 0)
    inhwInvPt = -inhwInvPt;
  //   absInhwRinv = -inhwRinv;
    
  if (debug) {
    std::cout << " FIRMWARE : pT calculation " << std::endl
	      << " -- 1/pT           = " << inhwInvPt << std::endl;
	      // << " -- pT             = " << absPt << std::endl
	      // << " -- in.hwPt        = " << in.hwPt << std::endl
	      // << " -- pT (converted) = " << in.hwPt*INVPT_CONVERSION << std::endl;
  }

  // Do the calculations!
  if (debug) std::cout << " FIRMWARE : Eta calculation " << std::endl;
  if (abshwEta < boundary){
    // barrel
    if (debug) std::cout << " FIRMWARE : -- Barrel " << std::endl;
    dzCorrPhi = unity;            // convert 1.0 to eta_t
    etaProp   = boundary;         // 1.1;

    // 2, 1DLUTs: [z0/550] * [1/cosh(|eta|)]
    deta_LUT(absInhwZ0,deta);                             // only takes positive values
    if (inhwZ0<0)
      deta *= -1;

    invCosh(abshwEta,invCoshEta_EtaBarrel);               // LUT: 1/cosh(|eta|)
    deta *= invCoshEta_EtaBarrel;

    if (debug) std::cout << " FIRMWARE :       deta       = " << deta << std::endl;
    if (debug) std::cout << " FIRMWARE :       invCoshEta = " << invCoshEta_EtaBarrel << std::endl;
    if (debug) std::cout << " FIRMWARE :       dzCorrPhi  = " << dzCorrPhi << std::endl;
  }
  else {
    // endcap
    if (debug) std::cout << " FIRMWARE : -- Endcap " << std::endl;
    etaProp = abshwEta;

    // LUT: z0/850.
    if (inhwZ0<0){
      delta_LUT(absInhwZ0,delta);     // only takes positive values
      delta *= -1;
    }
    else{
      delta_LUT(inhwZ0,delta);
    }

    if (debug) std::cout << " FIRMWARE :       delta   = " << delta << std::endl;

    // Only LUTs for
    // 1)  z0 / (850+z0)
    // 2)  z0 / (850-z0)
    // depending on sign(z0): multiply by '-1', if necessary, & choose the correct LUT!
    deta = 0;
    if (in.hwEta>0){
      dzCorrPhi = unity-delta;
      // Check z0 value, call the correct LUT!
      if (inhwZ0<0){
	delta_plus_LUT(absInhwZ0,deta);                                 // LUT: delta / (1-delta)
	deta *= -1;
      }
      else
	delta_minus_LUT(inhwZ0,deta);                                   // LUT: delta / (1-delta)
    }
    else{
      dzCorrPhi = unity+delta;
      if (inhwZ0<0){
	delta_minus_LUT(absInhwZ0,deta);                                // LUT: delta / (1+delta)
	deta *= -1;
      }
      else
	delta_plus_LUT(inhwZ0,deta);                                    // LUT: delta / (1+delta)
    }

    if (debug) std::cout << " FIRMWARE :       deta      = " << deta << std::endl;

    feta_t tanhEta;
    tanh(inhwEta,tanhEta);  // handles the sign internally
    deta*=tanhEta;

    if (debug) std::cout << " FIRMWARE :       deta      = " << deta << std::endl;
    if (debug) std::cout << " FIRMWARE :       tanhEta   = " << tanhEta << std::endl;
    if (debug) std::cout << " FIRMWARE :       dzCorrPhi = " << dzCorrPhi << std::endl;
  }
  
  // ** calculate the propagated eta ** //
  if (debug) std::cout << " FIRMWARE : -- ETA calculation " << inhwEta + deta << std::endl;
  eta_t etaconv(ETA_CONVERSION);
  eta_t hwPropEta = (inhwEta + deta)*etaconv;


  // ** calculate the propagated phi ** //
  if (debug) std::cout << " FIRMWARE : -- Phi calculation " << std::endl;

  invCosh(etaProp,invCoshEta_Phi);            // LUT: 1/cosh(x)

  // Include two constants used in calculation (1.464*cosh(1.7))
  tmp_A *= inhwInvPt;       // 1.464 * 1/pT
  tmp_B *= dzCorrPhi;       // cosh(1.7) * dzCorrPhi
  fphiglobal_t tmp_val4   = tmp_A * tmp_B * invCoshEta_Phi;
  fphiglobal_t outPropPhi = fphiglobal_t(in.hwPhiGlobal*INVPHI_CONVERSION) - tmp_val4 - M_PI_144;

  phiglobal_t hwPropPhi = outPropPhi*PHI_CONVERSION;

  // Print results to screen for debugging
  if (debug) {
    std::cout << " FIRMWARE : in.hwPhiGlobal                                 = " << in.hwPhiGlobal*INVPHI_CONVERSION << std::endl;
    std::cout << " FIRMWARE : in.hwPhiGlobal (fphiglobal_t)                  = " << fphiglobal_t(in.hwPhiGlobal*INVPHI_CONVERSION) << std::endl;
    std::cout << " FIRMWARE : invCoshEta                                     = " << invCoshEta_Phi << std::endl;
    std::cout << " FIRMWARE : 1.464/pT                                       = " << tmp_A << std::endl;
    std::cout << " FIRMWARE : cosh(1.7)*dzcorrphi                            = " << tmp_B << std::endl;
    std::cout << " FIRMWARE : 1.464*cosh(1.7)*dzcorrphi / (pT*cosh(etaProp)) = " << tmp_val4 << std::endl;
    std::cout << " FIRMWARE : outPropPhi                                     = " << outPropPhi << std::endl;
    std::cout << " FIRMWARE : hwPropPhi*INVPHI                           = " << fphiglobal_t(hwPropPhi*INVPHI_CONVERSION) << std::endl;
  }
  return etaphiglobal_t(hwPropEta, hwPropPhi);
}

HwTrackMuon match_hw(HwTrack& inTrack, const HwMuon& inMuon)
{
  HwTrackMuon outTrackMuon;
  
  // assign the pT
  //  assign_pt_hw(inTrack);

  // calculate eta and phi
  feta_t tkEta = inTrack.hwEta*INVETA_CONVERSION;
  fphiglobal_t tkPhi = normalizePhi(inTrack.hwPhiGlobal*INVPHI_CONVERSION);
  
  feta_m muEta = (1- 2*std::bitset<9>(inMuon.hwEta)[8]) * from_twos_complement<9>(inMuon.hwEta) * MUONETA_CONVERSION;
  fphi_m muPhi = normalizePhi(inMuon.hwPhi * MUONPHI_CONVERSION);
  
  // dR calculation
  fphiglobal_t dR2_tk_mu = dr2_int(tkEta, tkPhi, muEta, muPhi);

  bool debug(false);
  if (debug){
    std::cout << "Track eta " << tkEta << std::endl;
    std::cout << "Track phi " << tkPhi << std::endl;
    std::cout << "muon eta " << muEta << std::endl;
    std::cout << "muon phi " << muPhi << std::endl;
    std::cout << "dR2 " << dR2_tk_mu << std::endl << std::endl;
  }

  // need to allow for negative values, since we do not take 
  // the square root
  if (dR2_tk_mu < 0.2*0.2 and dR2_tk_mu > - 0.2*0.2) {     
    if (debug) std::cout << ">>> Matched!" << std::endl << std::endl;
    
    outTrackMuon.hwPt = inTrack.hwPt;
    outTrackMuon.hwEta = inMuon.hwEta;
    outTrackMuon.hwPhi = inMuon.hwPhi;
    outTrackMuon.hwQ = inTrack.hwQ;
    outTrackMuon.hwValid = 1; 
    outTrackMuon.hwBX = inTrack.hwBX;
  } else {
    outTrackMuon.hwValid = 0;     
  }
  return outTrackMuon;
}

HwTrackMuon match_prop_hw(HwPropTrack& inTrack, const HwMuon& inMuon)
{
  HwTrackMuon outTrackMuon;
  
  // assign the pT
  assign_pt_hw(inTrack);

  // calculate eta and phi
  feta_t tkEta = inTrack.hwPropEta*INVETA_CONVERSION;
  fphiglobal_t tkPhi = normalizePhi(inTrack.hwPropPhi*INVPHI_CONVERSION);
  
  feta_m muEta = (1- 2*std::bitset<9>(inMuon.hwEta)[8]) * from_twos_complement<9>(inMuon.hwEta) * MUONETA_CONVERSION;
  fphi_m muPhi = normalizePhi(inMuon.hwPhi * MUONPHI_CONVERSION);

  // dR calculation
  fphiglobal_t dR2_tk_mu = dr2_int(tkEta, tkPhi, muEta, muPhi);

  bool debug(false);
  if (debug){
    std::cout << "Prop track eta " << tkEta << std::endl;
    std::cout << "Prop track phi " << tkPhi << std::endl;
    std::cout << "muon eta " << muEta << std::endl;
    std::cout << "muon phi " << muPhi << std::endl;
    std::cout << "dR2 " << dR2_tk_mu << std::endl << std::endl;
    }
  
  // need to allow for negative values, since we do not take 
  // the square root
  if (dR2_tk_mu < 0.2*0.2 and dR2_tk_mu > - 0.2*0.2) {     
    if (debug)     std::cout << ">>> Matched!" << std::endl << std::endl;

    outTrackMuon.hwPt = inTrack.hwPt;
    outTrackMuon.hwEta = inMuon.hwEta;
    outTrackMuon.hwPhi = inMuon.hwPhi;
    outTrackMuon.hwQ = inTrack.hwQ;
    outTrackMuon.hwValid = 1; 
    outTrackMuon.hwBX = inTrack.hwBX;
  } else {
    outTrackMuon.hwValid = 0;     
  }
  return outTrackMuon;
}

void assign_pt_hw(HwTrack& inTrack)
{
  // Rinv (16384 = 2^14; number of unsigned bits)
  invpt_t absInvRinv;
  if (inTrack.hwRinv < 0){ 
    absInvRinv = (inTrack.hwRinv+16384);
  } 
  else{               
    absInvRinv = inTrack.hwRinv;
  }
  finvpt_t inhwRinv;
  inhwRinv = absInvRinv * INVRINV_CONVERSION;
  
  // calculate pT from 1/pT with a LUT
  fpt_t absPt;
  rinvToPt(inhwRinv, absPt);
  inTrack.hwPt = absPt * fpt_t(PT_CONVERSION);
}

pt_t calc_pt_hw(invpt_t hwRinv)
{
  // Rinv (16384 = 2^14; number of unsigned bits)
  invpt_t absInvRinv;
  if (hwRinv < 0){ 
    absInvRinv = (hwRinv+16384);
  } 
  else{               
    absInvRinv = hwRinv;
  }
  finvpt_t inhwRinv;
  inhwRinv = absInvRinv * INVRINV_CONVERSION;
  
  // calculate pT from 1/pT with a LUT
  fpt_t absPt;
  rinvToPt(inhwRinv, absPt);

  return absPt * fpt_t(PT_CONVERSION);
}

eta_t calc_eta_hw(eta_t hwSinhEta)
{
  feta_t inhwEta;
  feta_t absSinhEta;
  if (hwSinhEta<0){
    absSinhEta = (hwSinhEta+8192)*INVETA_CONVERSION;     // use ap_fixed<>
  }
  else{
    absSinhEta = hwSinhEta*INVETA_CONVERSION;
  }
  arcsinh(absSinhEta, inhwEta);
  feta_t abshwEta = inhwEta;
  if (hwSinhEta<0) inhwEta*=-1;

  return inhwEta*ETA_CONVERSION;
}

phiglobal_t calc_phi_hw(phi_t hwPhi, sector_t hwSector)
{
  // Phi0 (262144 = 2^18; number of unsigned bits)
  fphi_t inhwPhi;
  if (hwPhi<0) {
    inhwPhi = -1*(hwPhi+262144)*INVPHI_CONVERSION;
  }
  else {
    inhwPhi = hwPhi*INVPHI_CONVERSION;
  }
  // convert to a global phi value -- 1+22 bits (for 360 degrees)
  return (inhwPhi + phiOffSetValues[hwSector-1]) * PHI_CONVERSION;
}

// THE END
