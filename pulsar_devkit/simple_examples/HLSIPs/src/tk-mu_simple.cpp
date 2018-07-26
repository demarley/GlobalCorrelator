/*
Created: 
Last Updated:   17 February 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

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



//void tkmu_simple_hw(  TkObj_tkmu& in, PropTkObj_tkmu& out ){
PropTkObj_tkmu tkmu_simple_hw( TkObj_tkmu& in ){
    /* Hardware implementation of the track propagation */
    PropTkObj_tkmu out;               // propagated track

    // constants
    // Conversions between binary and floating point (using example file to derive)
    finvpt_t INVRINV_CONVERSION(76090E-11);
    feta_t ETA_CONVERSION(512);
    feta_t INVETA_CONVERSION(0.001953125);
    fphi_t PHI_CONVERSION(219037);
    fphi_t INVPHI_CONVERSION(0.0000045654);
    fz0_t Z_CONVERSION(18);
    fz0_t INVZ_CONVERSION(0.0556);


    feta_t m_boundary(1.1);           // barrel/endcap boundary
    feta_t m_unity(1.0);
    fphi_t M_PI_144(0.0218);          // used in phi propagation
    fphi_t m_phi_A(1.464);            // constant in phi extrapolation
    fphi_t m_phi_B(2.82832);          // cosh(1.7): constant in phi extrapolation

    fphi_t m_phi_add_conv(-0.0387851);  // conversion in phi, addition (sector-dependent)
    fphi_t m_phi_mult_conv(0.232711);   // conversion in phi, multiplication (sector-dependent)

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
    if (DEBUG) std::cout << " FIRMWARE : SinhEta calculation " << in.hwSinhEta << std::endl;
    feta_t sinhEta;
    feta_t absSinhEta;
    if (in.hwSinhEta<0){
        absSinhEta = (in.hwSinhEta+8192)*INVETA_CONVERSION;     // use ap_fixed<>
    }
    else{
        absSinhEta = in.hwSinhEta*INVETA_CONVERSION;
    }
    if (DEBUG) std::cout << " -- |sinheta| = " << absSinhEta << std::endl;
    arcsinh(absSinhEta, inhwEta);
    feta_t abshwEta = inhwEta;
    if (in.hwSinhEta<0) inhwEta*=-1;
    in.hwEta = inhwEta*ETA_CONVERSION;             // set input eta ap_int

    if (DEBUG) std::cout << " -- eta     = " << inhwEta << std::endl;

    // Z0 (1024 = 2^10; number of unsigned bits)
    if (DEBUG) std::cout << " FIRMWARE : z0 calculation " << std::endl;
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

    if (DEBUG) std::cout << " -- inz0 = " << in.hwZ0 << std::endl;
    if (DEBUG) std::cout << " -- z0   = " << inhwZ0 << std::endl;
    if (DEBUG) std::cout << " -- |z0| = " << absInhwZ0 << std::endl;

    // Phi0 (262144 = 2^18; number of unsigned bits)
    fphi_t inhwPhi;
    if (in.hwPhi<0) inhwPhi = -1*(in.hwPhi+262144)*INVPHI_CONVERSION;
    else inhwPhi = in.hwPhi*INVPHI_CONVERSION;
    // sector-dependent conversion
    inhwPhi = inhwPhi - m_phi_add_conv + (in.hwSector-1)*m_phi_mult_conv;

    std::cout << " -- phi  = " << inhwPhi << std::endl;


    // Rinv -> 1/pT (16384 = 2^14; number of unsigned bits)
    finvpt_t inhwRinv;
    if (in.hwRinv<0)
        inhwRinv = (in.hwRinv+16384)*76090E-11;  //INVRINV_CONVERSION;
    else
        inhwRinv = in.hwRinv*76090E-11; //*INVRINV_CONVERSION;

//    fphi_t inhwInvPt(inhwRinv*87.719297);
    fphi_t inhwInvPt(inhwRinv*87);
    if (in.hwQ==1) inhwInvPt*=-1;

    if (DEBUG) std::cout << " -- inrinv = " << in.hwRinv << std::endl;
    if (DEBUG) std::cout << " -- rinv   = " << inhwRinv << std::endl;
    if (DEBUG) std::cout << " -- invpt  = " << inhwInvPt << std::endl;

    // Do the calculations!
    if (DEBUG) std::cout << " FIRMWARE : Eta calculation " << std::endl;
    if (abshwEta < m_boundary){
        // barrel
        if (DEBUG) std::cout << " FIRMWARE : -- Barrel " << std::endl;
        dzCorrPhi = m_unity;            // convert 1.0 to eta_t
        etaProp   = m_boundary;         // 1.1;

        // 2, 1DLUTs: [z0/550] * [1/cosh(|eta|)]
        deta_LUT(absInhwZ0,deta);                             // only takes positive values
        if (inhwZ0<0)
            deta *= -1;

        invCosh(abshwEta,invCoshEta_EtaBarrel);               // LUT: 1/cosh(|eta|)
        deta *= invCoshEta_EtaBarrel;

        if (DEBUG) std::cout << " FIRMWARE :       deta       = " << deta << std::endl;
        if (DEBUG) std::cout << " FIRMWARE :       invCoshEta = " << invCoshEta_EtaBarrel << std::endl;
        if (DEBUG) std::cout << " FIRMWARE :       dzCorrPhi  = " << dzCorrPhi << std::endl;
    }
    else {
        // endcap
        if (DEBUG) std::cout << " FIRMWARE : -- Endcap " << std::endl;
        etaProp = abshwEta;

        // LUT: z0/850.
        if (inhwZ0<0){
            delta_LUT(absInhwZ0,delta);     // only takes positive values
            delta *= -1;
        }
        else{
            delta_LUT(inhwZ0,delta);
        }

        if (DEBUG) std::cout << " FIRMWARE :       delta   = " << delta << std::endl;

        // Only LUTs for
        // 1)  z0 / (850+z0)
        // 2)  z0 / (850-z0)
        // depending on sign(z0): multiply by '-1', if necessary, & choose the correct LUT!
        deta = 0;
        if (in.hwEta>0){
            dzCorrPhi = m_unity-delta;
            // Check z0 value, call the correct LUT!
            if (inhwZ0<0){
                delta_plus_LUT(absInhwZ0,deta);                                 // LUT: delta / (1-delta)
                deta *= -1;
            }
            else
                delta_minus_LUT(inhwZ0,deta);                                   // LUT: delta / (1-delta)
        }
        else{
            dzCorrPhi = m_unity+delta;
            if (inhwZ0<0){
                delta_minus_LUT(absInhwZ0,deta);                                // LUT: delta / (1+delta)
                deta *= -1;
            }
            else
                delta_plus_LUT(inhwZ0,deta);                                    // LUT: delta / (1+delta)
        }

        if (DEBUG) std::cout << " FIRMWARE :       deta      = " << deta << std::endl;

        feta_t tanhEta;
        tanh(inhwEta,tanhEta);  // handles the sign internally
        deta*=tanhEta;

        if (DEBUG) std::cout << " FIRMWARE :       deta      = " << deta << std::endl;
        if (DEBUG) std::cout << " FIRMWARE :       tanhEta   = " << tanhEta << std::endl;
        if (DEBUG) std::cout << " FIRMWARE :       dzCorrPhi = " << dzCorrPhi << std::endl;
    }


    // ** calculate the propagated eta ** //
    if (DEBUG) std::cout << " FIRMWARE : -- ETA calculation " << inhwEta + deta << std::endl;
    eta_t etaconv(ETA_CONVERSION);
    out.hwPropEta = (inhwEta + deta)*etaconv;


    // ** calculate the propagated phi ** //
    if (DEBUG) std::cout << " FIRMWARE : -- Phi calculation " << std::endl;

    invCosh(etaProp,invCoshEta_Phi);            // LUT: 1/cosh(x)

    // Include two constants used in calculation (1.464*cosh(1.7))
    fphi_t tmp_A = m_phi_A;
    fphi_t tmp_B = m_phi_B;

    tmp_A *= inhwInvPt;       // 1.464 * 1/pT
    tmp_B *= dzCorrPhi;       // cosh(1.7) * dzCorrPhi
    fphi_t tmp_val4   = tmp_A * tmp_B * invCoshEta_Phi;
    fphi_t outPropPhi = inhwPhi - tmp_val4 - M_PI_144;

    out.hwPropPhi = outPropPhi*PHI_CONVERSION;

    std::cout << " OUTPROPPHI " << outPropPhi << std::endl;

    // Print results to screen for debugging
    if (DEBUG) std::cout << " FIRMWARE :    invCoshEta = " << invCoshEta_Phi << std::endl;
    if (DEBUG) std::cout << " FIRMWARE :    1.464/pT   = " << tmp_A << std::endl;
    if (DEBUG) std::cout << " FIRMWARE :    cosh(1.7)*dzcorrphi = " << tmp_B << std::endl;
    if (DEBUG) std::cout << " FIRMWARE :    1.464*cosh(1.7)*dzcorrphi / (pT*cosh(etaProp)) = " << tmp_val4 << std::endl;
    if (DEBUG) std::cout << " FIRMWARE : hwPropPhi     = " << outPropPhi << std::endl;

    if (DEBUG) std::cout << " FIRMWARE : in.hwPhi      = " << in.hwPhi << std::endl;
    if (DEBUG) std::cout << " FIRMWARE : out.hwPropPhi = " << out.hwPropPhi << std::endl;

    return out;
}


// THE END
