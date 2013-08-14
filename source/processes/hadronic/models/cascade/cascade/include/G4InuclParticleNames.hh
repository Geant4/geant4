//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4InuclParticleNames.hh 72074 2013-07-06 06:28:53Z mkelsey $
//
// Defines enums to map G4InuclElementaryParticle type codes to human
// readable names.  Meant to replace similar local enums scattered through
// the code.
//
// 20101029  M. Kelsey -- Move antinucleons to 50-series, add deuteron
// 20111007  M. Kelsey -- Change photon code to 9, for use in initial states
// 20130508  D. Wright -- Add leptons and electroweak bosons
// 20130627  M. Kelsey -- Add functions to convert enum to strings for printing
// 20130702  M. Kelsey -- Add type classifiers ported from G4InuclElemPart

#ifndef G4INUCL_PARTICLE_NAMES_HH
#define G4INUCL_PARTICLE_NAMES_HH

#include "globals.hh"

namespace G4InuclParticleNames {
  enum Long { nuclei=0, proton=1, neutron=2,
	      pionPlus=3, pionMinus=5, pionZero=7, photon=9,
	      kaonPlus=11, kaonMinus=13, kaonZero=15, kaonZeroBar=17, 
	      lambda=21, sigmaPlus=23, sigmaZero=25, sigmaMinus=27, 
	      xiZero=29, xiMinus=31, omegaMinus=33, 
	      deuteron=41, triton=43, He3=45, alpha=47,
	      antiProton=51, antiNeutron=53,
	      antiDeuteron=61, antiTriton=63, antiHe3=65, antiAlpha=67,
	      diproton=111, unboundPN=112, dineutron=122,
              electronNu=-1, muonNu=-3, tauNu=-5,
              antiElectronNu=-7, antiMuonNu=-9, antiTauNu=-11,
              WMinus=-13, WPlus=-15, Zzero=-17,
              electron=-21, muonMinus=-23, tauMinus=-25,
              positron=-27, muonPlus=-29, tauPlus=-31};

  // NOTE:  "km" cannot be used as conflicts with "kilometers" unit!
  enum Short { nuc=nuclei, pro=proton, neu=neutron,
	       pip=pionPlus, pim=pionMinus, pi0=pionZero, gam=photon,
	       kpl=kaonPlus, kmi=kaonMinus, k0=kaonZero, k0b=kaonZeroBar,
	       lam=lambda, sp=sigmaPlus, s0=sigmaZero, sm=sigmaMinus,
	       xi0=xiZero, xim=xiMinus, om=omegaMinus, deu=deuteron,
	       ap=antiProton, an=antiNeutron, 
	       ade=antiDeuteron, atr=antiTriton, ahe=antiHe3, aal=antiAlpha,
	       pp=diproton, pn=unboundPN, nn=dineutron,
               enu=electronNu, mnu=muonNu, tnu=tauNu,
               aenu=antiElectronNu, amnu=antiMuonNu, atnu=antiTauNu,
               wm=WMinus, wp=WPlus, z0=Zzero,
               ele=electron, mum=muonMinus, tm=tauMinus,
               pos=positron, mup=muonPlus, tp=tauPlus};

  // Convert enum value to enum strings above for printing
  const char* nameLong(G4int ptype);
  const char* nameShort(G4int ptype);
  inline const char* name(G4int ptype) { return nameLong(ptype); }

  // Classify particle types to reduce if-blocks in client code
  inline G4bool isPhoton(G4int ityp) { return (ityp==photon); }

  inline G4bool isMuon(G4int ityp) { return (ityp==muonMinus||ityp==muonPlus); }

  inline G4bool isElectron(G4int ityp) { return (ityp==electron||ityp==positron); }

  inline G4bool isNeutrino(G4int ityp) {
    return (ityp==electronNu || ityp==muonNu || ityp==tauNu ||
	    ityp==antiElectronNu || ityp==antiMuonNu || ityp==antiTauNu);
  }
  
  inline G4bool pion(G4int ityp) {
    return (ityp==pionPlus || ityp==pionMinus || ityp==pionZero);
  }

  inline G4bool nucleon(G4int ityp) { return (ityp==proton || ityp==neutron); }

  inline G4bool antinucleon(G4int ityp) {
    return (ityp==antiProton || ityp==antiNeutron);
  }

  inline G4bool quasi_deutron(G4int ityp) { return (ityp > 100); }

  // Emulates G4PD::GetBaryonNumber(), rather than a simple bool
  inline G4int baryon(G4int ityp) {
    return ((ityp==pro || ityp==neu || ityp==lam || ityp==sp || ityp==s0 ||
	     ityp==sm  || ityp==xi0 || ityp==xim || ityp==om) ? 1 :
	    (ityp==deu || ityp==pp  || ityp==pn  || ityp==nn) ? 2 :
	    (ityp==ap  || ityp==an) ? -1 : 
	    (ityp==ade) ? -2 :
	    (ityp==atr || ityp==ahe) ? -3 :
	    (ityp==aal) ? -4 : 0);
  }

  inline G4bool antibaryon(G4int ityp) { return baryon(ityp) < 0; }

  inline G4bool hyperon(G4int ityp) {
    return (ityp==lam || ityp==sp || ityp==s0 || ityp==sm  || ityp==xi0 ||
	    ityp==xim || ityp==om);
  }
}

#endif	/* G4INUCL_PARTICLE_NAMES_HH */
