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
// $Id$
// Geant4 tag: $Name:  $
//
// Defines enums to map G4InuclElementaryParticle type codes to human
// readable names.  Meant to replace similar local enums scattered through
// the code.
//
// 20101029  M. Kelsey -- Move antinucleons to 50-series, add deuteron
// 20111007  M. Kelsey -- Change photon code to 9, for use in initial states

#ifndef G4INUCL_PARTICLE_NAMES_HH
#define G4INUCL_PARTICLE_NAMES_HH

namespace G4InuclParticleNames {
  enum Long { nuclei=0, proton=1, neutron=2,
	      pionPlus=3, pionMinus=5, pionZero=7, photon=9,
	      kaonPlus=11, kaonMinus=13, kaonZero=15, kaonZeroBar=17, 
	      lambda=21, sigmaPlus=23, sigmaZero=25, sigmaMinus=27, 
	      xiZero=29, xiMinus=31, omegaMinus=33, 
	      deuteron=41, triton=43, He3=45, alpha=47,
	      antiProton=51, antiNeutron=53,
	      antiDeuteron=61, antiTriton=63, antiHe3=65, antiAlpha=67,
	      diproton=111, unboundPN=112, dineutron=122 };

  // NOTE:  "km" cannot be used as conflicts with "kilometers" unit!
  enum Short { nuc=nuclei, pro=proton, neu=neutron,
	       pip=pionPlus, pim=pionMinus, pi0=pionZero, gam=photon,
	       kpl=kaonPlus, kmi=kaonMinus, k0=kaonZero, k0b=kaonZeroBar,
	       lam=lambda, sp=sigmaPlus, s0=sigmaZero, sm=sigmaMinus,
	       xi0=xiZero, xim=xiMinus, om=omegaMinus, deu=deuteron,
	       ap=antiProton, an=antiNeutron, 
	       ade=antiDeuteron, atr=antiTriton, ahe=antiHe3, aal=antiAlpha,
	       pp=diproton, pn=unboundPN, nn=dineutron };
}

#endif	/* G4INUCL_PARTICLE_NAMES_HH */
