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
// $Id: G4InuclParticleNames.hh,v 1.3 2010-06-25 09:43:30 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Defines enums to map G4InuclElementaryParticle type codes to human
// readable names.  Meant to replace similar local enums scattered through
// the code.

#ifndef G4INUCL_PARTICLE_NAMES_HH
#define G4INUCL_PARTICLE_NAMES_HH

namespace G4InuclParticleNames {
  enum Long { nuclei=0, proton=1, neutron=2,
	      pionPlus=3, pionMinus=5, pionZero=7, photon=10,
	      kaonPlus=11, kaonMinus=13, kaonZero=15, kaonZeroBar=17, 
	      lambda=21, sigmaPlus=23, sigmaZero=25, sigmaMinus=27, 
	      xiZero=29, xiMinus=31, omegaMinus=33, antiProton=35, 
	      antiNeutron=37, diproton=111, unboundPN=112, dineutron=122 };

  // NOTE:  "km" cannot be used as conflicts with "kilometers" unit!
  enum Short { nuc=nuclei, pro=proton, neu=neutron,
	       pip=pionPlus, pim=pionMinus, pi0=pionZero, gam=photon,
	       kpl=kaonPlus, kmi=kaonMinus, k0=kaonZero, k0b=kaonZeroBar,
	       lam=lambda, sp=sigmaPlus, s0=sigmaZero, sm=sigmaMinus,
	       xi0=xiZero, xim=xiMinus, om=omegaMinus, ap=antiProton,
	       an=antiNeutron, pp=diproton, pn=unboundPN, nn=dineutron };
}

#endif	/* G4INUCL_PARTICLE_NAMES_HH */
