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
// $Id: G4CascadeChannel.cc,v 1.8 2010-06-25 09:44:02 gunter Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100514  M. Kelsey -- All functionality removed except quantum-number
//		validation functions.

#include "G4CascadeChannel.hh"
#include "G4ParticleDefinition.hh"
#include <vector>


std::vector<G4int> G4CascadeChannel::getQnums(G4int type) {
  G4int bary=0, str=0, ch=0;
  std::vector<G4int> Qnums(3);

    switch(type) {
    case 1: // proton
      bary = 1;
      str = 0;
      ch = 1;
      break;
    case 2: // neutron
      bary = 1;
      str = 0;
      ch = 0;
      break;
    case 3: // pi+
      bary = 0;
      str = 0;
      ch = 1;
      break;
    case 5: // pi-
      bary = 0;
      str = 0;
      ch = -1;
      break;
    case 7: // pi0
      bary = 0;
      str = 0;
      ch = 0;
      break;
    case 11: // k+
      bary = 0;
      str = 1;
      ch = 1;
      break;
    case 13: // k-
      bary = 0;
      str = -1;
      ch = -1;
      break;
    case 15: // k0
      bary = 0;
      str = 1;
      ch = 0;
      break;
    case 17: // k0bar
      bary = 0;
      str = -1;
      ch = 0;
      break;
    case 21: // lambda
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 23: // sigma+
      bary = 1;
      str = -1;
      ch = 1;
      break;
    case 25: // sigma0
      bary = 1;
      str = -1;
      ch = 0;
      break;
    case 27: // sigma-
      bary = 1;
      str = -1;
      ch = -1;
      break;
    case 29: // xi0
      bary = 1;
      str = -2;
      ch = 0;
      break;
    case 31: // xi-
      bary = 1;
      str = -2;
      ch = -1;
      break;
    default:
      G4cout << " Unknown particle type " << type << G4endl;
    };

    Qnums[0] = bary;
    Qnums[1] = str;
    Qnums[2] = ch;
    return Qnums;
}

void 
G4CascadeChannel::CheckQnums(const G4FastVector<G4ReactionProduct,256> &vec,
			     G4int &vecLen,
			     G4ReactionProduct &currentParticle,
			     G4ReactionProduct &targetParticle,
			     G4double Q, G4double B, G4double S) {
  G4ParticleDefinition* projDef = currentParticle.GetDefinition();
  G4ParticleDefinition* targDef = targetParticle.GetDefinition();
  G4double chargeSum = projDef->GetPDGCharge() + targDef->GetPDGCharge();
  G4double baryonSum = projDef->GetBaryonNumber() + targDef->GetBaryonNumber();
  G4double strangenessSum = projDef->GetQuarkContent(3) - 
                            projDef->GetAntiQuarkContent(3) + 
                            targDef->GetQuarkContent(3) -
                            targDef->GetAntiQuarkContent(3);

  G4ParticleDefinition* secDef = 0;
  for (G4int i = 0; i < vecLen; i++) {
    secDef = vec[i]->GetDefinition();
    chargeSum += secDef->GetPDGCharge();
    baryonSum += secDef->GetBaryonNumber();
    strangenessSum += secDef->GetQuarkContent(3) 
                    - secDef->GetAntiQuarkContent(3);
  }

  G4bool OK = true;
  if (chargeSum != Q) {
    G4cout << " Charge not conserved " << G4endl;
    OK = false;
  }
  if (baryonSum != B) {
    G4cout << " Baryon number not conserved " << G4endl;
    OK = false;
  }
  if (strangenessSum != S) {
    G4cout << " Strangeness not conserved " << G4endl;
    OK = false;
  } 

  if (!OK) {
    G4cout << " projectile: " << projDef->GetParticleName() 
           << "  target: " << targDef->GetParticleName() << G4endl;
    for (G4int i = 0; i < vecLen; i++) {
      secDef = vec[i]->GetDefinition();
      G4cout << secDef->GetParticleName() << " " ;
    }
    G4cout << G4endl;
  }
}
