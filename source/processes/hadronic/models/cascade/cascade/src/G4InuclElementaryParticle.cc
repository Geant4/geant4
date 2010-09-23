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
// $Id: G4InuclElementaryParticle.cc,v 1.10 2010-09-23 05:33:56 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100428  M. Kelsey -- Use G4InuclParticleNames enums instead of numbers,
//		add Omega and antinucleons.
// 20100429  M. Kelsey -- Change "case gamma:" to "case photon:"
// 20100923  M. Kelsey -- Drop "uups" message when converting G4PartDef to code

#include "G4InuclElementaryParticle.hh"

#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Gamma.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaMinus.hh"
#include "G4XiZero.hh"
#include "G4XiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4Diproton.hh"
#include "G4UnboundPN.hh"
#include "G4Dineutron.hh"

#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;


G4ParticleDefinition* 
G4InuclElementaryParticle::makeDefinition(G4int ityp) {
  switch(ityp) {
  case proton:      return G4Proton::Definition(); break;
  case neutron:     return G4Neutron::Definition(); break;
  case pionPlus:    return G4PionPlus::Definition(); break;
  case pionMinus:   return G4PionMinus::Definition(); break;
  case pionZero:    return G4PionZero::Definition(); break;
  case photon:      return G4Gamma::Definition(); break;
  case kaonPlus:    return G4KaonPlus::Definition(); break;
  case kaonMinus:   return G4KaonMinus::Definition(); break;
  case kaonZero:    return G4KaonZero::Definition(); break;
  case kaonZeroBar: return G4AntiKaonZero::Definition(); break;
  case lambda:      return G4Lambda::Definition(); break;
  case sigmaPlus:   return G4SigmaPlus::Definition(); break;
  case sigmaZero:   return G4SigmaZero::Definition(); break;
  case sigmaMinus:  return G4SigmaMinus::Definition(); break;
  case xiZero:      return G4XiZero::Definition(); break;
  case xiMinus:     return G4XiMinus::Definition(); break;
  case omegaMinus:  return G4OmegaMinus::Definition(); break;
  case antiProton:  return G4AntiProton::Definition(); break;
  case antiNeutron: return G4AntiNeutron::Definition(); break;
  case diproton:    return G4Diproton::Definition(); break;  // Bertini class!
  case unboundPN:   return G4UnboundPN::Definition(); break; // Bertini class!
  case dineutron:   return G4Dineutron::Definition(); break; // Bertini class!
  default:
    G4cerr << " uups, unknown particle type " << ityp << G4endl;
  }
  
  return 0;
}

// This is the inverse mapping to makeDefinition above

G4int G4InuclElementaryParticle::type(const G4ParticleDefinition *pd) {
  if (pd == 0) return 0;
  if (pd == G4Proton::Definition())       return proton;
  if (pd == G4Neutron::Definition())      return neutron;
  if (pd == G4PionPlus::Definition())     return pionPlus;
  if (pd == G4PionMinus::Definition())    return pionMinus;
  if (pd == G4PionZero::Definition())     return pionZero;
  if (pd == G4Gamma::Definition())        return photon;
  if (pd == G4KaonPlus::Definition())     return kaonPlus;
  if (pd == G4KaonMinus::Definition())    return kaonMinus;
  if (pd == G4KaonZero::Definition())     return kaonZero;
  if (pd == G4AntiKaonZero::Definition()) return kaonZeroBar;
  if (pd == G4Lambda::Definition())       return lambda;
  if (pd == G4SigmaPlus::Definition())    return sigmaPlus;
  if (pd == G4SigmaZero::Definition())    return sigmaZero;
  if (pd == G4SigmaMinus::Definition())   return sigmaMinus;
  if (pd == G4XiZero::Definition())       return xiZero;
  if (pd == G4XiMinus::Definition())      return xiMinus;
  if (pd == G4OmegaMinus::Definition())   return omegaMinus;
  if (pd == G4AntiProton::Definition())   return antiProton;
  if (pd == G4AntiNeutron::Definition())  return antiNeutron;
  if (pd == G4Diproton::Definition())     return diproton;  // Bertini class!
  if (pd == G4UnboundPN::Definition())    return unboundPN; // Bertini class!
  if (pd == G4Dineutron::Definition())    return dineutron; // Bertini class!

  return 0;	// Unknown objects return zero (e.g., nuclei)
}

void G4InuclElementaryParticle::setType(G4int ityp) {
  setDefinition(makeDefinition(ityp));
}


// Assignment operator for use with std::sort()
G4InuclElementaryParticle& 
G4InuclElementaryParticle::operator=(const G4InuclElementaryParticle& right) {
  generation = right.generation;
  G4InuclParticle::operator=(right);
  return *this;
}


G4double G4InuclElementaryParticle::getStrangeness(G4int type) {
  G4ParticleDefinition* pd = makeDefinition(type);
  return pd ? (pd->GetQuarkContent(3) - pd->GetAntiQuarkContent(3)) : 0.;
}

G4double G4InuclElementaryParticle::getParticleMass(G4int type) {
  G4ParticleDefinition* pd = makeDefinition(type);
  return pd ? pd->GetPDGMass()*MeV/GeV : 0.0;	// From G4 to Bertini units
}


// Print particle parameters

void G4InuclElementaryParticle::printParticle() const {
  G4InuclParticle::printParticle();
  G4cout << " Particle: " << getDefinition()->GetParticleName() 
	 << " type " << type() << " mass " << getMass()
	 << " ekin " << getKineticEnergy() << G4endl; 
}

