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
// $Id: G4InuclElementaryParticle.cc 69638 2013-05-09 04:26:00Z mkelsey $
// Geant4 tag: $Name:  $
//
// 20100428  M. Kelsey -- Use G4InuclParticleNames enums instead of numbers,
//		add Omega and antinucleons.
// 20100429  M. Kelsey -- Change "case gamma:" to "case photon:"
// 20100923  M. Kelsey -- Drop "uups" message when converting G4PartDef to code
// 20101029  M. Kelsey -- Add instantiation of new particles, antiparticles
// 20110214  M. Kelsey -- Drop unused "generation"
// 20110307  M. Kelsey -- Add random K0 mixing if K0S/K0L passed to type()
// 20110321  M. Kelsey -- Fix getStrangeness to return int
// 20110801  M. Kelsey -- Add fill() functions to replicate ctors, allowing
//		reuse of objects as buffers; c.f. G4InuclNuclei.
// 20110922  M. Kelsey -- Add stream argument to printParticle() => print()
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20130508  D. Wright -- Add lepton construction, use wrapper header
// 20140310  M. Kelsey -- Fix constness in G4PD* passing

#include "G4InuclElementaryParticle.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Diproton.hh"
#include "G4UnboundPN.hh"
#include "G4Dineutron.hh"
#include "Randomize.hh"

#include "G4InuclParticleNames.hh"
using namespace G4InuclParticleNames;


const G4ParticleDefinition* 
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
    // NOTE:  The four light nuclei "particles" are actually G4Ions
  case deuteron:    return G4Deuteron::Definition(); break;
  case triton:      return G4Triton::Definition(); break;
  case He3:	    return G4He3::Definition(); break;
  case alpha:	    return G4Alpha::Definition(); break;
  case antiProton:  return G4AntiProton::Definition(); break;
  case antiNeutron: return G4AntiNeutron::Definition(); break;
    // NOTE:  The the four light antinuclei "particles" are actually G4Ions
  case antiDeuteron: return G4AntiDeuteron::Definition(); break;
  case antiTriton:  return G4AntiTriton::Definition(); break;
  case antiHe3:     return G4AntiHe3::Definition(); break;
  case antiAlpha:   return G4AntiAlpha::Definition(); break;
    // NOTE:  The three unbound dibaryons are local Bertini classes
  case diproton:    return G4Diproton::Definition(); break;
  case unboundPN:   return G4UnboundPN::Definition(); break;
  case dineutron:   return G4Dineutron::Definition(); break;
    // Leptons are included for muon capture and future tau/neutrino physics
  case electron:    return G4Electron::Definition(); break;
  case positron:    return G4Positron::Definition(); break;
  case electronNu:  return G4NeutrinoE::Definition(); break;
  case antiElectronNu: return G4AntiNeutrinoE::Definition(); break;
  case muonMinus:   return G4MuonMinus::Definition(); break;
  case muonPlus:    return G4MuonPlus::Definition(); break;
  case muonNu:      return G4NeutrinoMu::Definition(); break;
  case antiMuonNu:  return G4AntiNeutrinoMu::Definition(); break;
  case tauMinus:    return G4TauMinus::Definition(); break;
  case tauPlus:     return G4TauPlus::Definition(); break;
  case tauNu:       return G4NeutrinoTau::Definition(); break;
  case antiTauNu:   return G4AntiNeutrinoTau::Definition(); break;
  default:
    G4cerr << "G4InuclElementaryParticle::makeDefinition: unknown particle type "
           << ityp << G4endl;
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
  // NOTE:  The four light nuclei "particles" are actually G4Ions
  if (pd == G4Deuteron::Definition())     return deuteron;
  if (pd == G4Triton::Definition())       return triton;
  if (pd == G4He3::Definition())          return He3;
  if (pd == G4Alpha::Definition())        return alpha;
  if (pd == G4AntiProton::Definition())   return antiProton;
  if (pd == G4AntiNeutron::Definition())  return antiNeutron;
  // NOTE:  The the four light antinuclei "particles" are actually G4Ions
  if (pd == G4AntiDeuteron::Definition()) return antiDeuteron;
  if (pd == G4AntiTriton::Definition())   return antiTriton;
  if (pd == G4AntiHe3::Definition())      return antiHe3;
  if (pd == G4AntiAlpha::Definition())    return antiAlpha;
  // NOTE:  The three unbound dibaryons are local Bertini classes
  if (pd == G4Diproton::Definition())     return diproton;
  if (pd == G4UnboundPN::Definition())    return unboundPN;
  if (pd == G4Dineutron::Definition())    return dineutron;

  if (pd == G4Electron::Definition())     return electron;
  if (pd == G4Positron::Definition())     return positron;
  if (pd == G4NeutrinoE::Definition())    return electronNu;
  if (pd == G4AntiNeutrinoE::Definition()) return antiElectronNu;
  if (pd == G4MuonMinus::Definition())    return muonMinus;
  if (pd == G4MuonPlus::Definition())     return muonPlus;
  if (pd == G4NeutrinoMu::Definition())   return muonNu;
  if (pd == G4AntiNeutrinoMu::Definition()) return antiMuonNu;
  if (pd == G4TauMinus::Definition())     return tauMinus;
  if (pd == G4TauPlus::Definition())      return tauPlus;
  if (pd == G4NeutrinoTau::Definition())  return tauNu;
  if (pd == G4AntiNeutrinoTau::Definition()) return antiTauNu;

  // Weak neutral kaons must be mixed back to strong (strangeness states)
  if (pd==G4KaonZeroShort::Definition() || pd==G4KaonZeroLong::Definition()) {
    return ((G4UniformRand() > 0.5) ? kaonZero : kaonZeroBar);
  }

  return 0;	// Unknown objects return zero (e.g., nuclei)
}

void G4InuclElementaryParticle::setType(G4int ityp) {
  setDefinition(makeDefinition(ityp));
}


// Overwrite data structure (avoids creating/copying temporaries)

void G4InuclElementaryParticle::fill(const G4LorentzVector& mom, G4int ityp,
				     G4InuclParticle::Model model) {
  setType(ityp);
  setMomentum(mom);
  setModel(model);
}

void G4InuclElementaryParticle::fill(G4double ekin, G4int ityp,
				     G4InuclParticle::Model model) {
  setType(ityp);
  setKineticEnergy(ekin);
  setModel(model);
}

void G4InuclElementaryParticle::fill(const G4LorentzVector& mom,
				     const G4ParticleDefinition* pd,
				     G4InuclParticle::Model model) {
  setDefinition(pd);
  setMomentum(mom);
  setModel(model);
}


// Assignment operator for use with std::sort()
G4InuclElementaryParticle& 
G4InuclElementaryParticle::operator=(const G4InuclElementaryParticle& right) {
  G4InuclParticle::operator=(right);
  return *this;
}


G4int G4InuclElementaryParticle::getStrangeness(G4int ityp) {
  const G4ParticleDefinition* pd = makeDefinition(ityp);
  return pd ? (pd->GetQuarkContent(3) - pd->GetAntiQuarkContent(3)) : 0;
}

G4double G4InuclElementaryParticle::getParticleMass(G4int ityp) {
  const G4ParticleDefinition* pd = makeDefinition(ityp);
  return pd ? pd->GetPDGMass()*MeV/GeV : 0.0;	// From G4 to Bertini units
}


// Print particle parameters

void G4InuclElementaryParticle::print(std::ostream& os) const {
  G4InuclParticle::print(os);
  os << G4endl << " Particle: " << getDefinition()->GetParticleName() 
     << " type " << type() << " mass " << getMass()
     << " ekin " << getKineticEnergy(); 
}

