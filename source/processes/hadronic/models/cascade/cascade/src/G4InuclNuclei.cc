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
// $Id: G4InuclNuclei.cc,v 1.1 2010-01-12 06:27:15 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $

#include "G4InuclNuclei.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include <assert.h>


G4Allocator<G4InuclNuclei> anInuclNucleiAllocator;

// Convert nuclear configuration to standard GEANT4 pointer

// WARNING:  Opposite conventions!  G4InuclNuclei uses (A,Z) everywhere, while
//	  G4ParticleTable::GetIon() uses (Z,A)!

G4ParticleDefinition* 
G4InuclNuclei::makeDefinition(G4double a, G4double z, G4double exc) {
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition *pd = pTable->GetIon(G4int(z), G4int(a), exc);
  assert(0 != pd);
  return pd;
}

G4double G4InuclNuclei::getNucleiMass(G4double a, G4double z) {
  G4ParticleDefinition* pd = makeDefinition(a,z);
  return pd ? pd->GetPDGMass()*MeV/GeV : 0.;	// From G4 to Bertini units
}

// Assignment operator for use with std::sort()
G4InuclNuclei& G4InuclNuclei::operator=(const G4InuclNuclei& right) {
  exitationEnergy = right.exitationEnergy;
  theExitonConfiguration = right.theExitonConfiguration;
  G4InuclParticle::operator=(right);
  return *this;
}
