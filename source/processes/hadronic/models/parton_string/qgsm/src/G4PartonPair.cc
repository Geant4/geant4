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
#include "G4PartonPair.hh"
#include "G4HadronicException.hh"

//#define DEBUG_PartonPair 1

G4PartonPair::G4PartonPair(G4Parton* P1, G4Parton* P2, G4int Type, G4int aDirection)
{
  CollisionType = Type;
  Parton1 = P1;
  Parton2 = P2;
  Direction = aDirection;

  #ifdef DEBUG_PartonPair
  G4cout << "ctor G4PartonPair - " 
	 << (aDirection==PROJECTILE ? "Projectile": "Target") 
	 << (CollisionType==SOFT ? " Soft " : " Diffractive " ) << G4endl
	 << "  Parton 1 name, type, spin-3, colour, 4-mom " 
	 << P1->GetDefinition()->GetParticleName() << ", " 
	 << P1->GetDefinition()->GetParticleType() << ", " 
	 << P1->GetSpinZ() << ", "
	 << P1->GetColour() << ", "
	 << P1->Get4Momentum() << " " << G4endl
	 << "  Parton 2 name, type, spin-3, colour, 4-mom " 
	 << P2->GetDefinition()->GetParticleName() << ", " 
	 << P2->GetDefinition()->GetParticleType() << ", " 
	 << P2->GetSpinZ() << ", "
	 << P2->GetColour() << ", "
	 << P2->Get4Momentum() << G4endl
	 << " string mass, 4mom " 
	 << (P1->Get4Momentum()+P2->Get4Momentum()).m() << " "
	 << (P1->Get4Momentum()+P2->Get4Momentum()) << G4endl;
  #endif      
}

G4PartonPair::~G4PartonPair()
{
}

