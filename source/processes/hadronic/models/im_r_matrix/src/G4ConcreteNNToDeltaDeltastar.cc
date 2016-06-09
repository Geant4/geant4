//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ConcreteNNToDeltaDeltastar.cc,v 1.4 2004/12/07 13:48:47 gunter Exp $ //

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ConcreteNNToDeltaDeltastar.hh"
#include "G4DeltaDeltastarBuilder.hh"
#include <typeinfo>

G4XDeltaDeltastarTable G4ConcreteNNToDeltaDeltastar::theSigmaTable;

G4ConcreteNNToDeltaDeltastar::G4ConcreteNNToDeltaDeltastar(const G4ParticleDefinition* aPrimary,
					   const G4ParticleDefinition* bPrimary,
					   const G4ParticleDefinition* aSecondary,
					   const G4ParticleDefinition* bSecondary)
  : G4ConcreteNNTwoBodyResonance(aPrimary, bPrimary, aSecondary, bSecondary,
                                 G4DeltaDeltastarBuilder(bSecondary->GetParticleName(), theSigmaTable))
{
  G4double chargeBalance = aPrimary->GetPDGCharge()+bPrimary->GetPDGCharge();
  chargeBalance -= aSecondary->GetPDGCharge();
  chargeBalance -= bSecondary->GetPDGCharge();
  if(std::abs(chargeBalance) >.1)
  {
    G4cout << "Charge conservation problem in G4ConcreteNNToDeltaDeltastar"<<G4endl;
    G4cout << "Initial charges in "<<typeid(*this).name()<<G4endl;
    G4cout << aPrimary->GetPDGCharge()<<" "<<aPrimary->GetParticleName()
           << bPrimary->GetPDGCharge()<<" "<<bPrimary->GetParticleName()
	   << aSecondary->GetPDGCharge()<<" "<<aSecondary->GetParticleName()
	   << bSecondary->GetPDGCharge()<<" "<<bSecondary->GetParticleName()<<G4endl;
  }
}

G4ConcreteNNToDeltaDeltastar::~G4ConcreteNNToDeltaDeltastar()
{ 
}
