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
//
 // Hadronic Process: Energy Range Manager
 // original by H.P. Wellisch
 // modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Last modified: 24-Mar-1997
 // fix in the counter-hndling: H.P. Wellisch 04-Apr-97
 // throw an exception if no model found:  J.L. Chuma  04-Apr-97
 
#include "G4EnergyRangeManager.hh"
#include "Randomize.hh"
#include "G4HadronicException.hh"
#include "G4SystemOfUnits.hh"

G4EnergyRangeManager::G4EnergyRangeManager()
 : theHadronicInteractionCounter(0)
{}

G4EnergyRangeManager::~G4EnergyRangeManager()
{}

void G4EnergyRangeManager::RegisterMe(G4HadronicInteraction* a)
{
  if(!a) { return; }
  if(0 < theHadronicInteractionCounter) {
    for(G4int i=0; i<theHadronicInteractionCounter; ++i) {
      if(a == theHadronicInteraction[i]) { return; }
    }
  }
  theHadronicInteraction.push_back(a);
  ++theHadronicInteractionCounter;
}

G4HadronicInteraction*
G4EnergyRangeManager::GetHadronicInteraction(const G4HadProjectile & aHadProjectile, 
                                             G4Nucleus & aTargetNucleus,
                                             const G4Material* aMaterial,
                                             const G4Element* anElement) const
{
  // VI shortcut: if only one interaction is registered skip all checks
  if(1 == theHadronicInteractionCounter) { return theHadronicInteraction[0]; }
  else if(0 == theHadronicInteractionCounter) {
    G4cout << "G4EnergyRangeManager::GetHadronicInteraction: "
	   << "no models defined for a process" << G4endl;
    return nullptr;
  }

  G4double kineticEnergy = aHadProjectile.GetKineticEnergy();
  // For ions, get kinetic energy per nucleon
  if ( std::abs( aHadProjectile.GetDefinition()->GetBaryonNumber() ) > 1 ) {
    kineticEnergy /= static_cast< G4double >( std::abs( aHadProjectile.GetDefinition()->GetBaryonNumber() ) );
  }

  G4int cou = 0, memory = 0, memor2 = 0;
  G4double emi1 = 0.0, ema1 = 0.0, emi2 = 0.0, ema2 = 0.0;

  for (G4int i = 0; i<theHadronicInteractionCounter; ++i) {
    if ( theHadronicInteraction[i]->IsApplicable( aHadProjectile, aTargetNucleus ) ) {
      G4double low  = theHadronicInteraction[i]->GetMinEnergy( aMaterial, anElement );
      G4double high = theHadronicInteraction[i]->GetMaxEnergy( aMaterial, anElement );
      if (low <= kineticEnergy && high >= kineticEnergy) {
        ++cou;
        emi2 = emi1;
        ema2 = ema1;
        emi1 = low;
        ema1 = high;
        memor2 = memory;
        memory = i;
      }
    }
  }

  G4HadronicInteraction* hi = nullptr;
  switch (cou) {
    case 0:
      G4cout << "No model found out of " << theHadronicInteractionCounter << G4endl;
      for( G4int j=0; j<theHadronicInteractionCounter; ++j) {
	G4HadronicInteraction* hint=theHadronicInteraction[j];
	G4cout << "   "<< j << ".  Elow= " << hint->GetMinEnergy(aMaterial,anElement)
	       <<", Ehigh= " << hint->GetMaxEnergy(aMaterial,anElement)
	       <<"   " << hint->GetModelName() << G4endl;
      }
      break;

    case 1:
      hi = theHadronicInteraction[memory];
      break;

    case 2:
      if( (emi2<=emi1 && ema2>=ema1) || (emi2>=emi1 && ema2<=ema1) ) {
	G4cout << "Energy ranges of two models fully overlapping " << G4endl;
	for( G4int j=0; j<theHadronicInteractionCounter; ++j) {
	  G4HadronicInteraction* hint=theHadronicInteraction[j];
	  G4cout << "   "<< j << ".  Elow= " << hint->GetMinEnergy(aMaterial,anElement)
		 <<", Ehigh= " << hint->GetMaxEnergy(aMaterial,anElement)
		 <<"   " << hint->GetModelName() << G4endl;
	}
      } else {
	G4double rand = G4UniformRand();
	G4int mem;
	if( emi1 < emi2 ) {
	  if( (ema1-kineticEnergy) < rand*(ema1-emi2) ) {
	    mem = memor2;
	  } else {
	    mem = memory;
	  }
	} else {
	  if( (ema2-kineticEnergy) < rand*(ema2-emi1) ) {
	    mem = memory;
	  } else {
	    mem = memor2;
	  }
	}
	hi = theHadronicInteraction[mem];
      }
      break;

    default:
      G4cout << "More than two competing models for this energy" << G4endl;
      for( G4int j=0; j<theHadronicInteractionCounter; ++j) {
	G4HadronicInteraction* hint=theHadronicInteraction[j];
	G4cout << "   "<< j << ".  Elow= " << hint->GetMinEnergy(aMaterial,anElement)
	       <<", Ehigh= " << hint->GetMaxEnergy(aMaterial,anElement)
	       <<"   " << hint->GetModelName() << G4endl;
      }
      break;
  }
  return hi;
} 

std::vector<G4HadronicInteraction*>& 
G4EnergyRangeManager::GetHadronicInteractionList()
{
  return theHadronicInteraction;
}

void G4EnergyRangeManager::Dump( G4int verbose )
{
  G4cout << "G4EnergyRangeManager " << this << G4endl;
  for (G4int i = 0 ; i < theHadronicInteractionCounter; i++) {
    G4cout << "   HadronicModel " << i <<":"
           << theHadronicInteraction[i]->GetModelName() << G4endl;
    if (verbose > 0) {
      G4cout << "      Minimum Energy " 
	     << theHadronicInteraction[i]->GetMinEnergy()/GeV << " [GeV], "
             << "Maximum Energy " 
	     << theHadronicInteraction[i]->GetMaxEnergy()/GeV << " [GeV]"
             << G4endl;
    }
  }
}

void
G4EnergyRangeManager::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  for (auto & hadi : theHadronicInteraction) {
    hadi->BuildPhysicsTable( aParticleType );
  }
}

 
