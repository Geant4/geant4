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
#include "G4AtomicDeexcitation.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"
#include "G4DynamicParticle.hh"


int main() { 
  G4int Z;
  G4int vacancyId;
  G4int numberOfRun;

  G4cout << "Enter Z " << G4endl;
  G4cin >> Z;
  G4cout << "Enter the id of the vacancy" << G4endl;
  G4cin >> vacancyId;
  G4cout<<"Enter the number of runs "<<G4endl;
  G4cin>> numberOfRun;

  G4AtomicDeexcitation* deexcitation = new G4AtomicDeexcitation;

  deexcitation->ActivateAugerElectronProduction(1);

  for(G4int i = 0; i<numberOfRun;i++)
    {G4cout<<"begin of run "<<i<<G4endl;
      G4std::vector<G4DynamicParticle*>* vectorOfParticles;
      
      vectorOfParticles = deexcitation-> GenerateParticles(Z,vacancyId);
      
      G4cout<<  vectorOfParticles->size()<<" particles in the vector "<<G4endl;
  
  for (G4int k=0; k< vectorOfParticles->size();k++)
    {
      G4DynamicParticle* newParticle = (*vectorOfParticles)[k];
      if ( newParticle->GetDefinition()->GetParticleName() == "e-" )
	{
	  G4cout <<" An auger has been generated"<<G4endl;
	G4cout<<" vectorOfParticles ["<<k<<"]:"<<G4endl;
	G4cout<<"Non zero particle. Index: "<<k<<G4endl;

      G4DynamicParticle* newElectron = (*vectorOfParticles)[k];

      
      G4ThreeVector augerDirection =newElectron ->GetMomentum();

      G4double  augerEnergy =newElectron ->GetKineticEnergy();
      G4cout<< "The Auger electron has a kinetic energy = "<<augerEnergy
          <<" MeV " <<G4endl;
	}
      else{
	G4cout<<" vectorOfParticles ["<<k<<"]:"<<G4endl;
	G4cout<<"Non zero particle. Index: "<<k<<G4endl;
      
	G4ThreeVector photonDirection = newParticle ->GetMomentum();
      
      G4double  photonEnergy =newParticle ->GetKineticEnergy();
      G4cout<< "The photon has a kinetic energy = "<<photonEnergy
          <<" MeV " <<G4endl;
      }
    }
  delete vectorOfParticles;
    }
  
  delete deexcitation;
  G4cout<<"END OF THE MAIN PROGRAM"<<G4endl;
}
