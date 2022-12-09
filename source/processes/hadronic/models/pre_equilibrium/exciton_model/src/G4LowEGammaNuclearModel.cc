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
// Physics model class G4LowEGammaNuclearModel 
// Created:  15 May 2019
// Author  V.Ivanchenko
//  
//

#include "G4LowEGammaNuclearModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Fragment.hh"
#include "G4NucleiProperties.hh"
#include "G4DynamicParticle.hh"
#include "G4HadSecondary.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4HadronicParameters.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4PreCompoundModel.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicsModelCatalog.hh"

G4LowEGammaNuclearModel::G4LowEGammaNuclearModel() 
  : G4HadronicInteraction("GammaNPreco"),lab4mom(0.,0.,0.,0.), secID(-1)
{
  secID = G4PhysicsModelCatalog::GetModelID("model_" + GetModelName());
  SetMinEnergy( 0.0*CLHEP::GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  // reuse existing pre-compound model
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  fPreco = static_cast<G4PreCompoundModel*>(p);
  if(!fPreco) { fPreco = new G4PreCompoundModel(); }
}

G4LowEGammaNuclearModel::~G4LowEGammaNuclearModel()
{}

void G4LowEGammaNuclearModel::InitialiseModel()
{}

G4HadFinalState* G4LowEGammaNuclearModel::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  theParticleChange.Clear();

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();

  // Create initial state
  lab4mom.set(0.,0.,0.,G4NucleiProperties::GetNuclearMass(A, Z));
  lab4mom += aTrack.Get4Momentum();

  G4Fragment frag(A, Z, lab4mom);

  frag.SetCreatorModelID(secID);

  if (verboseLevel > 1) {
    G4cout << "G4LowEGammaNuclearModel::ApplyYourself initial G4Fragmet:" 
	   << G4endl;
    G4cout << frag << G4endl;
  }
  G4ReactionProductVector* res = fPreco->DeExcite(frag);

  // secondaries produced
  if(res) {

    theParticleChange.SetStatusChange(stopAndKill);
    std::size_t nsec = res->size();
    if (verboseLevel > 1) {
      G4cout << "G4LowEGammaNuclearModel: " << nsec << " secondaries" << G4endl;
    }
    for(std::size_t i=0; i<nsec; ++i) {
      if((*res)[i]) {
	G4double ekin = (*res)[i]->GetKineticEnergy();
	G4ThreeVector dir(0.,0.,1.);
	if(ekin > 0.0) { dir = (*res)[i]->GetMomentum().unit(); }
	G4HadSecondary* news = new G4HadSecondary(
          new G4DynamicParticle((*res)[i]->GetDefinition(), dir, ekin));
	news->SetTime((*res)[i]->GetTOF());
	news->SetCreatorModelID(secID);
	theParticleChange.AddSecondary(*news);
	if (verboseLevel > 1) {
	  G4cout << i << ". " << (*res)[i]->GetDefinition()->GetParticleName()
		 << " Ekin(MeV)= " << ekin/MeV
		 << " dir: " << dir << G4endl;
	}
	delete (*res)[i];
        delete news;
      }
    } 
    delete res;
  } 
  return &theParticleChange;
}

