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
#include "G4LENDorBERTModel.hh"
#include "G4LENDCombinedModel.hh"
#include "G4CascadeInterface.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsModelCatalog.hh"

G4LENDorBERTModel::G4LENDorBERTModel( G4ParticleDefinition* pd )
  :G4LENDModel( "LENDorBERTModel" ), secID( -1 ) {
   proj = pd;
   lend = new G4LENDCombinedModel( proj ); 
   bert = new G4CascadeInterface;
   secID = G4PhysicsModelCatalog::GetModelID( "model_" + GetModelName() );
}

 #include "G4UImanager.hh"
 #include "G4SystemOfUnits.hh"

 Em10TrackingAction::Em10TrackingAction() 
 : G4UserTrackingAction()
 { }

   G4int iZ = aTarg.GetZ_asInt();
   G4int iA = aTarg.GetA_asInt();
   G4int iM = 0;
   if ( aTarg.GetIsotope() != nullptr ) iM = aTarg.GetIsotope()->Getm();

   G4DynamicParticle* dp = new G4DynamicParticle( aTrack.GetDefinition() , G4ThreeVector(0.,0.,1.) , aTrack.GetKineticEnergy() );
   G4bool lendIsOK = lend->HasData( dp , iZ , iA , iM , aTarg.GetIsotope() , nullptr , aTrack.GetMaterial() );
   delete dp;

   G4HadronicInteraction* model = nullptr;
   if ( lendIsOK ) { 
      //G4cout << "LEND is selected" << G4endl;
      model = lend;
   } else { 
      //G4cout << "BERT is selected" << G4endl;
      model = bert;
   }

   G4HadFinalState* result = model->ApplyYourself(aTrack,aTarg);
   
   // Assign the creator model ID to the secondaries
   if ( result != nullptr ) {
     for ( size_t i = 0; i < result->GetNumberOfSecondaries(); ++i ) {
       result->GetSecondary( i )->SetCreatorModelID( secID );
     }
   }

   return result;
}
