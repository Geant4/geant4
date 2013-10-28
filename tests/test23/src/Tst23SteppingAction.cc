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
#include "Tst23SteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"

#include "Tst23Histo.hh"
#include "Tst23ParticleChange.hh"

#include "G4RunManager.hh"
#include "TstPrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"


Tst23SteppingAction::Tst23SteppingAction( TstHisto* hptr ) 
{

   fHistoPtr  = hptr;
   
   fFirstInter = new Tst23ParticleChange( true );
   fOtherInter = new Tst23ParticleChange();

}

Tst23SteppingAction::~Tst23SteppingAction() 
{

   if ( fFirstInter ) delete fFirstInter;
   if ( fOtherInter ) delete fOtherInter;

}


void Tst23SteppingAction::UserSteppingAction( const G4Step * theStep ) {

   G4VPhysicalVolume* vol = theStep->GetTrack()->GetVolume();
   
   // G4cout << "Volume name = " << vol->GetName() << G4endl;
   
   if ( vol != fTargetPtr ) return; // we're outside primary target
      
   int nsc = fFirstInter->GetNumberOfSecondaries();
   if ( nsc > 0  )
   {
         for ( int i=0; i<nsc; i++) 
         {   
            delete fFirstInter->GetSecondary(i);
         } 
         fFirstInter->Clear();    
   }
   
   nsc = fOtherInter->GetNumberOfSecondaries();
   if ( nsc > 0 )
   {
         for ( int i=0; i<nsc; i++) 
         {   
            delete fOtherInter->GetSecondary(i);
         } 
         fOtherInter->Clear();    
   }
   
   const std::vector<const G4Track*>* secs = theStep->GetSecondaryInCurrentStep();
   int nsec = secs->size();
   
/*
   // the beam particle can go through without interactions or at least without changing its identity;
   // in such case most likely procName=Transportation
   // in case on "read" interaction with a nucleus, it's likely to be procName=ProtonInelastic,
   // but even then there isn't much way to say what exactly caused it because there's a "model" that
   // attaches to G4ProtonInelaticProcess object (via RegisterMe in the builder)
   
   G4String procName = "";
   if ( theStep->GetTrack()->GetParentID() == 0 && theStep->GetTrack()->GetTrackStatus() != fAlive )
   {
      procName = theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
      std::cout << " process name = " << procName << std::endl;
   }
*/
   
   for ( int i=0; i<nsec; i++ )
   {

      G4Track* tr = new G4Track( *((*secs)[i]) );
      
      if ( theStep->GetTrack()->GetTrackStatus() != fAlive ) // track losses identity
      {
         if ( theStep->GetTrack()->GetParentID() == 0 ) // primary track
         {
	    fFirstInter->AddSecondary( tr );
         }      
         else // secondary track, and it's also looses identity (re-interaction)
         {
	    fOtherInter->AddSecondary( tr );
         }
      }

   } //end loop over secondaries
   
   const TstPrimaryGeneratorAction* beam = 
         dynamic_cast<const TstPrimaryGeneratorAction*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
   
   fHistoPtr->FillEvt( fFirstInter, beam->GetLabV(), beam->GetLabP() );
   fHistoPtr->FillEvt( fOtherInter, beam->GetLabV(), beam->GetLabP() );
      
   return;
   
}
