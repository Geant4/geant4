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
// $Id: G4AdjointSteppingAction.cc 102436 2017-01-27 08:28:53Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////
//      Class Name:	G4AdjointSteppingAction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////

#include "G4AdjointSteppingAction.hh"
#include "G4Track.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4AffineTransform.hh"
#include "G4AdjointCrossSurfChecker.hh"
////////////////////////////////////////////////////////////////////////////////////////////////
//

G4AdjointSteppingAction::G4AdjointSteppingAction()
  : ext_sourceEMax(0.), start_event(false),
    did_adj_part_reach_ext_source(false), last_ekin(0.), last_weight(0.),
    prim_weight(1.), last_part_def(0), theUserAdjointSteppingAction(0),
    theUserFwdSteppingAction(0)
{ 
  theG4AdjointCrossSurfChecker = G4AdjointCrossSurfChecker::GetInstance();
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSteppingAction::~G4AdjointSteppingAction()
{;}

////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* aTrack =aStep->GetTrack();
  //G4cout<<"Ste weight "<<<<std:.endl;
  //forward tracking mode
   if(!is_adjoint_tracking_mode){
     if (!did_one_adj_part_reach_ext_source_during_event) {
        aTrack->SetTrackStatus(fStopAndKill);
        return;
     }
     if(theUserFwdSteppingAction)
       theUserFwdSteppingAction->UserSteppingAction(aStep);
     return;
   }
  //Apply first the user adjoint stepping action
  //---------------------------
  did_adj_part_reach_ext_source=false;
  if (theUserAdjointSteppingAction) theUserAdjointSteppingAction->UserSteppingAction(aStep);


  G4double nb_nuc=1.;
  G4ParticleDefinition* thePartDef = aTrack->GetDefinition();
 
  if (thePartDef->GetParticleType() == "adjoint_nucleus"){
  	nb_nuc=double(thePartDef->GetBaryonNumber());
  }
  //Kill conditions for adjoint particles reaching the maximum energy
  //-----------------------------------------------------------------
  if(aTrack->GetKineticEnergy() >= ext_sourceEMax*nb_nuc){
	aTrack->SetTrackStatus(fStopAndKill);
	did_adj_part_reach_ext_source=false;
	return;
  }
  //G4cout<<"Weight adjoint "<<aTrack->GetWeight()<<std::endl;
  /*G4double weight_factor = aTrack->GetWeight()/prim_weight;
  if ( (weight_factor>0 && weight_factor<=0) || weight_factor<= 1e-290 || weight_factor>1.e200)
  {
	//std::cout<<"Weight_factor problem! Value = "<<weight_factor<<std::endl;
	aTrack->SetTrackStatus(fStopAndKill);
	did_adj_part_reach_ext_source=false;
	//G4cout<<thePartDef->GetParticleName()<<std::endl;
	return;	
  }
  */
  
  
  //Kill conditions for surface crossing
  //--------------------------------------
  
  G4String surface_name;
  G4double cos_to_surface;
  G4bool GoingIn;
  G4ThreeVector crossing_pos;
  if (theG4AdjointCrossSurfChecker->CrossingOneOfTheRegisteredSurface(aStep, surface_name, crossing_pos, cos_to_surface, GoingIn) ){

	//G4cout<<"Test_step11"<<std::endl;
	if (surface_name == "ExternalSource") {
		//Registering still needed
		did_adj_part_reach_ext_source=true;
		did_one_adj_part_reach_ext_source_during_event=true;
		aTrack->SetTrackStatus(fStopAndKill);
		//now register the adjoint particles reaching the external surface
		last_momentum =aTrack->GetMomentum();
		last_ekin=aTrack->GetKineticEnergy();
		last_weight = aTrack->GetWeight();
		last_part_def = aTrack->GetDefinition();
		last_pos = crossing_pos;
    		return;
	}	
	else if (surface_name == "AdjointSource" && GoingIn) {
		did_adj_part_reach_ext_source=false;
		aTrack->SetTrackStatus(fStopAndKill);
		return;
	}  
  }
  //Check for reaching out of world
  //G4cout<<aStep->GetPostStepPoint()->GetStepStatus()<<std::endl;
  if (aStep->GetPostStepPoint()->GetStepStatus() == fWorldBoundary) {
	  did_adj_part_reach_ext_source=true;
	  did_one_adj_part_reach_ext_source_during_event=true;
	  last_momentum =aTrack->GetMomentum();
	  last_ekin=aTrack->GetKineticEnergy();
	  last_weight = aTrack->GetWeight();
	  last_part_def = aTrack->GetDefinition();
	  last_pos = crossing_pos;
  		return;
  }
}

