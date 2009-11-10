#include "G4AdjointSteppingAction.hh"
#include "G4Track.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4AffineTransform.hh"
#include "G4AdjointCrossSurfChecker.hh"
////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSteppingAction::G4AdjointSteppingAction()
{ 
  theG4AdjointCrossSurfChecker = G4AdjointCrossSurfChecker::GetInstance();
  did_adj_part_reach_external_source =false;
  theUserAdjointSteppingAction =0;
}
////////////////////////////////////////////////////////////////////////////////////////////////
//
G4AdjointSteppingAction::~G4AdjointSteppingAction()
{;}

////////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //Apply first the user adjoint stepping action
  //---------------------------
  if (theUserAdjointSteppingAction) theUserAdjointSteppingAction->UserSteppingAction(aStep);
  
  
  
  
  G4Track* aTrack =aStep->GetTrack();
  G4double nb_nuc=1.;
  G4ParticleDefinition* thePartDef = aTrack->GetDefinition();
 
  if (thePartDef->GetParticleType() == "adjoint_nucleus"){
  	nb_nuc=double(thePartDef->GetBaryonNumber());
  }
  //Kill conditions for adjoint particles reaching the maximum energy
  //-----------------------------------------------------------------
  if(aTrack->GetKineticEnergy() >= ext_sourceEMax*nb_nuc){
	aTrack->SetTrackStatus(fStopAndKill);
	did_adj_part_reach_external_source=false;
	return;
  }

  G4double weight_factor = aTrack->GetWeight()/prim_weight;

#ifdef WIN32
  if (!!_isnan(weight_factor) || weight_factor<= 1e-290 || weight_factor>1.e200){
#else
  if (std::isnan(weight_factor) || weight_factor<= 1e-290 || weight_factor>1.e200){
#endif
	
	//std::cout<<"Weight_factor problem! Value = "<<weight_factor<<std::endl;
	aTrack->SetTrackStatus(fStopAndKill);
	did_adj_part_reach_external_source=false;
	return;
		
  }
  
  
  //Kill conditions for surface crossing
  //--------------------------------------
  
  G4String surface_name;
  G4double cos_to_surface;
  G4bool GoingIn;
  G4ThreeVector crossing_pos;
  if (theG4AdjointCrossSurfChecker->CrossingOneOfTheRegisteredSurface(aStep, surface_name, crossing_pos, cos_to_surface, GoingIn) ){
  	
	//G4cout<<"Test_step11"<<std::endl;
	//G4cout<<surface_name<<std::endl;
	if (surface_name == "ExternalSource") {
		//Registering still needed
		did_adj_part_reach_external_source=true;
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
		did_adj_part_reach_external_source=false;
		return;
	}	
	
  
  }
  
  		
  
}
