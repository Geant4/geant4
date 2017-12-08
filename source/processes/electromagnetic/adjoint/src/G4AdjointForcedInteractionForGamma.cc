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
// $Id: G4AdjointForcedInteractionForGamma.cc 87443 2014-12-04 12:26:31Z gunter $
//
#include "G4AdjointForcedInteractionForGamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4AdjointCSManager.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4VEmAdjointModel.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointGamma.hh"


G4AdjointForcedInteractionForGamma::
	G4AdjointForcedInteractionForGamma(G4String process_name):
			G4VContinuousDiscreteProcess(process_name),theAdjointComptonModel(0),theAdjointBremModel(0)
{ theAdjointCSManager = G4AdjointCSManager::GetAdjointCSManager();
  fParticleChange=new G4ParticleChange();
  lastAdjCS=0.;
  trackid = nstep = 0;
  is_free_flight_gamma = false;
  copy_gamma_for_forced_interaction = false;
  last_free_flight_trackid=1000;

  theAdjointComptonModel =0;
  theAdjointBremModel=0;

  acc_track_length=0.;
  acc_nb_adj_interaction_length=0.;
  acc_nb_fwd_interaction_length=0.;
  total_acc_nb_adj_interaction_length=0.;
  total_acc_nb_fwd_interaction_length=0.;
  continue_gamma_as_new_free_flight =false;



}
//////////////////////////////////////////////////////////////////////////////
//
G4AdjointForcedInteractionForGamma::
	~G4AdjointForcedInteractionForGamma()
{ if (fParticleChange) delete fParticleChange;
}			
//////////////////////////////////////////////////////////////////////////////
//
void G4AdjointForcedInteractionForGamma::PreparePhysicsTable(const G4ParticleDefinition&)
{;
}
//////////////////////////////////////////////////////////////////////////////
//
void G4AdjointForcedInteractionForGamma::BuildPhysicsTable(const G4ParticleDefinition&)
{
 theAdjointCSManager->BuildCrossSectionMatrices(); //do not worry it will be done just once
 theAdjointCSManager->BuildTotalSigmaTables();
}
//Note on weight correction for forced interaction
//For the forced interaction applied here we do use a truncated exponential law for the probability of survival
//over a fixed total length. This is done by using a linear transformation of the non biased probability survival
//In mathematic this writes
//P'(x)=C1P(x)+C2
//With P(x)=exp(-sum(sigma_ixi))  x and L can cross different volumes with different cross section sigma.
//For forced interaction  we get the following limit conditions
//P'(L)=0 P'(0)=1       (L can be used over different volumes)
//From simple solving of linear equation we get
//C1=1/(1-P(L)) et C2=-P(L)/(1-P(L))
//P'(x)=(P(x)-P(L))/(1-P(L))
//For the probability  over a step x1 to x2
//P'(x1->x2)=P'(x2)/P'(x1)
//The effective cross section is defined -d(P'(x))/dx/P'(x)
//We get therefore
//sigma_eff=C1sigmaP(x)/(C1P(x)+C2)=sigmaP(x)/(P(x)+C2/C1)=sigmaP(x)/(P(x)-P(L))=sigma/(1-P(L)/P(x))
//////////////////////////////////////////////////////////////////////////////
//
G4VParticleChange* G4AdjointForcedInteractionForGamma::PostStepDoIt(const G4Track& track, const G4Step&  )
{  fParticleChange->Initialize(track);
  //For the free flight gamma no  interaction occur but a gamma with same property is
  //produces for further forced interaction
  //It is done at the very beginning of the track such that the weight can be the same
  if  (copy_gamma_for_forced_interaction) {
	  G4ThreeVector theGammaMomentum = track.GetMomentum();
	  fParticleChange->AddSecondary(new G4DynamicParticle(G4AdjointGamma::AdjointGamma(),theGammaMomentum));
	  fParticleChange->SetParentWeightByProcess(false);
	  fParticleChange->SetSecondaryWeightByProcess(false);
  }
  else { //Occurrence of forced interaction

	  //Selection of the model to be called
	  G4VEmAdjointModel* theSelectedModel =0;
	  G4bool is_scat_proj_to_proj_case=false;
	  if (!theAdjointComptonModel && !theAdjointBremModel) return fParticleChange;
	  if (!theAdjointComptonModel) {
		  theSelectedModel = theAdjointBremModel;
		  is_scat_proj_to_proj_case=false;
		  //This is needed because the results of it will be used in the post step do it weight correction inside the model
		  theAdjointBremModel->AdjointCrossSection(
		                    track.GetMaterialCutsCouple(),track.GetKineticEnergy(), false);

	  }
	  else if (!theAdjointBremModel) {
		  theSelectedModel = theAdjointComptonModel;
		  is_scat_proj_to_proj_case=true;
	  }
	  else { //Choose the model according to cross sections

        G4double bremAdjCS = theAdjointBremModel->AdjointCrossSection(
                  track.GetMaterialCutsCouple(),track.GetKineticEnergy(), false);
        if (G4UniformRand()*lastAdjCS<bremAdjCS) {
        	theSelectedModel = theAdjointBremModel;
        	is_scat_proj_to_proj_case=false;
        }
        else {
        	theSelectedModel = theAdjointComptonModel;
        	is_scat_proj_to_proj_case=true;
        }
	  }

	  //Compute the weight correction factor
	  G4double one_over_effectiveAdjointCS= (1.-std::exp(acc_nb_adj_interaction_length-total_acc_nb_adj_interaction_length))/lastAdjCS;
      G4double  weight_correction_factor = lastAdjCS*one_over_effectiveAdjointCS;
      //G4cout<<"Weight correction factor start "<<weight_correction_factor<<std::endl;
	  //Call the  selected model without correction of the weight in the model
	  theSelectedModel->SetCorrectWeightForPostStepInModel(false);
	  theSelectedModel->SetAdditionalWeightCorrectionFactorForPostStepOutsideModel(weight_correction_factor);
	  theSelectedModel->SampleSecondaries(track,is_scat_proj_to_proj_case,fParticleChange);
	  theSelectedModel->SetCorrectWeightForPostStepInModel(true);

	  continue_gamma_as_new_free_flight =true;
  }
  return fParticleChange;
}
//////////////////////////////////////////////////////////////////////////////
//
G4VParticleChange* G4AdjointForcedInteractionForGamma::AlongStepDoIt(const G4Track& track, const G4Step&  )
{ fParticleChange->Initialize(track);
  //Compute nb of interactions length over step length
  G4ThreeVector position = track.GetPosition();
  G4double stepLength = track.GetStep()->GetStepLength();
  G4double ekin = track.GetKineticEnergy();
  G4double nb_fwd_interaction_length_over_step=0.;
  G4double nb_adj_interaction_length_over_step=0.;
  lastAdjCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalAdjointCS(track.GetDefinition(), ekin, track.GetMaterialCutsCouple());
  lastFwdCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalForwardCS(G4AdjointGamma::AdjointGamma(),
          ekin,track.GetMaterialCutsCouple());
  nb_fwd_interaction_length_over_step = stepLength*lastFwdCS;
  nb_adj_interaction_length_over_step = stepLength*lastAdjCS;
  G4double fwd_survival_probability=std::exp(-nb_fwd_interaction_length_over_step);
  G4double mc_induced_survival_probability=1.;

  if (is_free_flight_gamma) { //for free_flight survival probability stays 1
      //Accumulate the number of interaction lengths during free flight of gamma
	 total_acc_nb_fwd_interaction_length+=nb_fwd_interaction_length_over_step;
     total_acc_nb_adj_interaction_length+=nb_adj_interaction_length_over_step;
     acc_track_length+=stepLength;
  }
  else {
	  G4double previous_acc_nb_adj_interaction_length =acc_nb_adj_interaction_length;
	  acc_nb_fwd_interaction_length+=nb_fwd_interaction_length_over_step;
	  acc_nb_adj_interaction_length+=nb_adj_interaction_length_over_step;
	  theNumberOfInteractionLengthLeft-=nb_adj_interaction_length_over_step;

	  //Following condition to remove very rare FPE issue
	  //if (total_acc_nb_adj_interaction_length <=  1.e-50 &&  theNumberOfInteractionLengthLeft<=1.e-50) { //condition added to avoid FPE issue
          // VI 06.11.2017 - new condition
	  if (std::abs(total_acc_nb_adj_interaction_length - previous_acc_nb_adj_interaction_length) <= 1.e-15) {
	    mc_induced_survival_probability =  1.e50;
	    /*
	    G4cout << "FPE protection: " << total_acc_nb_adj_interaction_length << " " 
		   << previous_acc_nb_adj_interaction_length << " "
		   << acc_nb_fwd_interaction_length << " " 
		   << acc_nb_adj_interaction_length << " " 
		   << theNumberOfInteractionLengthLeft 
		   << G4endl;
	    */
	  }
      else {
       mc_induced_survival_probability= std::exp(-acc_nb_adj_interaction_length)-std::exp(-total_acc_nb_adj_interaction_length);
       mc_induced_survival_probability=mc_induced_survival_probability/(std::exp(-previous_acc_nb_adj_interaction_length)-std::exp(-total_acc_nb_adj_interaction_length));
      }
  }
 G4double weight_correction = fwd_survival_probability/mc_induced_survival_probability;

 //weight_correction = 1.;
 //Caution!!!
   // It is important  to select the weight of the post_step_point
   // as the current weight and not the weight of the track, as t
   // the  weight of the track is changed after having applied all
   // the along_step_do_it.
 G4double new_weight=weight_correction*track.GetStep()->GetPostStepPoint()->GetWeight();
/*
 G4cout<<"New weight "<<new_weight<<std::endl;
 G4cout<<"Weight correction "<<weight_correction<<std::endl;
 */

 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);

 return fParticleChange;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointForcedInteractionForGamma::PostStepGetPhysicalInteractionLength(
                  const G4Track& track,
 			     G4double   ,
 			     G4ForceCondition* condition)
{ G4int step_id = track.GetCurrentStepNumber();
 *condition = NotForced;
 copy_gamma_for_forced_interaction = false;
 G4int track_id=track.GetTrackID();
 is_free_flight_gamma =  (track_id != last_free_flight_trackid+1 || continue_gamma_as_new_free_flight);
 if (is_free_flight_gamma) {
  if (step_id == 1 || continue_gamma_as_new_free_flight) {
    *condition=Forced;
    //A gamma with same conditions will be generate at next post_step do it for the forced interaction
    copy_gamma_for_forced_interaction = true;
	last_free_flight_trackid = track_id;
	acc_track_length=0.;
	total_acc_nb_adj_interaction_length=0.;
	total_acc_nb_fwd_interaction_length=0.;
	continue_gamma_as_new_free_flight=false;
	return 1.e-90;
  }
  else {
	  //Computation of accumulated length for
	  return DBL_MAX;
  }
 }
 else { //compute the interaction length for forced interaction
   if (step_id ==1) {
	   G4double min_val= std::exp(-total_acc_nb_adj_interaction_length);
	   theNumberOfInteractionLengthLeft =  -std::log( min_val+G4UniformRand()*(1.-min_val));
	   theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft;
	   acc_nb_adj_interaction_length=0.;
	   acc_nb_fwd_interaction_length=0.;
   }
   G4VPhysicalVolume* thePostPhysVolume = track.GetStep()->GetPreStepPoint()->GetPhysicalVolume();
   G4double ekin =track.GetKineticEnergy();
   G4double postCS=0.;
   if (thePostPhysVolume){
      postCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalAdjointCS(G4AdjointGamma::AdjointGamma(),
                                                    ekin,thePostPhysVolume->GetLogicalVolume()->GetMaterialCutsCouple());
   }
   if (postCS>0.) return theNumberOfInteractionLengthLeft/postCS;
   else return DBL_MAX;
  }
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointForcedInteractionForGamma::GetContinuousStepLimit(const G4Track& ,
                                G4double  ,
                                G4double  ,
                                G4double& )
{return DBL_MAX;
}
////////////////////////////////////////////////////////////////////////////////
//Not used in this process but should be implemented as virtual method
G4double G4AdjointForcedInteractionForGamma::GetMeanFreePath(const G4Track& ,
                                     G4double ,
                                      G4ForceCondition*)
{ return 0.;
}


