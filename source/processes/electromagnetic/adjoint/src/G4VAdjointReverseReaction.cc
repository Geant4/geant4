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
// $Id: G4VAdjointReverseReaction.cc 87443 2014-12-04 12:26:31Z gunter $
//
#include "G4VAdjointReverseReaction.hh"
#include "G4SystemOfUnits.hh"
#include "G4AdjointCSManager.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4AdjointInterpolator.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4VEmAdjointModel.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4AdjointCSManager.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"


G4VAdjointReverseReaction::
	G4VAdjointReverseReaction(G4String process_name, G4bool whichScatCase):
			G4VDiscreteProcess(process_name)
{theAdjointCSManager = G4AdjointCSManager::GetAdjointCSManager();
 IsScatProjToProjCase=whichScatCase;
 fParticleChange=new G4ParticleChange();
 IsFwdCSUsed=false;
 IsIntegralModeUsed=false;
 lastCS=0.;
 trackid = nstep = 0;
}
//////////////////////////////////////////////////////////////////////////////
//
G4VAdjointReverseReaction::
	~G4VAdjointReverseReaction()
{ if (fParticleChange) delete fParticleChange;
}			
//////////////////////////////////////////////////////////////////////////////
//
void G4VAdjointReverseReaction::PreparePhysicsTable(const G4ParticleDefinition&)
{;
}
//////////////////////////////////////////////////////////////////////////////
//
void G4VAdjointReverseReaction::BuildPhysicsTable(const G4ParticleDefinition&)
{

 theAdjointCSManager->BuildCrossSectionMatrices(); //do not worry it will be done just once
 theAdjointCSManager->BuildTotalSigmaTables();

}
//////////////////////////////////////////////////////////////////////////////
//
G4VParticleChange* G4VAdjointReverseReaction::PostStepDoIt(const G4Track& track, const G4Step&  )
{ 
  
  fParticleChange->Initialize(track);
 
 /* if (IsFwdCSUsed && IsIntegralModeUsed){ //INtegral mode still unstable
  	 G4double Tkin = step.GetPostStepPoint()->GetKineticEnergy();
  	 G4double fwdCS = theAdjointCSManager->GetTotalForwardCS(track.GetDefinition(), Tkin, track.GetMaterialCutsCouple());
	 //G4cout<<"lastCS "<<lastCS<<G4endl;
	 if (fwdCS<lastCS*G4UniformRand()) { // the reaction does not take place, same integral method as the one used for forward ionisation in  G4 
	 	ClearNumberOfInteractionLengthLeft();
  		return fParticleChange;
	 } 
	 
  }
 */
 
  theAdjointEMModel->SampleSecondaries(track,
                                       IsScatProjToProjCase,
					fParticleChange);
  
  ClearNumberOfInteractionLengthLeft();
  return fParticleChange;
  			
   
  
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4VAdjointReverseReaction::GetMeanFreePath(const G4Track& track,
                                         					     G4double ,
                                         					     G4ForceCondition* condition)
{ *condition = NotForced;
  G4double preStepKinEnergy = track.GetKineticEnergy();

  if(track.GetTrackID() != trackid) {
    trackid = track.GetTrackID();
    nstep = 0;
  }
  ++nstep;  


  
  /*G4double Sigma =
  		theAdjointEMModel->AdjointCrossSection(track.GetMaterialCutsCouple(),preStepKinEnergy,IsScatProjToProjCase);*/
  		
  G4double Sigma =
  		theAdjointEMModel->GetAdjointCrossSection(track.GetMaterialCutsCouple(),preStepKinEnergy,IsScatProjToProjCase);	

  //G4double sig = Sigma;

  G4double fwd_TotCS;
  G4double corr =  theAdjointCSManager->GetCrossSectionCorrection(track.GetDefinition(),preStepKinEnergy,track.GetMaterialCutsCouple(),IsFwdCSUsed, fwd_TotCS);

  if(std::fabs(corr) > 100.) { Sigma = 0.0; }
  else { Sigma *= corr; }

  //G4cout<<fwd_TotCS<<G4endl;
  /*if (IsFwdCSUsed && IsIntegralModeUsed){ //take the maximum cross section only for charged particle		
  	G4double e_sigma_max, sigma_max;
	theAdjointCSManager->GetMaxFwdTotalCS(track.GetDefinition(),
	 		   	     track.GetMaterialCutsCouple(), e_sigma_max, sigma_max);
	if (e_sigma_max > preStepKinEnergy){
		Sigma*=sigma_max/fwd_TotCS;
	}		     
  }
  */		

  G4double mean_free_path = 1.e60 *mm; 
  if (Sigma>0) mean_free_path = 1./Sigma;
  lastCS=Sigma;
  /*
  if(nstep > 100) {
  
    G4cout << "#* " << track.GetDefinition()->GetParticleName()
	   << " " << GetProcessName() 
	   << " Nstep " << nstep 
	   << " E(MeV)= " <<  preStepKinEnergy << "  Sig0= " << sig 
	   << " sig1= " << Sigma << " mfp= " << mean_free_path << G4endl;

  } 
  if (nstep > 20000) {
    exit(1);
  }
  */
  /*G4cout<<"Sigma  "<<Sigma<<G4endl;
  G4cout<<"mean_free_path [mm] "<<mean_free_path/mm<<G4endl;
  */
  

  return mean_free_path;
}					 
