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
#include "G4VAdjointInverseScattering.hh"
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


G4VAdjointInverseScattering::
	G4VAdjointInverseScattering(G4String process_name, G4bool whichScatCase):
			G4VDiscreteProcess(process_name)
{theAdjointCSManager = G4AdjointCSManager::GetAdjointCSManager();
 IsScatProjToProjCase=whichScatCase;
 /*theAdjointEMModel=aModel;
 IsScatProjToProjCase=whichScatCase;*/
 fParticleChange=new G4ParticleChange();
 
}
//////////////////////////////////////////////////////////////////////////////
//
G4VAdjointInverseScattering::
	~G4VAdjointInverseScattering()
{;
}			
//////////////////////////////////////////////////////////////////////////////
//
void G4VAdjointInverseScattering::PreparePhysicsTable(const G4ParticleDefinition&)
{;
}
//////////////////////////////////////////////////////////////////////////////
//
void G4VAdjointInverseScattering::BuildPhysicsTable(const G4ParticleDefinition&)
{

 theAdjointCSManager->BuildCrossSectionMatrices(); //do not worry it will be done just once
 theAdjointCSManager->BuildTotalSigmaTables();

}
//////////////////////////////////////////////////////////////////////////////
//
G4VParticleChange* G4VAdjointInverseScattering::PostStepDoIt(const G4Track& track, const G4Step& )
{ 
  fParticleChange->Initialize(track);
 
  theAdjointEMModel->SampleSecondaries(track,
                                       IsScatProjToProjCase,
					fParticleChange);
  
  ClearNumberOfInteractionLengthLeft();
  return fParticleChange;
  			
   
  
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4VAdjointInverseScattering::GetMeanFreePath(const G4Track& track,
                                         					     G4double ,
                                         					     G4ForceCondition* condition)
{ *condition = NotForced;
  G4double preStepKinEnergy = track.GetKineticEnergy();
  
  G4double Sigma =
  		theAdjointEMModel->AdjointCrossSection(track.GetMaterialCutsCouple(),preStepKinEnergy,IsScatProjToProjCase);
  		

  G4double mean_free_path = 1./Sigma;
  //G4cout<<"mean_free_path [mm] "<<mean_free_path/mm<<std::endl;
  return mean_free_path;
}					 
				 
