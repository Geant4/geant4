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
// $Id: G4AdjointPhotoElectricModel.cc 91870 2015-08-07 15:21:40Z gcosmo $
//
#include "G4AdjointPhotoElectricModel.hh"
#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include  "G4Gamma.hh"
#include "G4AdjointGamma.hh"


////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPhotoElectricModel::G4AdjointPhotoElectricModel():
 G4VEmAdjointModel("AdjointPEEffect")

{ SetUseMatrix(false);
  SetApplyCutInRange(false);

  //Initialization
  current_eEnergy =0.;
  totAdjointCS=0.;
  factorCSBiasing =1.;
  post_step_AdjointCS =0.;
  pre_step_AdjointCS =0.;
  totBiasedAdjointCS =0.;

  index_element=0;

  theAdjEquivOfDirectPrimPartDef =G4AdjointGamma::AdjointGamma();
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  theDirectPrimaryPartDef=G4Gamma::Gamma();
  second_part_of_same_type=false;
  theDirectPEEffectModel = new G4PEEffectFluoModel();
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointPhotoElectricModel::~G4AdjointPhotoElectricModel()
{;}

////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPhotoElectricModel::SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange)
{ if (IsScatProjToProjCase) return ;

  //Compute the totAdjointCS vectors if not already done for the current couple and electron energy
  //-----------------------------------------------------------------------------------------------
  const G4MaterialCutsCouple* aCouple = aTrack.GetMaterialCutsCouple();
  const G4DynamicParticle* aDynPart =  aTrack.GetDynamicParticle() ;
  G4double electronEnergy = aDynPart->GetKineticEnergy();
  G4ThreeVector electronDirection= aDynPart->GetMomentumDirection() ;
  pre_step_AdjointCS = totAdjointCS; //The last computed CS was  at pre step point
  post_step_AdjointCS = AdjointCrossSection(aCouple, electronEnergy,IsScatProjToProjCase);
  post_step_AdjointCS = totAdjointCS; 
				



  //Sample element
  //-------------
   const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
   size_t nelm =  currentMaterial->GetNumberOfElements();
   G4double rand_CS= G4UniformRand()*xsec[nelm-1];
   for (index_element=0; index_element<nelm-1; index_element++){
	if (rand_CS<xsec[index_element]) break;
   }
	
   //Sample shell and binding energy
   //-------------
   G4int nShells = (*theElementVector)[index_element]->GetNbOfAtomicShells();
   rand_CS= shell_prob[index_element][nShells-1]*G4UniformRand();
   G4int i  = 0;  
   for (i=0; i<nShells-1; i++){
	if (rand_CS<shell_prob[index_element][i]) break;
   }
   G4double gammaEnergy= electronEnergy+(*theElementVector)[index_element]->GetAtomicShell(i);
	
  //Sample cos theta
  //Copy of the G4PEEfectFluoModel cos theta sampling method ElecCosThetaDistribution.   
  //This method cannot be used directly from G4PEEfectFluoModel because it is a friend method. I should ask Vladimir to change that  
  //------------------------------------------------------------------------------------------------	
  //G4double cos_theta = theDirectPEEffectModel->ElecCosThetaDistribution(electronEnergy);
	
   G4double  cos_theta = 1.;
   G4double gamma   = 1. + electronEnergy/electron_mass_c2;
   if (gamma <= 5.) {
  	G4double beta  = std::sqrt(gamma*gamma-1.)/gamma;
 	G4double b     = 0.5*gamma*(gamma-1.)*(gamma-2);

  	G4double rndm,term,greject,grejsup;
  	if (gamma < 2.) grejsup = gamma*gamma*(1.+b-beta*b);
  	else            grejsup = gamma*gamma*(1.+b+beta*b);

  	do {	rndm = 1.-2*G4UniformRand();
       	 	cos_theta = (rndm+beta)/(rndm*beta+1.);
       		term = 1.-beta*cos_theta;
       		greject = (1.-cos_theta*cos_theta)*(1.+b*term)/(term*term);
          // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  	} while(greject < G4UniformRand()*grejsup);
  }
	
  // direction of the adjoint gamma electron
  //---------------------------------------
  
     
  G4double sin_theta = std::sqrt(1.-cos_theta*cos_theta);
  G4double Phi     = twopi * G4UniformRand();
  G4double dirx = sin_theta*std::cos(Phi),diry = sin_theta*std::sin(Phi),dirz = cos_theta;
  G4ThreeVector adjoint_gammaDirection(dirx,diry,dirz);
  adjoint_gammaDirection.rotateUz(electronDirection);
  
  
  
  //Weight correction
 //-----------------------					   
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), electronEnergy,gammaEnergy,IsScatProjToProjCase);	
 
  
  
  //Create secondary and modify fParticleChange 
  //--------------------------------------------
  G4DynamicParticle* anAdjointGamma = new G4DynamicParticle (
                       G4AdjointGamma::AdjointGamma(),adjoint_gammaDirection, gammaEnergy);
  
  
  
 
  
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->AddSecondary(anAdjointGamma);    
     	

  

}

////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointPhotoElectricModel::CorrectPostStepWeight(G4ParticleChange* fParticleChange, 
							    G4double old_weight,  
							    G4double adjointPrimKinEnergy, 
							    G4double projectileKinEnergy ,
							    G4bool  ) 
{
 G4double new_weight=old_weight;

 G4double w_corr =G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection()/factorCSBiasing;
 w_corr*=post_step_AdjointCS/pre_step_AdjointCS; 


 new_weight*=w_corr; 
 new_weight*=projectileKinEnergy/adjointPrimKinEnergy;
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);	
}	

////////////////////////////////////////////////////////////////////////////////
//			

G4double G4AdjointPhotoElectricModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double electronEnergy,
				G4bool IsScatProjToProjCase)
{ 
  

  if (IsScatProjToProjCase)  return 0.;

  	
  if (aCouple !=currentCouple || current_eEnergy !=electronEnergy) {
  	totAdjointCS = 0.;
	DefineCurrentMaterialAndElectronEnergy(aCouple, electronEnergy);
  	const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
  	const double* theAtomNumDensityVector = currentMaterial->GetVecNbOfAtomsPerVolume();
  	size_t nelm =  currentMaterial->GetNumberOfElements();
  	for (index_element=0;index_element<nelm;index_element++){
		
		totAdjointCS +=AdjointCrossSectionPerAtom((*theElementVector)[index_element],electronEnergy)*theAtomNumDensityVector[index_element];
		xsec[index_element] = totAdjointCS;
	} 

	totBiasedAdjointCS=std::min(totAdjointCS,0.01);
//	totBiasedAdjointCS=totAdjointCS;
	factorCSBiasing = totBiasedAdjointCS/totAdjointCS;
	lastCS=totBiasedAdjointCS;
  	
	
  }
  return totBiasedAdjointCS;

  
}
////////////////////////////////////////////////////////////////////////////////
//			

G4double G4AdjointPhotoElectricModel::GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double electronEnergy,
				G4bool IsScatProjToProjCase)
{	return AdjointCrossSection(aCouple,electronEnergy,IsScatProjToProjCase);
}
////////////////////////////////////////////////////////////////////////////////
//			

G4double G4AdjointPhotoElectricModel::AdjointCrossSectionPerAtom(const G4Element* anElement,G4double electronEnergy)
{ 
  G4int nShells = anElement->GetNbOfAtomicShells();
  G4double Z= anElement->GetZ();
  G4int i  = 0;  
  G4double B0=anElement->GetAtomicShell(0);
  G4double gammaEnergy = electronEnergy+B0;
  G4double CS= theDirectPEEffectModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),gammaEnergy,Z,0.,0.,0.);
  G4double adjointCS =0.;
  if (CS >0) adjointCS += CS/gammaEnergy; 
  shell_prob[index_element][0] = adjointCS;                                          
  for (i=1;i<nShells;i++){
  	//G4cout<<i<<G4endl;
  	G4double Bi_= anElement->GetAtomicShell(i-1);
	G4double Bi = anElement->GetAtomicShell(i);
	//G4cout<<Bi_<<'\t'<<Bi<<G4endl;
	if (electronEnergy <Bi_-Bi) {
		gammaEnergy = electronEnergy+Bi;
		
		CS=theDirectPEEffectModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),gammaEnergy,Z,0.,0.,0.);
		if (CS>0) adjointCS +=CS/gammaEnergy;
	}
	shell_prob[index_element][i] = adjointCS;	
  
  }
  adjointCS*=electronEnergy;
  return adjointCS;
  
}				
////////////////////////////////////////////////////////////////////////////////
//			

void G4AdjointPhotoElectricModel::DefineCurrentMaterialAndElectronEnergy(const G4MaterialCutsCouple* couple, G4double anEnergy)
{ currentCouple   = const_cast<G4MaterialCutsCouple*> (couple);
  currentMaterial = const_cast<G4Material*> (couple->GetMaterial());
  currentCoupleIndex = couple->GetIndex();
  currentMaterialIndex = currentMaterial->GetIndex();
  current_eEnergy = anEnergy;	
  theDirectPEEffectModel->SetCurrentCouple(couple);
}
