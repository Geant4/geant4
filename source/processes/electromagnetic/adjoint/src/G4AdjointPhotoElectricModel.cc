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
#include "G4AdjointPhotoElectricModel.hh"
#include "G4AdjointCSManager.hh"


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
  current_eEnergy =0.;
  totAdjointCS=0.;
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
  const G4MaterialCutsCouple* aCouple = aTrack.GetMaterialCutsCouple();
  const G4DynamicParticle* aDynPart =  aTrack.GetDynamicParticle() ;
  G4double electronEnergy = aDynPart->GetKineticEnergy();
  G4ThreeVector electronDirection= aDynPart->GetMomentumDirection() ;
  totAdjointCS = AdjointCrossSection(aCouple, electronEnergy,IsScatProjToProjCase);
				
 
  //Sample gamma energy
  //-------------
  	/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4ContinuousGainOfEnergy.hh
//	Author:       	L. Desorgher
//	Date:		1 September 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1 September 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Modell for the adjoint compton scattering
//

  //Sample element
  //-------------
   const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
   const G4double* theAtomNumDensityVector = currentMaterial->GetVecNbOfAtomsPerVolume();
   size_t nelm =  currentMaterial->GetNumberOfElements();
   G4double rand_CS= totAdjointCS*G4UniformRand();
   for (index_element=0; index_element<nelm-1; index_element++){
	if (rand_CS<xsec[index_element]) break;
   }
	
   //Sample shell and binding energy
   //-------------
   rand_CS= totAdjointCS*G4UniformRand()/theAtomNumDensityVector[index_element];
   G4int nShells = (*theElementVector)[index_element]->GetNbOfAtomicShells();
   G4int i  = 0;  
   for (i=0; i<nShells-1; i++){
	if (rand_CS<shell_prob[index_element][i]) break;
   }
   G4double gammaEnergy= electronEnergy+(*theElementVector)[index_element]->GetAtomicShell(i);
	
  //Sample cos theta
  //Copy of the G4PEEffectModel cos theta sampling method ElecCosThetaDistribution.   
  //This method cannot be used directly from G4PEEffectModel because it is a friend method. I should ask Vladimir to change that  
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
  CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), electronEnergy,gammaEnergy);	
 
  
  
  //Create secondary and modify fParticleChange 
  //--------------------------------------------
  G4DynamicParticle* anAdjointGamma = new G4DynamicParticle (
                       G4AdjointGamma::AdjointGamma(),adjoint_gammaDirection, gammaEnergy);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->AddSecondary(anAdjointGamma);
     
     	

  

}	

////////////////////////////////////////////////////////////////////////////////
//			

G4double G4AdjointPhotoElectricModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double electronEnergy,
				G4bool IsScatProjToProjCase)
{ if (IsScatProjToProjCase) return 0.;
  if (aCouple !=currentCouple || current_eEnergy !=electronEnergy) {
  	totAdjointCS = 0.;
	DefineCurrentMaterialAndElectronEnergy(aCouple, electronEnergy);
  	const G4ElementVector* theElementVector = currentMaterial->GetElementVector();
  	const G4double* theAtomNumDensityVector = currentMaterial->GetVecNbOfAtomsPerVolume();
  	size_t nelm =  currentMaterial->GetNumberOfElements();
  	for (index_element=0;index_element<nelm;index_element++){
		
		totAdjointCS +=AdjointCrossSectionPerAtom((*theElementVector)[index_element],electronEnergy)*theAtomNumDensityVector[index_element];
		xsec[index_element] = totAdjointCS;
	} 
  }	
  return totAdjointCS;
  
  
}
////////////////////////////////////////////////////////////////////////////////
//			

G4double G4AdjointPhotoElectricModel::AdjointCrossSectionPerAtom(const G4Element* anElement,G4double electronEnergy)
{ 
  G4int nShells = anElement->GetNbOfAtomicShells();
  G4double Z= anElement->GetZ();
  G4double N= anElement->GetN();
  G4int i  = 0;  
  G4double B0=anElement->GetAtomicShell(0);
  G4double gammaEnergy = electronEnergy+B0;
  G4double adjointCS = theDirectPEEffectModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),gammaEnergy,Z,N,0.,0.)*electronEnergy/gammaEnergy; 
  shell_prob[index_element][0] = adjointCS;                                          
  for (i=1;i<nShells;i++){
  	//G4cout<<i<<std::endl;
  	G4double Bi_= anElement->GetAtomicShell(i-1);
	G4double Bi = anElement->GetAtomicShell(i);
	//G4cout<<Bi_<<'\t'<<Bi<<std::endl;
	if (electronEnergy <Bi_-Bi) {
		gammaEnergy = electronEnergy+Bi;
		adjointCS +=theDirectPEEffectModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),gammaEnergy,anElement->GetZ(),N,0.,0.)*electronEnergy/gammaEnergy;
	}
	shell_prob[index_element][i] = adjointCS;	
  
  }
  
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
}
