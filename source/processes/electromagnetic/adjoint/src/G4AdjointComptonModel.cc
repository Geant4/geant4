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
#include "G4AdjointComptonModel.hh"
#include "G4AdjointCSManager.hh"


#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Gamma.hh"


////////////////////////////////////////////////////////////////////////////////
//
G4AdjointComptonModel::G4AdjointComptonModel():
 G4VEmAdjointModel("AdjointCompton")

{ SetApplyCutInRange(false);
  SetUseMatrix(true);
  SetUseMatrixPerElement(true);
  SetIsIonisation(false);
  SetUseOnlyOneMatrixForAllElements(true);
  theAdjEquivOfDirectPrimPartDef =G4AdjointGamma::AdjointGamma();
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  theDirectPrimaryPartDef=G4Gamma::Gamma();
  second_part_of_same_type=false;
 
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointComptonModel::~G4AdjointComptonModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointComptonModel::SampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 
   
   
   //A recall of the compton scattering law is 
   //Egamma2=Egamma1/(1+(Egamma1/E0_electron)(1.-cos_th))
   //Therefore Egamma2_max= Egamma2(cos_th=1) = Egamma1
   //Therefore Egamma2_min= Egamma2(cos_th=-1) = Egamma1/(1+2.(Egamma1/E0_electron))
   
   
  const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
  //DefineCurrentMaterial(aTrack->GetMaterialCutsCouple());
  size_t ind= 0;
  
  
 
 
 //Elastic inverse scattering //not correct in all the cases 
 //---------------------------------------------------------
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 
 //G4cout<<adjointPrimKinEnergy<<std::endl;
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
 
 //Sample secondary energy
 //-----------------------
 G4double gammaE1;

 gammaE1 = SampleAdjSecEnergyFromCSMatrix(ind,
 						  adjointPrimKinEnergy,
						  IsScatProjToProjCase);
 
 
 //gammaE2
 //-----------
 
 G4double gammaE2 = adjointPrimKinEnergy;
 if (!IsScatProjToProjCase) gammaE2 = gammaE1 - adjointPrimKinEnergy;	
 
 
 
 
 
 
 //Cos th
 //-------
// G4cout<<"Compton scattering "<<gammaE1<<'\t'<<gammaE2<<std::endl;
 G4double cos_th = 1.+ electron_mass_c2*(1./gammaE1 -1./gammaE2);
 if (!IsScatProjToProjCase) {
 	G4double p_elec=theAdjointPrimary->GetTotalMomentum();
	cos_th = (gammaE1 - gammaE2*cos_th)/p_elec;
 }
 G4double sin_th = 0.;
 if (std::abs(cos_th)>1){
 	//G4cout<<"Problem in compton scattering with cos_th "<<cos_th<<std::endl;
	if (cos_th>0) {
		cos_th=1.;
	}
	else 	cos_th=-1.;
	sin_th=0.;
 }
 else  sin_th = std::sqrt(1.-cos_th*cos_th);

 
 
 
 //gamma0 momentum
 //--------------------

 
 G4ThreeVector dir_parallel=theAdjointPrimary->GetMomentumDirection();
 G4double phi =G4UniformRand()*2.*3.1415926;
 G4ThreeVector gammaMomentum1 = gammaE1*G4ThreeVector(std::cos(phi)*sin_th,std::sin(phi)*sin_th,cos_th);
 gammaMomentum1.rotateUz(dir_parallel);
// G4cout<<gamma0Energy<<'\t'<<gamma0Momentum<<std::endl;
 
 
 //It is important to correct the weight of particles before adding the secondary
 //------------------------------------------------------------------------------
 CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), adjointPrimKinEnergy,gammaE1);
 
 if (!IsScatProjToProjCase && CorrectWeightMode){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,gammaMomentum1));
	//G4cout<<"gamma0Momentum "<<gamma0Momentum<<std::endl;
 }
 else {
 	fParticleChange->ProposeEnergy(gammaE1);
	fParticleChange->ProposeMomentumDirection(gammaMomentum1.unit());
 }
 
 	
} 			
////////////////////////////////////////////////////////////////////////////////
//
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine  
G4double G4AdjointComptonModel::DiffCrossSectionPerAtomPrimToSecond(
                                      G4double gamEnergy0, 
                                      G4double kinEnergyElec,
				      G4double Z, 
                                      G4double A)
{ 
  G4double gamEnergy1 =  gamEnergy0 - kinEnergyElec;
  G4double dSigmadEprod=0.;
  if (gamEnergy1>0.) dSigmadEprod=DiffCrossSectionPerAtomPrimToScatPrim(gamEnergy0,gamEnergy1,Z,A);
  return dSigmadEprod;	  
}


////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointComptonModel::DiffCrossSectionPerAtomPrimToScatPrim(
                                      G4double gamEnergy0, 
                                      G4double gamEnergy1,
				      G4double Z, 
                                      G4double )
{ //Based on Klein Nishina formula
 //* In the forward case (see G4KleinNishinaModel)  the  cross section is parametrised while the secondaries are sampled from the 
 // Klein Nishida differential cross section
 // The used diffrential cross section here is therefore the cross section multiplied by the normalidsed differential Klein Nishida cross section
 
 
 //Klein Nishida Cross Section
 //-----------------------------
 G4double epsilon = gamEnergy0 / electron_mass_c2 ;
 G4double one_plus_two_epsi =1.+2.*epsilon;
 G4double gamEnergy1_max = gamEnergy0;
 G4double gamEnergy1_min = gamEnergy0/one_plus_two_epsi;
 if (gamEnergy1 >gamEnergy1_max ||  gamEnergy1<gamEnergy1_min) {
 	/*G4cout<<"the differential CS is null"<<std::endl;
	G4cout<<gamEnergy0<<std::endl;
	G4cout<<gamEnergy1<<std::endl;
	G4cout<<gamEnergy1_min<<std::endl;*/
 	return 0.;
 }
 	
 
 G4double epsi2 = epsilon *epsilon ;
 G4double one_plus_two_epsi_2=one_plus_two_epsi*one_plus_two_epsi;
 
 
 G4double CS=std::log(one_plus_two_epsi)*(1.- 2.*(1.+epsilon)/epsi2);
 CS+=4./epsilon +0.5*(1.-1./one_plus_two_epsi_2);
 CS/=epsilon;
 //Note that the pi*re2*Z factor is neglected because it is suppresed when computing dCS_dE1/CS;
 // in the differential cross section
 
 
 //Klein Nishida Differential Cross Section
 //-----------------------------------------
 G4double epsilon1 = gamEnergy1 / electron_mass_c2 ;
 G4double v= epsilon1/epsilon;
 G4double term1 =1.+ 1./epsilon -1/epsilon1;
 G4double dCS_dE1= 1./v +v + term1*term1 -1.;
 dCS_dE1 *=1./epsilon/gamEnergy0;
 
 
 //Normalised to the CS used in G4
 //-------------------------------
 
 G4double G4direct_CS = theDirectEMModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),
                                             gamEnergy0,
                                             Z, 0., 0.,0.);
 
 dCS_dE1 *= G4direct_CS/CS;
/* G4cout<<"the differential CS is not null"<<std::endl;
 G4cout<<gamEnergy0<<std::endl;
 G4cout<<gamEnergy1<<std::endl;*/
 
 return dCS_dE1;


}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointComptonModel::GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy)
{ G4double inv_e_max =  1./PrimAdjEnergy - 2./electron_mass_c2;
  G4double e_max = HighEnergyLimit;
  if (inv_e_max > 0. ) e_max=std::min(1./inv_e_max,HighEnergyLimit);
  return  e_max; 
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointComptonModel::GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy)
{ G4double half_e=PrimAdjEnergy/2.;
  G4double term=std::sqrt(half_e*(electron_mass_c2+half_e));
  G4double emin=half_e+term;
  return  emin; 
}
