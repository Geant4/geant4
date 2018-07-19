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
// $Id: G4AdjointComptonModel.cc 100666 2016-10-31 10:27:00Z gcosmo $
//
#include "G4AdjointComptonModel.hh"
#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Gamma.hh"
#include "G4KleinNishinaCompton.hh"


////////////////////////////////////////////////////////////////////////////////
//
G4AdjointComptonModel::G4AdjointComptonModel():
 G4VEmAdjointModel("AdjointCompton")

{ SetApplyCutInRange(false);
  SetUseMatrix(false);
  SetUseMatrixPerElement(true);
  SetUseOnlyOneMatrixForAllElements(true);
  theAdjEquivOfDirectPrimPartDef =G4AdjointGamma::AdjointGamma();
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  theDirectPrimaryPartDef=G4Gamma::Gamma();
  second_part_of_same_type=false;
  theDirectEMModel=new G4KleinNishinaCompton(G4Gamma::Gamma(),"ComptonDirectModel");
  G4direct_CS = 0.;
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
   if (!UseMatrix) return RapidSampleSecondaries(aTrack,IsScatProjToProjCase,fParticleChange); 
   
   //A recall of the compton scattering law is 
   //Egamma2=Egamma1/(1+(Egamma1/E0_electron)(1.-cos_th))
   //Therefore Egamma2_max= Egamma2(cos_th=1) = Egamma1
   //Therefore Egamma2_min= Egamma2(cos_th=-1) = Egamma1/(1+2.(Egamma1/E0_electron))
   
   
  const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
  G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
  if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
  }
 
 
 //Sample secondary energy
 //-----------------------
 G4double gammaE1;
 gammaE1 = SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy,
						  IsScatProjToProjCase);
 
 
 //gammaE2
 //-----------
 
 G4double gammaE2 = adjointPrimKinEnergy;
 if (!IsScatProjToProjCase) gammaE2 = gammaE1 - adjointPrimKinEnergy;	
 
 
 
 
 
 
 //Cos th
 //-------
// G4cout<<"Compton scattering "<<gammaE1<<'\t'<<gammaE2<<G4endl;
 G4double cos_th = 1.+ electron_mass_c2*(1./gammaE1 -1./gammaE2);
 if (!IsScatProjToProjCase) {
 	G4double p_elec=theAdjointPrimary->GetTotalMomentum();
	cos_th = (gammaE1 - gammaE2*cos_th)/p_elec;
 }
 G4double sin_th = 0.;
 if (std::abs(cos_th)>1){
 	//G4cout<<"Problem in compton scattering with cos_th "<<cos_th<<G4endl;
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
// G4cout<<gamma0Energy<<'\t'<<gamma0Momentum<<G4endl;
 
 
 //It is important to correct the weight of particles before adding the secondary
 //------------------------------------------------------------------------------
 CorrectPostStepWeight(fParticleChange, 
 			aTrack.GetWeight(), 
			adjointPrimKinEnergy,
			gammaE1,
			IsScatProjToProjCase);
 
 if (!IsScatProjToProjCase){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,gammaMomentum1));
	//G4cout<<"gamma0Momentum "<<gamma0Momentum<<G4endl;
 }
 else {
 	fParticleChange->ProposeEnergy(gammaE1);
	fParticleChange->ProposeMomentumDirection(gammaMomentum1.unit());
 }
 
 	
} 
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointComptonModel::RapidSampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 

 const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
 DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());
 
 
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();

 
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
  
 
 
 G4double diffCSUsed=0.1*currentMaterial->GetElectronDensity()*twopi_mc2_rcl2;
 G4double gammaE1=0.;
 G4double gammaE2=0.;
 if (!IsScatProjToProjCase){
 	
	G4double Emax = GetSecondAdjEnergyMaxForProdToProjCase(adjointPrimKinEnergy);
        G4double Emin=  GetSecondAdjEnergyMinForProdToProjCase(adjointPrimKinEnergy);;
	if (Emin>=Emax) return;
	G4double f1=(Emin-adjointPrimKinEnergy)/Emin;
	G4double f2=(Emax-adjointPrimKinEnergy)/Emax/f1;
	gammaE1=adjointPrimKinEnergy/(1.-f1*std::pow(f2,G4UniformRand()));;
	gammaE2=gammaE1-adjointPrimKinEnergy;
	diffCSUsed= diffCSUsed*(1.+2.*std::log(1.+electron_mass_c2/adjointPrimKinEnergy))*adjointPrimKinEnergy/gammaE1/gammaE2;
	
 	
 }
 else {	G4double Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(adjointPrimKinEnergy);
	G4double Emin = GetSecondAdjEnergyMinForScatProjToProjCase(adjointPrimKinEnergy,currentTcutForDirectSecond);
	if (Emin>=Emax) return;
	gammaE2 =adjointPrimKinEnergy;
	gammaE1=Emin*std::pow(Emax/Emin,G4UniformRand());
	diffCSUsed= diffCSUsed/gammaE1;
 }
  
  
  
						  	
 //Weight correction
 //-----------------------
 //First w_corr is set to the ratio between adjoint total CS and fwd total CS
 G4double w_corr=additional_weight_correction_factor_for_post_step_outside_model;
  if (correct_weight_for_post_step_in_model) {
      w_corr=G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
  }
 //Then another correction is needed due to the fact that a biaised differential CS has been used rather than the 
 //one consistent with the direct model
 
 
 G4double diffCS = DiffCrossSectionPerAtomPrimToScatPrim(gammaE1, gammaE2,1,0.);
 if (diffCS >0)  diffCS /=G4direct_CS;  // here we have the normalised diffCS
 //An we remultiply by the lambda of the forward process
 diffCS*=theDirectEMProcess->GetLambda(gammaE1,currentCouple);
 //diffCS*=theDirectEMModel->CrossSectionPerVolume(currentMaterial,G4Gamma::Gamma(),gammaE1,0.,2.*gammaE1);
 //G4cout<<"diffCS/diffCSUsed "<<diffCS/diffCSUsed<<'\t'<<gammaE1<<'\t'<<gammaE2<<G4endl;                                 
 
 w_corr*=diffCS/diffCSUsed;
	   
 G4double new_weight = aTrack.GetWeight()*w_corr;
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight); 
 
 
 
 //Cos th
 //-------

 G4double cos_th = 1.+ electron_mass_c2*(1./gammaE1 -1./gammaE2);
 if (!IsScatProjToProjCase) {
 	G4double p_elec=theAdjointPrimary->GetTotalMomentum();
	cos_th = (gammaE1 - gammaE2*cos_th)/p_elec;
 }
 G4double sin_th = 0.;
 if (std::abs(cos_th)>1){
 	//G4cout<<"Problem in compton scattering with cos_th "<<cos_th<<G4endl;
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
 
 

 
 if (!IsScatProjToProjCase){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,gammaMomentum1));
	//G4cout<<"gamma0Momentum "<<gamma0Momentum<<G4endl;
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
 // In the forward case (see G4KleinNishinaCompton)  the  cross section is parametrised
 // but the secondaries are sampled from the
 // Klein Nishida differential cross section
 // The used diffrential cross section here is therefore the cross section multiplied by the normalised 
 //differential Klein Nishida cross section
 
 
 //Klein Nishida Cross Section
 //-----------------------------
 G4double epsilon = gamEnergy0 / electron_mass_c2 ;
 G4double one_plus_two_epsi =1.+2.*epsilon;
 G4double gamEnergy1_max = gamEnergy0;
 G4double gamEnergy1_min = gamEnergy0/one_plus_two_epsi;
 if (gamEnergy1 >gamEnergy1_max ||  gamEnergy1<gamEnergy1_min) {
 	/*G4cout<<"the differential CS is null"<<G4endl;
	G4cout<<gamEnergy0<<G4endl;
	G4cout<<gamEnergy1<<G4endl;
	G4cout<<gamEnergy1_min<<G4endl;*/
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
 
 G4direct_CS = theDirectEMModel->ComputeCrossSectionPerAtom(G4Gamma::Gamma(),
                                             gamEnergy0,
                                             Z, 0., 0.,0.);
 
 dCS_dE1 *= G4direct_CS/CS;
/* G4cout<<"the differential CS is not null"<<G4endl;
 G4cout<<gamEnergy0<<G4endl;
 G4cout<<gamEnergy1<<G4endl;*/
 
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
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointComptonModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{ 
  if (UseMatrix) return G4VEmAdjointModel::AdjointCrossSection(aCouple,primEnergy,IsScatProjToProjCase);
  DefineCurrentMaterial(aCouple);
  
  
  float Cross=0.;
  float Emax_proj =0.;
  float Emin_proj =0.;
  if (!IsScatProjToProjCase ){
  	Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(primEnergy);
  	Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(primEnergy);
	if (Emax_proj>Emin_proj ){
		 Cross= 0.1*std::log((Emax_proj-float (primEnergy))*Emin_proj/Emax_proj/(Emin_proj-primEnergy))
		 						*(1.+2.*std::log(float(1.+electron_mass_c2/primEnergy)));
	}	 
  }
  else {
        Emax_proj = GetSecondAdjEnergyMaxForScatProjToProjCase(primEnergy);
	Emin_proj = GetSecondAdjEnergyMinForScatProjToProjCase(primEnergy,0.);
	if (Emax_proj>Emin_proj) {
		Cross = 0.1*std::log(Emax_proj/Emin_proj);
		//+0.5*primEnergy*primEnergy(1./(Emin_proj*Emin_proj) - 1./(Emax_proj*Emax_proj)); neglected at the moment
	}
  	
  	
  }
  
  Cross*=currentMaterial->GetElectronDensity()*twopi_mc2_rcl2;
  lastCS=Cross;
  return double(Cross);	
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointComptonModel::GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{ return AdjointCrossSection(aCouple, primEnergy,IsScatProjToProjCase);
  //return G4VEmAdjointModel::GetAdjointCrossSection(aCouple, primEnergy,IsScatProjToProjCase);
}
