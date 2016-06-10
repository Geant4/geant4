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
// $Id: G4AdjointeIonisationModel.cc 66892 2013-01-17 10:57:59Z gunter $
//
#include "G4AdjointeIonisationModel.hh"
#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4Gamma.hh"
#include "G4AdjointGamma.hh"


////////////////////////////////////////////////////////////////////////////////
//
G4AdjointeIonisationModel::G4AdjointeIonisationModel():
 G4VEmAdjointModel("Inv_eIon_model")

{
  
  UseMatrix =true;
  UseMatrixPerElement = true;
  ApplyCutInRange = true;
  UseOnlyOneMatrixForAllElements = true;
  CS_biasing_factor =1.;
  WithRapidSampling = false; 
  
  theAdjEquivOfDirectPrimPartDef =G4AdjointElectron::AdjointElectron();
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  theDirectPrimaryPartDef=G4Electron::Electron();
  second_part_of_same_type=true;
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointeIonisationModel::~G4AdjointeIonisationModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointeIonisationModel::SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange)
{ 


 const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
  
 //Elastic inverse scattering 
 //---------------------------------------------------------
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimP =theAdjointPrimary->GetTotalMomentum();

 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
 
 //Sample secondary energy
 //-----------------------
 G4double projectileKinEnergy;
 if (!WithRapidSampling ) { //used by default
 	projectileKinEnergy = SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy, IsScatProjToProjCase);
	
	CorrectPostStepWeight(fParticleChange, 
			      aTrack.GetWeight(), 
			      adjointPrimKinEnergy,
			      projectileKinEnergy,
			      IsScatProjToProjCase); //Caution !!!this weight correction should be always applied
 }
  else { //only  for test at the moment
  	
	G4double Emin,Emax;
	if (IsScatProjToProjCase) {
		Emin=GetSecondAdjEnergyMinForScatProjToProjCase(adjointPrimKinEnergy,currentTcutForDirectSecond);
		Emax=GetSecondAdjEnergyMaxForScatProjToProjCase(adjointPrimKinEnergy);
	}
	else {
		Emin=GetSecondAdjEnergyMinForProdToProjCase(adjointPrimKinEnergy);
		Emax=GetSecondAdjEnergyMaxForProdToProjCase(adjointPrimKinEnergy);
	}
	projectileKinEnergy = Emin*std::pow(Emax/Emin,G4UniformRand());
	
	
	
	lastCS=lastAdjointCSForScatProjToProjCase;
  	if ( !IsScatProjToProjCase) lastCS=lastAdjointCSForProdToProjCase;
	
	G4double new_weight=aTrack.GetWeight();
	G4double used_diffCS=lastCS*std::log(Emax/Emin)/projectileKinEnergy;
	G4double needed_diffCS=adjointPrimKinEnergy/projectileKinEnergy;
	if (!IsScatProjToProjCase) needed_diffCS *=DiffCrossSectionPerVolumePrimToSecond(currentMaterial,projectileKinEnergy,adjointPrimKinEnergy);
	else needed_diffCS *=DiffCrossSectionPerVolumePrimToScatPrim(currentMaterial,projectileKinEnergy,adjointPrimKinEnergy);
	new_weight*=needed_diffCS/used_diffCS;
	fParticleChange->SetParentWeightByProcess(false);
 	fParticleChange->SetSecondaryWeightByProcess(false);
 	fParticleChange->ProposeParentWeight(new_weight);
	
  
 }	
					  
 
 		 
 //Kinematic: 
 //we consider a two body elastic scattering for the forward processes where the projectile knock on an e- at rest  and gives
 // him part of its  energy
 //----------------------------------------------------------------------------------------
 
 G4double projectileM0 = theAdjEquivOfDirectPrimPartDef->GetPDGMass();
 G4double projectileTotalEnergy = projectileM0+projectileKinEnergy;
 G4double projectileP2 = projectileTotalEnergy*projectileTotalEnergy - projectileM0*projectileM0;	
 
 
 
 //Companion
 //-----------
 G4double companionM0 = theAdjEquivOfDirectPrimPartDef->GetPDGMass();
 if (IsScatProjToProjCase) {
  	companionM0=theAdjEquivOfDirectSecondPartDef->GetPDGMass();
 } 
 G4double companionTotalEnergy =companionM0+ projectileKinEnergy-adjointPrimKinEnergy;
 G4double companionP2 = companionTotalEnergy*companionTotalEnergy - companionM0*companionM0;	
 
 
 //Projectile momentum
 //--------------------
 G4double  P_parallel = (adjointPrimP*adjointPrimP +  projectileP2 - companionP2)/(2.*adjointPrimP);
 G4double  P_perp = std::sqrt( projectileP2 -  P_parallel*P_parallel);
 G4ThreeVector dir_parallel=theAdjointPrimary->GetMomentumDirection();
 G4double phi =G4UniformRand()*2.*3.1415926;
 G4ThreeVector projectileMomentum = G4ThreeVector(P_perp*std::cos(phi),P_perp*std::sin(phi),P_parallel);
 projectileMomentum.rotateUz(dir_parallel);
 
 
 
 if (!IsScatProjToProjCase ){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,projectileMomentum));
	//G4cout<<"projectileMomentum "<<projectileMomentum<<G4endl;
 }
 else {
 	fParticleChange->ProposeEnergy(projectileKinEnergy);
	fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
 }	
     	

  

}
////////////////////////////////////////////////////////////////////////////////
//
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine  
G4double G4AdjointeIonisationModel::DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyProd,
				      G4double Z, 
                                      G4double )
{
 G4double dSigmadEprod=0;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){ //the produced particle should have a kinetic energy smaller than the projectile 
	dSigmadEprod=Z*DiffCrossSectionMoller(kinEnergyProj,kinEnergyProd);
 }
 return dSigmadEprod;	
 
 
 
}

//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointeIonisationModel::DiffCrossSectionMoller(G4double kinEnergyProj,G4double kinEnergyProd){

  G4double energy = kinEnergyProj + electron_mass_c2;
  G4double x   = kinEnergyProd/kinEnergyProj;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  
  G4double gg = (2.0*gam - 1.0)/gamma2;
  G4double y = 1.0 - x;
  G4double fac=twopi_mc2_rcl2/electron_mass_c2;
  G4double dCS = fac*( 1.-gg + ((1.0 - gg*x)/(x*x)) 
		       + ((1.0 - gg*y)/(y*y)))/(beta2*(gam-1));
  return dCS/kinEnergyProj;
  
 

}  

