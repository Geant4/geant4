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
// $Id: G4AdjointhIonisationModel.cc 66892 2013-01-17 10:57:59Z gunter $
//
#include "G4AdjointhIonisationModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AdjointCSManager.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointProton.hh"
#include "G4AdjointInterpolator.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggModel.hh"
#include "G4Proton.hh"
#include "G4NistManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4AdjointhIonisationModel::G4AdjointhIonisationModel(G4ParticleDefinition* projectileDefinition):
  G4VEmAdjointModel("Adjoint_hIonisation")
{ 



  UseMatrix =true;
  UseMatrixPerElement = true;
  ApplyCutInRange = true;
  UseOnlyOneMatrixForAllElements = true; 
  CS_biasing_factor =1.;
  second_part_of_same_type =false;
  
  //The direct EM Modfel is taken has BetheBloch it is only used for the computation 
  // of the differential cross section.
  //The Bragg model could be used  as an alternative as  it offers the same differential cross section
  
  theDirectEMModel = new G4BetheBlochModel(projectileDefinition);
  theBraggDirectEMModel = new G4BraggModel(projectileDefinition); 
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  
  theDirectPrimaryPartDef = projectileDefinition;
  theAdjEquivOfDirectPrimPartDef = 0;
  if (projectileDefinition == G4Proton::Proton()) {
  	theAdjEquivOfDirectPrimPartDef = G4AdjointProton::AdjointProton();
	
  }
  
  DefineProjectileProperty();
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointhIonisationModel::~G4AdjointhIonisationModel()
{;}


////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointhIonisationModel::SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange)
{ 
 if (!UseMatrix) return RapidSampleSecondaries(aTrack,IsScatProjToProjCase,fParticleChange);
 
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
 G4double projectileKinEnergy = SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy, IsScatProjToProjCase);
 CorrectPostStepWeight(fParticleChange, 
			      aTrack.GetWeight(), 
			      adjointPrimKinEnergy,
			      projectileKinEnergy,
			      IsScatProjToProjCase); //Caution !!!this weight correction should be always applied

 		 
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
void G4AdjointhIonisationModel::RapidSampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 

 const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
 DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());
 
 
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimP =theAdjointPrimary->GetTotalMomentum();
 
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
  
 G4double projectileKinEnergy =0.;
 G4double eEnergy=0.;
 G4double newCS=currentMaterial->GetElectronDensity()*twopi_mc2_rcl2*mass;
 if (!IsScatProjToProjCase){//1/E^2 distribution
 	
	eEnergy=adjointPrimKinEnergy;
	G4double Emax = GetSecondAdjEnergyMaxForProdToProjCase(adjointPrimKinEnergy);
        G4double Emin=  GetSecondAdjEnergyMinForProdToProjCase(adjointPrimKinEnergy);
	if (Emin>=Emax) return;
	G4double a=1./Emax;
	G4double b=1./Emin;
	newCS=newCS*(b-a)/eEnergy;
	
	projectileKinEnergy =1./(b- (b-a)*G4UniformRand()); 
	
 	
 }
 else {	G4double Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(adjointPrimKinEnergy);
	G4double Emin = GetSecondAdjEnergyMinForScatProjToProjCase(adjointPrimKinEnergy,currentTcutForDirectSecond);
	if (Emin>=Emax) return;
	G4double diff1=Emin-adjointPrimKinEnergy;
	G4double diff2=Emax-adjointPrimKinEnergy;
	
	G4double t1=adjointPrimKinEnergy*(1./diff1-1./diff2);
	G4double t2=adjointPrimKinEnergy*(1./Emin-1./Emax);
	/*G4double f31=diff1/Emin;
	G4double f32=diff2/Emax/f31;*/
	G4double t3=2.*std::log(Emax/Emin);
	G4double sum_t=t1+t2+t3;
	newCS=newCS*sum_t/adjointPrimKinEnergy/adjointPrimKinEnergy;
	G4double t=G4UniformRand()*sum_t;
	if (t <=t1 ){
		G4double q= G4UniformRand()*t1/adjointPrimKinEnergy ;
		projectileKinEnergy =adjointPrimKinEnergy +1./(1./diff1-q);
		
	}
	else if (t <=t2 )  {
		G4double q= G4UniformRand()*t2/adjointPrimKinEnergy;
		projectileKinEnergy =1./(1./Emin-q);
	}
	else {	
		projectileKinEnergy=Emin*std::pow(Emax/Emin,G4UniformRand());
		
	}
	eEnergy=projectileKinEnergy-adjointPrimKinEnergy;
	
	
 }

 

 G4double diffCS_perAtom_Used=twopi_mc2_rcl2*mass*adjointPrimKinEnergy/projectileKinEnergy/projectileKinEnergy/eEnergy/eEnergy; 
  
  
						  	
 //Weight correction
 //-----------------------
 //First w_corr is set to the ratio between adjoint total CS and fwd total CS
 G4double w_corr=G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
 
 //G4cout<<w_corr<<G4endl;
 w_corr*=newCS/lastCS;
 //G4cout<<w_corr<<G4endl;
 //Then another correction is needed due to the fact that a biaised differential CS has been used rather than the one consistent with the direct model
 //Here we consider the true  diffCS as the one obtained by the numerical differentiation over Tcut of the direct CS
 
 G4double diffCS = DiffCrossSectionPerAtomPrimToSecond(projectileKinEnergy, eEnergy,1,1);
 w_corr*=diffCS/diffCS_perAtom_Used;
 //G4cout<<w_corr<<G4endl;
	   
 G4double new_weight = aTrack.GetWeight()*w_corr;
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);
 
 
 
 
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
G4double G4AdjointhIonisationModel::DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyProd,
				      G4double Z, 
                                      G4double A)
{//Probably that here the Bragg Model should be also used for kinEnergyProj/nuc<2MeV
  

 
 G4double dSigmadEprod=0;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){ //the produced particle should have a kinetic energy smaller than the projectile 
	G4double Tmax=kinEnergyProj;
	
        G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.000001;
 	G4double dE=(E2-E1);
	G4double sigma1,sigma2;
 	if (kinEnergyProj >2.*MeV){
	    sigma1=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E1,1.e20);
	    sigma2=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E2,1.e20);
        }
	else {
	    sigma1=theBraggDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E1,1.e20);
	    sigma2=theBraggDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E2,1.e20);
	}
	
 	
	dSigmadEprod=(sigma1-sigma2)/dE;
	if (dSigmadEprod>1.) {
		G4cout<<"sigma1 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma1<<G4endl;
		G4cout<<"sigma2 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma2<<G4endl;
		G4cout<<"dsigma "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<dSigmadEprod<<G4endl;
		
	}
	
	 
	
	 //correction of differential cross section at high energy to correct for the suppression of particle at secondary at high
	 //energy used in the Bethe Bloch Model. This correction consist to multiply by g the probability function used
	 //to test the rejection of a secondary
	 //-------------------------
	 
	 //Source code taken from    G4BetheBlochModel::SampleSecondaries
	 
	 G4double deltaKinEnergy = kinEnergyProd;
	 
	 //Part of the taken code
	 //----------------------
	 
	 
	 
	 // projectile formfactor - suppresion of high energy
         // delta-electron production at high energy
	 G4double x = formfact*deltaKinEnergy;
          if(x > 1.e-6) {
		
		
		G4double totEnergy     = kinEnergyProj + mass;
         	G4double etot2         = totEnergy*totEnergy;
	 	G4double beta2 = kinEnergyProj*(kinEnergyProj + 2.0*mass)/etot2;
	 	G4double f;
	 	G4double f1 = 0.0;
		f = 1.0 - beta2*deltaKinEnergy/Tmax;
    	 	if( 0.5 == spin ) {
      			f1 = 0.5*deltaKinEnergy*deltaKinEnergy/etot2;
      			f += f1;
    	 	}
		G4double x1 = 1.0 + x;
    		G4double gg  = 1.0/(x1*x1);
    		if( 0.5 == spin ) {
      			G4double x2 = 0.5*electron_mass_c2*deltaKinEnergy/(mass*mass);
      			gg *= (1.0 + magMoment2*(x2 - f1/f)/(1.0 + x2));
    		}
    		if(gg > 1.0) {
      			G4cout << "### G4BetheBlochModel in Adjoint Sim WARNING: g= " << g
	     			<< G4endl;
			gg=1.;
		}
		//G4cout<<"gg"<<gg<<G4endl;
		dSigmadEprod*=gg;		
    	 }
	 
  }

 return dSigmadEprod;	
}



//////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointhIonisationModel::DefineProjectileProperty()
{   
    //Slightly modified code taken from G4BetheBlochModel::SetParticle
    //------------------------------------------------
    G4String pname = theDirectPrimaryPartDef->GetParticleName();
    if (theDirectPrimaryPartDef->GetParticleType() == "nucleus" &&
	pname != "deuteron" && pname != "triton") {
      isIon = true;
    }
    
    mass = theDirectPrimaryPartDef->GetPDGMass();
    mass_ratio_projectile = proton_mass_c2/theDirectPrimaryPartDef->GetPDGMass();;
    spin = theDirectPrimaryPartDef->GetPDGSpin();
    G4double q = theDirectPrimaryPartDef->GetPDGCharge()/eplus;
    chargeSquare = q*q;
    ratio = electron_mass_c2/mass;
    ratio2 = ratio*ratio;
    one_plus_ratio_2=(1+ratio)*(1+ratio);
    one_minus_ratio_2=(1-ratio)*(1-ratio);
    G4double magmom = theDirectPrimaryPartDef->GetPDGMagneticMoment()
      *mass/(0.5*eplus*hbar_Planck*c_squared);
    magMoment2 = magmom*magmom - 1.0;
    formfact = 0.0;
    if(theDirectPrimaryPartDef->GetLeptonNumber() == 0) {
      G4double x = 0.8426*GeV;
      if(spin == 0.0 && mass < GeV) {x = 0.736*GeV;}
      else if(mass > GeV) {
        x /= G4NistManager::Instance()->GetZ13(mass/proton_mass_c2);
	//	tlimit = 51.2*GeV*A13[iz]*A13[iz];
      }
      formfact = 2.0*electron_mass_c2/(x*x);
      tlimit   = 2.0/formfact;
   }
}

////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointhIonisationModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{ 
  if (UseMatrix) return G4VEmAdjointModel::AdjointCrossSection(aCouple,primEnergy,IsScatProjToProjCase);
  DefineCurrentMaterial(aCouple);
  	
  
  G4double Cross=currentMaterial->GetElectronDensity()*twopi_mc2_rcl2*mass;
  
  if (!IsScatProjToProjCase ){
  	G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(primEnergy);
  	G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(primEnergy);
	if (Emax_proj>Emin_proj && primEnergy > currentTcutForDirectSecond) {
		Cross*=(1./Emin_proj -1./Emax_proj)/primEnergy;
	} 
	else Cross=0.;
		
	

	
	
	
  }
  else {
  	G4double Emax_proj = GetSecondAdjEnergyMaxForScatProjToProjCase(primEnergy);
	G4double Emin_proj = GetSecondAdjEnergyMinForScatProjToProjCase(primEnergy,currentTcutForDirectSecond);
	G4double diff1=Emin_proj-primEnergy;
	G4double diff2=Emax_proj-primEnergy;
	G4double t1=(1./diff1+1./Emin_proj-1./diff2-1./Emax_proj)/primEnergy;
	//G4double t2=2.*std::log(diff2*Emin_proj/Emax_proj/diff1)/primEnergy/primEnergy;
	G4double t2=2.*std::log(Emax_proj/Emin_proj)/primEnergy/primEnergy;
	Cross*=(t1+t2);

  	
  }
  lastCS =Cross;
  return Cross;	
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy)
{ 
  G4double Tmax=PrimAdjEnergy*one_plus_ratio_2/(one_minus_ratio_2-2.*ratio*PrimAdjEnergy/mass);
  return Tmax;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMinForScatProjToProjCase(G4double PrimAdjEnergy,G4double Tcut)
{ return PrimAdjEnergy+Tcut;
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMaxForProdToProjCase(G4double )
{ return HighEnergyLimit;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointhIonisationModel::GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy)
{  G4double Tmin= (2*PrimAdjEnergy-4*mass + std::sqrt(4.*PrimAdjEnergy*PrimAdjEnergy +16.*mass*mass + 8.*PrimAdjEnergy*mass*(1/ratio +ratio)))/4.;
   return Tmin;
}
