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
// $Id: G4AdjointIonIonisationModel.cc 66892 2013-01-17 10:57:59Z gunter $
//
#include "G4AdjointIonIonisationModel.hh"
#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointProton.hh"
#include "G4AdjointInterpolator.hh"
#include "G4BetheBlochModel.hh"
#include "G4BraggIonModel.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4NistManager.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4AdjointIonIonisationModel::G4AdjointIonIonisationModel():
  G4VEmAdjointModel("Adjoint_IonIonisation")
{ 


  UseMatrix =true;
  UseMatrixPerElement = true;
  ApplyCutInRange = true;
  UseOnlyOneMatrixForAllElements = true;
  CS_biasing_factor =1.;
  second_part_of_same_type =false;
  use_only_bragg = false; // for the Ion ionisation using the parametrised table model the cross sections and the sample of secondaries is done
  				// as in the BraggIonModel, Therefore the use of this flag;  
  
  //The direct EM Model is taken has BetheBloch it is only used for the computation 
  // of the differential cross section.
  //The Bragg model could be used  as an alternative as  it offers the same differential cross section
  
  theBetheBlochDirectEMModel = new G4BetheBlochModel(G4GenericIon::GenericIon());
  theBraggIonDirectEMModel = new G4BraggIonModel(G4GenericIon::GenericIon()); 
  theAdjEquivOfDirectSecondPartDef=G4AdjointElectron::AdjointElectron();
  theDirectPrimaryPartDef =0;
  theAdjEquivOfDirectPrimPartDef =0;
 /* theDirectPrimaryPartDef =fwd_ion;
  theAdjEquivOfDirectPrimPartDef =adj_ion;
 
  DefineProjectileProperty();*/




}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointIonIonisationModel::~G4AdjointIonIonisationModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointIonIonisationModel::SampleSecondaries(const G4Track& aTrack,
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
G4double G4AdjointIonIonisationModel::DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyProd,
				      G4double Z, 
                                      G4double A)
{//Probably that here the Bragg Model should be also used for kinEnergyProj/nuc<2MeV
  

 
 G4double dSigmadEprod=0;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 G4double kinEnergyProjScaled = massRatio*kinEnergyProj;
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){ //the produced particle should have a kinetic energy smaller than the projectile 
	G4double Tmax=kinEnergyProj;
	
        G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.000001;
 	G4double dE=(E2-E1);
	G4double sigma1,sigma2;
	theDirectEMModel =theBraggIonDirectEMModel;
 	if (kinEnergyProjScaled >2.*MeV && !use_only_bragg) theDirectEMModel = theBetheBlochDirectEMModel; //Bethe Bloch Model
	sigma1=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E1,1.e20);
	sigma2=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E2,1.e20);
	
	dSigmadEprod=(sigma1-sigma2)/dE;
	
	//G4double chargeSqRatio =  currentModel->GetChargeSquareRatio(theDirectPrimaryPartDef,currentMaterial,E);
	
	
	
	if (dSigmadEprod>1.) {
		G4cout<<"sigma1 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma1<<G4endl;
		G4cout<<"sigma2 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma2<<G4endl;
		G4cout<<"dsigma "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<dSigmadEprod<<G4endl;
		
	}
	
	 
	
	 
	 
	 
	 if (theDirectEMModel == theBetheBlochDirectEMModel ){
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
      				G4cout << "### G4BetheBlochModel in Adjoint Sim WARNING: gg= " << gg
	     			<< G4endl;
				gg=1.;
			}
			//G4cout<<"gg"<<gg<<G4endl;
			dSigmadEprod*=gg;		
    	 	}
	}	
	 
  }

 return dSigmadEprod;	
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointIonIonisationModel::SetIon(G4ParticleDefinition* adj_ion, G4ParticleDefinition* fwd_ion)
{ theDirectPrimaryPartDef =fwd_ion;
  theAdjEquivOfDirectPrimPartDef =adj_ion;
 
  DefineProjectileProperty();
}
//////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointIonIonisationModel::CorrectPostStepWeight(G4ParticleChange* fParticleChange, G4double old_weight,
							G4double adjointPrimKinEnergy, G4double projectileKinEnergy,G4bool )
{
 //It is needed because the direct cross section used to compute the differential cross section is not the one used in 
 // the direct model where the GenericIon stuff is considered with correction of effective charge.  In the direct model the samnepl of secondaries does
 // not reflect the integral cross section. The integral fwd cross section that we used to compute the differential CS
 // match the sample of secondaries in the forward case despite the fact that its is not the same total CS than in the FWD case. For this reasion an extra
 // weight correction is needed at the end.
 

 G4double new_weight=old_weight;
 
 //the correction of CS due to the problem explained above
 G4double kinEnergyProjScaled = massRatio*projectileKinEnergy;
 theDirectEMModel =theBraggIonDirectEMModel;
 if (kinEnergyProjScaled >2.*MeV && !use_only_bragg) theDirectEMModel = theBetheBlochDirectEMModel; //Bethe Bloch Model
 G4double UsedFwdCS=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,projectileKinEnergy,1,1 ,currentTcutForDirectSecond,1.e20);
 G4double chargeSqRatio =1.;
 if (chargeSquare>1.) chargeSqRatio =  theDirectEMModel->GetChargeSquareRatio(theDirectPrimaryPartDef,currentMaterial,projectileKinEnergy);
 G4double  CorrectFwdCS = chargeSqRatio*theDirectEMModel->ComputeCrossSectionPerAtom(G4GenericIon::GenericIon(),kinEnergyProjScaled,1,1 ,currentTcutForDirectSecond,1.e20);
 if (UsedFwdCS >0)  new_weight*= CorrectFwdCS/UsedFwdCS;//May be some check is needed if UsedFwdCS ==0 probably that then we should avoid a secondary to be produced, 
 
 
 //additional CS crorrection  needed for cross section biasing in general. 
 //May be wrong for ions!!! Most of the time not used!
  G4double w_corr =1./CS_biasing_factor;
  w_corr*=G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
  new_weight*=w_corr;
 
  new_weight*=projectileKinEnergy/adjointPrimKinEnergy;
 
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);
}


//////////////////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointIonIonisationModel::DefineProjectileProperty()
{   
    //Slightly modified code taken from G4BetheBlochModel::SetParticle
    //------------------------------------------------
    G4String pname = theDirectPrimaryPartDef->GetParticleName();
    if (theDirectPrimaryPartDef->GetParticleType() == "nucleus" &&
	pname != "deuteron" && pname != "triton") {
      isIon = true;
    }
    
    mass = theDirectPrimaryPartDef->GetPDGMass();
    massRatio= G4GenericIon::GenericIon()->GetPDGMass()/mass;
    mass_ratio_projectile = massRatio;
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


//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy)
{ 
  G4double Tmax=PrimAdjEnergy*one_plus_ratio_2/(one_minus_ratio_2-2.*ratio*PrimAdjEnergy/mass);
  return Tmax;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMinForScatProjToProjCase(G4double PrimAdjEnergy,G4double Tcut)
{ return PrimAdjEnergy+Tcut;
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMaxForProdToProjCase(G4double )
{ return HighEnergyLimit;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointIonIonisationModel::GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy)
{  G4double Tmin= (2*PrimAdjEnergy-4*mass + std::sqrt(4.*PrimAdjEnergy*PrimAdjEnergy +16.*mass*mass + 8.*PrimAdjEnergy*mass*(1/ratio +ratio)))/4.;
   return Tmin;
}
