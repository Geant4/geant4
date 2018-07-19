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
// $Id: G4AdjointBremsstrahlungModel.cc 100666 2016-10-31 10:27:00Z gcosmo $
//
#include "G4AdjointBremsstrahlungModel.hh"
#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4Electron.hh"
#include "G4Timer.hh"
#include "G4SeltzerBergerModel.hh"


////////////////////////////////////////////////////////////////////////////////
//
G4AdjointBremsstrahlungModel::G4AdjointBremsstrahlungModel(G4VEmModel* aModel):
 G4VEmAdjointModel("AdjointeBremModel")
{ 
  SetUseMatrix(false);
  SetUseMatrixPerElement(false);
  
  theDirectStdBremModel = aModel;
  theDirectEMModel=theDirectStdBremModel;
  theEmModelManagerForFwdModels = new G4EmModelManager();
  isDirectModelInitialised = false;
  G4VEmFluctuationModel* f=0;
  G4Region* r=0;
  theEmModelManagerForFwdModels->AddEmModel(1, theDirectStdBremModel, f, r);

  SetApplyCutInRange(true);
  highKinEnergy= 1.*GeV;
  lowKinEnergy = 1.0*keV;

  lastCZ =0.;

  
  theAdjEquivOfDirectPrimPartDef =G4AdjointElectron::AdjointElectron();
  theAdjEquivOfDirectSecondPartDef=G4AdjointGamma::AdjointGamma();
  theDirectPrimaryPartDef=G4Electron::Electron();
  second_part_of_same_type=false;

  
  CS_biasing_factor =1.;

  
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointBremsstrahlungModel::G4AdjointBremsstrahlungModel():
 G4VEmAdjointModel("AdjointeBremModel")
{
  SetUseMatrix(false);
  SetUseMatrixPerElement(false);

  theDirectStdBremModel = new G4SeltzerBergerModel();
  theDirectEMModel=theDirectStdBremModel;
  theEmModelManagerForFwdModels = new G4EmModelManager();
  isDirectModelInitialised = false;
  G4VEmFluctuationModel* f=0;
  G4Region* r=0;
  theEmModelManagerForFwdModels->AddEmModel(1, theDirectStdBremModel, f, r);
 // theDirectPenelopeBremModel =0;
  SetApplyCutInRange(true);
  highKinEnergy= 1.*GeV;
  lowKinEnergy = 1.0*keV;
  lastCZ =0.;
  theAdjEquivOfDirectPrimPartDef =G4AdjointElectron::AdjointElectron();
  theAdjEquivOfDirectSecondPartDef=G4AdjointGamma::AdjointGamma();
  theDirectPrimaryPartDef=G4Electron::Electron();
  second_part_of_same_type=false;
}
////////////////////////////////////////////////////////////////////////////////
//
G4AdjointBremsstrahlungModel::~G4AdjointBremsstrahlungModel()
{if (theDirectStdBremModel) delete theDirectStdBremModel;
 if (theEmModelManagerForFwdModels) delete theEmModelManagerForFwdModels;
}

////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointBremsstrahlungModel::SampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{
 if (!UseMatrix) return RapidSampleSecondaries(aTrack,IsScatProjToProjCase,fParticleChange); 

 const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
 DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());
 
 
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();
 
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
  
  G4double projectileKinEnergy = SampleAdjSecEnergyFromCSMatrix(adjointPrimKinEnergy,
						  	IsScatProjToProjCase);
 //Weight correction
 //-----------------------					   
 CorrectPostStepWeight(fParticleChange, 
 		       aTrack.GetWeight(), 
		       adjointPrimKinEnergy,
		       projectileKinEnergy,
		       IsScatProjToProjCase);	
 
 
 //Kinematic
 //---------
 G4double projectileM0 = theAdjEquivOfDirectPrimPartDef->GetPDGMass();
 G4double projectileTotalEnergy = projectileM0+projectileKinEnergy;
 G4double projectileP2 = projectileTotalEnergy*projectileTotalEnergy - projectileM0*projectileM0;	
 G4double projectileP = std::sqrt(projectileP2);
 
 
 //Angle of the gamma direction with the projectile taken from G4eBremsstrahlungModel
 //------------------------------------------------
  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - std::log(G4UniformRand()*G4UniformRand())/a1;
     else                          u = - std::log(G4UniformRand()*G4UniformRand())/a2;

  G4double theta = u*electron_mass_c2/projectileTotalEnergy;

  G4double sint = std::sin(theta);
  G4double cost = std::cos(theta);

  G4double phi = twopi * G4UniformRand() ;
  
  G4ThreeVector projectileMomentum;
  projectileMomentum=G4ThreeVector(std::cos(phi)*sint,std::sin(phi)*sint,cost)*projectileP; //gamma frame
  if (IsScatProjToProjCase) {//the adjoint primary is the scattered e-
  	G4ThreeVector gammaMomentum = (projectileTotalEnergy-adjointPrimTotalEnergy)*G4ThreeVector(0.,0.,1.);
	G4ThreeVector dirProd=projectileMomentum-gammaMomentum;
	G4double cost1 = std::cos(dirProd.angle(projectileMomentum));
	G4double sint1 =  std::sqrt(1.-cost1*cost1);
	projectileMomentum=G4ThreeVector(std::cos(phi)*sint1,std::sin(phi)*sint1,cost1)*projectileP;
  
  }
  
  projectileMomentum.rotateUz(theAdjointPrimary->GetMomentumDirection());
 
 
 
  if (!IsScatProjToProjCase ){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,projectileMomentum));
  }
  else {
 	fParticleChange->ProposeEnergy(projectileKinEnergy);
	fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
	
  }	
} 
////////////////////////////////////////////////////////////////////////////////
//
void G4AdjointBremsstrahlungModel::RapidSampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 

 const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
 DefineCurrentMaterial(aTrack.GetMaterialCutsCouple());
 
 
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();
 
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
  
 G4double projectileKinEnergy =0.;
 G4double gammaEnergy=0.;
 G4double diffCSUsed=0.; 
 if (!IsScatProjToProjCase){
 	gammaEnergy=adjointPrimKinEnergy;
	G4double Emax = GetSecondAdjEnergyMaxForProdToProjCase(adjointPrimKinEnergy);
        G4double Emin=  GetSecondAdjEnergyMinForProdToProjCase(adjointPrimKinEnergy);;
	if (Emin>=Emax) return;
	projectileKinEnergy=Emin*std::pow(Emax/Emin,G4UniformRand());
	diffCSUsed=CS_biasing_factor*lastCZ/projectileKinEnergy;
 	
 }
 else {	G4double Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(adjointPrimKinEnergy);
	G4double Emin = GetSecondAdjEnergyMinForScatProjToProjCase(adjointPrimKinEnergy,currentTcutForDirectSecond);
	if (Emin>=Emax) return;
	G4double f1=(Emin-adjointPrimKinEnergy)/Emin;
	G4double f2=(Emax-adjointPrimKinEnergy)/Emax/f1;
	projectileKinEnergy=adjointPrimKinEnergy/(1.-f1*std::pow(f2,G4UniformRand()));
	gammaEnergy=projectileKinEnergy-adjointPrimKinEnergy;
	diffCSUsed=lastCZ*adjointPrimKinEnergy/projectileKinEnergy/gammaEnergy;
	
 }
  
  
  
						  	
 //Weight correction
 //-----------------------
 //First w_corr is set to the ratio between adjoint total CS and fwd total CS
 //if this has to be done in the model
 //For the case of forced interaction this will be done in the PostStepDoIt of the
 //forced interaction
 //It is important to set the weight before the vreation of the  secondary
 //
 G4double w_corr=additional_weight_correction_factor_for_post_step_outside_model;
 if (correct_weight_for_post_step_in_model) {
     w_corr=G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
 }
 //G4cout<<"Correction factor start in brem model "<<w_corr<<std::endl;


 //Then another correction is needed due to the fact that a biaised differential CS has been used rather than the one consistent with the direct model
 //Here we consider the true  diffCS as the one obtained by the numericla differentiation over Tcut of the direct CS, corrected by the Migdal term.
 //Basically any other differential CS   diffCS could be used here (example Penelope). 

 G4double diffCS = DiffCrossSectionPerVolumePrimToSecond(currentMaterial, projectileKinEnergy, gammaEnergy);
 /*G4cout<<"diffCS "<<diffCS <<std::endl;
 G4cout<<"diffCS_Used "<<diffCSUsed <<std::endl;*/
 w_corr*=diffCS/diffCSUsed;


 G4double new_weight = aTrack.GetWeight()*w_corr;
 /*G4cout<<"New weight brem "<<new_weight<<std::endl;
 G4cout<<"Weight correction brem "<<w_corr<<std::endl;*/
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);
 
 //Kinematic
 //---------
 G4double projectileM0 = theAdjEquivOfDirectPrimPartDef->GetPDGMass();
 G4double projectileTotalEnergy = projectileM0+projectileKinEnergy;
 G4double projectileP2 = projectileTotalEnergy*projectileTotalEnergy - projectileM0*projectileM0;	
 G4double projectileP = std::sqrt(projectileP2);
 

 //Use the angular model of the forward model to generate the gamma direction
 //---------------------------------------------------------------------------
//Dum dynamic particle to use the model
 G4DynamicParticle * aDynPart = new G4DynamicParticle(G4Electron::Electron(),G4ThreeVector(0.,0.,1.)*projectileP);

 //Get the element from the direct model
 const G4Element* elm = theDirectEMModel->SelectRandomAtom(currentCouple,G4Electron::Electron(),
                                                        projectileKinEnergy,currentTcutForDirectSecond);
 G4int Z=elm->GetZasInt();
 G4double energy = aDynPart->GetTotalEnergy()-gammaEnergy;
 G4ThreeVector projectileMomentum =
   theDirectEMModel->GetAngularDistribution()->SampleDirection(aDynPart,energy,Z,currentMaterial)*projectileP;
  G4double phi = projectileMomentum.getPhi();

/*
 //Angle of the gamma direction with the projectile taken from G4eBremsstrahlungModel
 //------------------------------------------------
  G4double u;
  const G4double a1 = 0.625 , a2 = 3.*a1 , d = 27. ;

  if (9./(9.+d) > G4UniformRand()) u = - std::log(G4UniformRand()*G4UniformRand())/a1;
     else                          u = - std::log(G4UniformRand()*G4UniformRand())/a2;

  G4double theta = u*electron_mass_c2/projectileTotalEnergy;

  G4double sint = std::sin(theta);
  G4double cost = std::cos(theta);

  G4double phi = twopi * G4UniformRand() ;
  G4ThreeVector projectileMomentum;
  projectileMomentum=G4ThreeVector(std::cos(phi)*sint,std::sin(phi)*sint,cost)*projectileP; //gamma frame
*/
  if (IsScatProjToProjCase) {//the adjoint primary is the scattered e-
  	G4ThreeVector gammaMomentum = (projectileTotalEnergy-adjointPrimTotalEnergy)*G4ThreeVector(0.,0.,1.);
	G4ThreeVector dirProd=projectileMomentum-gammaMomentum;
	G4double cost1 = std::cos(dirProd.angle(projectileMomentum));
	G4double sint1 =  std::sqrt(1.-cost1*cost1);
	projectileMomentum=G4ThreeVector(std::cos(phi)*sint1,std::sin(phi)*sint1,cost1)*projectileP;
  }

  projectileMomentum.rotateUz(theAdjointPrimary->GetMomentumDirection());
 
  if (!IsScatProjToProjCase ){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,projectileMomentum));
  }
  else {
 	fParticleChange->ProposeEnergy(projectileKinEnergy);
	fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
  }	
} 
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecond(const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{if (!isDirectModelInitialised) {
	theEmModelManagerForFwdModels->Initialise(G4Electron::Electron(),G4Gamma::Gamma(),1.,0);
	isDirectModelInitialised =true;
 }
/*
 return  DiffCrossSectionPerVolumePrimToSecondApproximated2(aMaterial,
                                      			         kinEnergyProj, 
                                      			         kinEnergyProd);
*/
return G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToSecond(aMaterial,
                                      			         kinEnergyProj, 
                                      			         kinEnergyProd);								 								 
}				      

////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecondApproximated1(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 G4double dCrossEprod=0.;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 //In this approximation we consider that the secondary gammas are sampled with 1/Egamma energy distribution
 //This is what is applied in the discrete standard model before the  rejection test  that make a correction
 //The application of the same rejection function is not possible here.
 //The differentiation of the CS over Ecut does not produce neither a good differential CS. That is due to the 
 // fact that in the discrete model the differential CS and the integrated CS are both fitted but separatly and 
 // therefore do not allow a correct numerical differentiation of the integrated CS to get the differential one. 
 // In the future we plan to use the brem secondary spectra from the G4Penelope implementation 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){
 	G4double sigma=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,1.*keV);
	dCrossEprod=sigma/kinEnergyProd/std::log(kinEnergyProj/keV);
 }
 return dCrossEprod;
  
}

////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::DiffCrossSectionPerVolumePrimToSecondApproximated2(
  				      const G4Material* material,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      )
{
 //In this approximation we derive the direct cross section over Tcut=gamma energy, en after apply the Migdla correction factor 
  //used in the direct model
 
 G4double dCrossEprod=0.;
 
 const G4ElementVector* theElementVector = material->GetElementVector();
 const double* theAtomNumDensityVector = material->GetAtomicNumDensityVector();
 G4double dum=0.;
 G4double E1=kinEnergyProd,E2=kinEnergyProd*1.001;
 G4double dE=E2-E1;
 for (size_t i=0; i<material->GetNumberOfElements(); i++) { 
 	G4double C1=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,(*theElementVector)[i]->GetZ(),dum ,E1);
	G4double C2=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,(*theElementVector)[i]->GetZ(),dum,E2);
	dCrossEprod += theAtomNumDensityVector[i] * (C1-C2)/dE;
 }
 
 return dCrossEprod;
  
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointBremsstrahlungModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{ if (!isDirectModelInitialised) {
  	theEmModelManagerForFwdModels->Initialise(G4Electron::Electron(),G4Gamma::Gamma(),1.,0);
	isDirectModelInitialised =true;
  }
  if (UseMatrix) return G4VEmAdjointModel::AdjointCrossSection(aCouple,primEnergy,IsScatProjToProjCase);
  DefineCurrentMaterial(aCouple);
  G4double Cross=0.;
  lastCZ=theDirectEMModel->CrossSectionPerVolume(aCouple->GetMaterial(),theDirectPrimaryPartDef,100.*MeV,100.*MeV/std::exp(1.));//this give the constant above
  
  if (!IsScatProjToProjCase ){
  	G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(primEnergy);
  	G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(primEnergy);
	if (Emax_proj>Emin_proj && primEnergy > currentTcutForDirectSecond) Cross= CS_biasing_factor*lastCZ*std::log(Emax_proj/Emin_proj);
  }
  else {
  	G4double Emax_proj = GetSecondAdjEnergyMaxForScatProjToProjCase(primEnergy);
	G4double Emin_proj = GetSecondAdjEnergyMinForScatProjToProjCase(primEnergy,currentTcutForDirectSecond);
	if (Emax_proj>Emin_proj) Cross= lastCZ*std::log((Emax_proj-primEnergy)*Emin_proj/Emax_proj/(Emin_proj-primEnergy));
  	
  }
  return Cross;	
}					     

G4double G4AdjointBremsstrahlungModel::GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{ 
  return AdjointCrossSection(aCouple, primEnergy,IsScatProjToProjCase);
  lastCZ=theDirectEMModel->CrossSectionPerVolume(aCouple->GetMaterial(),theDirectPrimaryPartDef,100.*MeV,100.*MeV/std::exp(1.));//this give the constant above
  return G4VEmAdjointModel::GetAdjointCrossSection(aCouple, primEnergy,IsScatProjToProjCase);
  	
}



