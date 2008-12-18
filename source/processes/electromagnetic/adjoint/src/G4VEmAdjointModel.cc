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
#include "G4VEmAdjointModel.hh"
#include "G4AdjointCSManager.hh"


#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointInterpolator.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4VEmAdjointModel::G4VEmAdjointModel(const G4String& nam):
name(nam)
// lowLimit(0.1*keV), highLimit(100.0*TeV), fluc(0), name(nam), pParticleChange(0)
{ G4AdjointCSManager::GetAdjointCSManager()->RegisterEmAdjointModel(this);
  CorrectWeightMode =true;
  UseMatrix =true;
  UseMatrixPerElement = true;
  ApplyCutInRange = true;
  ApplyBiasing = true;
  UseOnlyOneMatrixForAllElements = true;
  IsIonisation =true; 
  CS_biasing_factor =1.;
  //ApplyBiasing = false;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VEmAdjointModel::~G4VEmAdjointModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::SampleSecondaries(const G4Track& aTrack,
                       G4bool IsScatProjToProjCase,
	               G4ParticleChange* fParticleChange)
{ 

  const G4DynamicParticle* theAdjointPrimary =aTrack.GetDynamicParticle();
  //DefineCurrentMaterial(aTrack->GetMaterialCutsCouple());
  size_t ind=0;
  if (!UseMatrixPerElement) ind = currentMaterialIndex;
  //G4cout<<theAdjointPrimary<<std::endl;
  else if (!UseOnlyOneMatrixForAllElements) { //Select Material
   	std::vector<double>* CS_Vs_Element = &CS_Vs_ElementForScatProjToProjCase;
  	if ( !IsScatProjToProjCase) CS_Vs_Element = &CS_Vs_ElementForProdToProjCase;
  	G4double rand_var= G4UniformRand();
  	G4double SumCS=0.;
  	for (size_t i=0;i<CS_Vs_Element->size();i++){
 		SumCS+=(*CS_Vs_Element)[i];
		if (rand_var<=SumCS/lastCS){
			ind=i;
			break;
		}
  	}
	ind = currentMaterial->GetElement(ind)->GetIndex();
  }
  
 
 
 //Elastic inverse scattering //not correct in all the cases 
 //---------------------------------------------------------
 G4double adjointPrimKinEnergy = theAdjointPrimary->GetKineticEnergy();
 G4double adjointPrimTotalEnergy = theAdjointPrimary->GetTotalEnergy();
 G4double adjointPrimP =theAdjointPrimary->GetTotalMomentum();
 //G4cout<<adjointPrimKinEnergy<<std::endl;
 if (adjointPrimKinEnergy>HighEnergyLimit*0.999){
 	return;
 }
 
 //Sample secondary energy
 //-----------------------
 G4double projectileKinEnergy;
// if (!IsIonisation  ) {
 	projectileKinEnergy = SampleAdjSecEnergyFromCSMatrix(ind,
 						  adjointPrimKinEnergy,
						  IsScatProjToProjCase);
 //}
 /*else {
	projectileKinEnergy = 	SampleAdjSecEnergyFromDiffCrossSectionPerAtom(adjointPrimKinEnergy,IsScatProjToProjCase);
 	//G4cout<<projectileKinEnergy<<std::endl;
 }*/	
 //Weight correction
 //-----------------------
 
 CorrectPostStepWeight(fParticleChange, aTrack.GetWeight(), adjointPrimKinEnergy,projectileKinEnergy);					   
 
 		 
 
 
 
 //Kinematic
 //---------
 
 G4double projectileM0 = theAdjEquivOfDirectPrimPartDef->GetPDGMass();
 G4double projectileTotalEnergy = projectileM0+projectileKinEnergy;
 G4double projectileP2 = projectileTotalEnergy*projectileTotalEnergy - projectileM0*projectileM0;	
 
 
 
 //Companion
 //-----------
 G4double companionM0;
 companionM0=(adjointPrimTotalEnergy-adjointPrimKinEnergy);
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
 
 
 
 if (!IsScatProjToProjCase && CorrectWeightMode){ //kill the primary and add a secondary
 	fParticleChange->ProposeTrackStatus(fStopAndKill);
 	fParticleChange->AddSecondary(new G4DynamicParticle(theAdjEquivOfDirectPrimPartDef,projectileMomentum));
	//G4cout<<"projectileMomentum "<<projectileMomentum<<std::endl;
 }
 else {
 	fParticleChange->ProposeEnergy(projectileKinEnergy);
	fParticleChange->ProposeMomentumDirection(projectileMomentum.unit());
 }	
}
////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::CorrectPostStepWeight(G4ParticleChange* fParticleChange, G4double old_weight,  G4double , G4double ) 
{
 G4double new_weight=old_weight;
 if (CorrectWeightMode) {
 	G4double w_corr =1./CS_biasing_factor;
	//G4cout<<w_corr<<std::endl;
	
	/*G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection(theAdjEquivOfDirectPrimPartDef, 
									       theAdjEquivOfDirectSecondPartDef,
					    		 		       adjointPrimKinEnergy,projectileKinEnergy,
	 		   	     			  		       aTrack.GetMaterialCutsCouple());
	w_corr = projectileKinEnergy;
	G4double Emin,Emax;
	if (IsScatProjToProjCase) {
		Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(adjointPrimKinEnergy);
		Emin = GetSecondAdjEnergyMinForScatProjToProjCase(adjointPrimKinEnergy, currentTcutForDirectSecond);
		
	}
	else {
		Emax = GetSecondAdjEnergyMaxForProdToProjCase(adjointPrimKinEnergy);
		Emin = GetSecondAdjEnergyMinForProdToProjCase(adjointPrimKinEnergy);	
	}
	w_corr *=std::log(Emax/Emin)/(Emax-Emin); */	
		
	new_weight*=w_corr;
 }
 G4cout<< "new weight"<<new_weight<<std::endl;
 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);	
}
////////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double primEnergy,
				G4bool IsScatProjToProjCase)
{ 
  DefineCurrentMaterial(aCouple);
  //G4double fwdCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalForwardCS(G4AdjointElectron::AdjointElectron(),primEnergy,aCouple);
  //G4double adjCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalAdjointCS(G4AdjointElectron::AdjointElectron(), primEnergy,aCouple);
  if (IsScatProjToProjCase){
  	lastCS = G4AdjointCSManager::GetAdjointCSManager()->ComputeAdjointCS(currentMaterial,
					    		 		this, 
					    		 		primEnergy,
					    		 		currentTcutForDirectSecond,
					    		 		true,
							 		CS_Vs_ElementForScatProjToProjCase);
	/*G4double fwdCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalForwardCS(theAdjEquivOfDirectPrimPartDef,primEnergy,aCouple);
  	G4double adjCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalAdjointCS(theAdjEquivOfDirectPrimPartDef, primEnergy,aCouple);
	*/
	//if (adjCS >0 )lastCS *=fwdCS/adjCS;
 
  }
  else {
  	lastCS = G4AdjointCSManager::GetAdjointCSManager()->ComputeAdjointCS(currentMaterial,
					    		 		this, 
					    		 		primEnergy,
					    		 		currentTcutForDirectSecond,
					    		 		false,
									CS_Vs_ElementForProdToProjCase);
	/*G4double fwdCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalForwardCS(theAdjEquivOfDirectSecondPartDef,primEnergy,aCouple);
  	G4double adjCS = G4AdjointCSManager::GetAdjointCSManager()->GetTotalAdjointCS(theAdjEquivOfDirectSecondPartDef, primEnergy,aCouple);
	*/
	//if (adjCS >0 )lastCS *=fwdCS/adjCS;						 		
	//lastCS=0.;								
  }
 
  
  return lastCS;
  									
}				
////////////////////////////////////////////////////////////////////////////////
//
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine  
G4double G4VEmAdjointModel::DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyProd,
				      G4double Z, 
                                      G4double A)
{
 G4double dSigmadEprod=0;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){ //the produced particle should have a kinetic energy smaller than the projectile 
	G4double Tmax=kinEnergyProj;
	if (second_part_of_same_type) Tmax = kinEnergyProj/2.;
	return Z*DiffCrossSectionMoller(kinEnergyProj,kinEnergyProd);
	//it could be thta Tmax here should be DBLMAX
	//Tmax=DBLMAX;
	
 	G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.000001;
 	G4double dE=(E2-E1);
 	G4double sigma1=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E1,1.e20);
 	G4double sigma2=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E2,1.e20);
 	
	dSigmadEprod=(sigma1-sigma2)/dE;
	if (dSigmadEprod>1.) {
		G4cout<<"sigma1 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma1<<std::endl;
		G4cout<<"sigma2 "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<sigma2<<std::endl;
		G4cout<<"dsigma "<<kinEnergyProj/MeV<<'\t'<<kinEnergyProd/MeV<<'\t'<<dSigmadEprod<<std::endl;
		
	}
	
	
 	
 }
 
 

 	
 return dSigmadEprod;	
 
 
 
}
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine 
////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionPerAtomPrimToScatPrim(
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyScatProj,
				      G4double Z, 
                                      G4double A)
{ G4double kinEnergyProd = kinEnergyProj - kinEnergyScatProj;
  G4double dSigmadEprod;
  if (kinEnergyProd <=0) dSigmadEprod=0;
  else dSigmadEprod=DiffCrossSectionPerAtomPrimToSecond(kinEnergyProj,kinEnergyProd,Z,A);
  return dSigmadEprod;	

}

////////////////////////////////////////////////////////////////////////////////
//
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine  
G4double G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToSecond(
				      const G4Material* aMaterial,
                                      G4double kinEnergyProj, 
                                      G4double kinEnergyProd)
{
 G4double dSigmadEprod=0;
 G4double Emax_proj = GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
 G4double Emin_proj = GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
 
 
 if (kinEnergyProj>Emin_proj && kinEnergyProj<=Emax_proj){ 
	G4double Tmax=kinEnergyProj;
	if (second_part_of_same_type) Tmax = kinEnergyProj/2.;
	//it could be thta Tmax here should be DBLMAX
	//Tmax=DBLMAX;
	
 	G4double E1=kinEnergyProd;
	
 	G4double E2=kinEnergyProd*1.0001;
 	G4double dE=(E2-E1);
 	G4double sigma1=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E1,E2);
        
	//G4double sigma2=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E2,1.e50);
 	dSigmadEprod=sigma1/dE;
	if (dSigmadEprod <0) { //could happen with bremstrahlung dur to suppression effect
		G4cout<<"Halllllllllllllllllllllllllllllllllllllllllllllllo "<<kinEnergyProj<<'\t'<<E1<<'\t'<<dSigmadEprod<<std::endl;
		E1=kinEnergyProd;
		E2=E1*1.1;
		dE=E2-E1;
		sigma1=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E1,1.e50);
        	G4double sigma2=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E2,1.e50);
		dSigmadEprod=(sigma1-sigma2)/dE;
		G4cout<<dSigmadEprod<<std::endl;
	}	
	
 	
 }


 	
 return dSigmadEprod;	
 
 
 
}
//The implementation here is correct for energy loss process, for the photoelectric and compton scattering the method should be redefine 
////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionPerVolumePrimToScatPrim(
                                      const G4Material* aMaterial,
				      G4double kinEnergyProj, 
                                      G4double kinEnergyScatProj)
{ G4double kinEnergyProd = kinEnergyProj - kinEnergyScatProj;
  G4double dSigmadEprod;
  if (kinEnergyProd <=0) dSigmadEprod=0;
  else dSigmadEprod=DiffCrossSectionPerVolumePrimToSecond(aMaterial,kinEnergyProj,kinEnergyProd);
  return dSigmadEprod;	

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionFunction1(G4double kinEnergyProj){
  //return kinEnergyProj*kinEnergyProj;
  //ApplyBiasing=false;
  G4double bias_factor = CS_biasing_factor*kinEnergyProdForIntegration/kinEnergyProj;	
  if (!ApplyBiasing) bias_factor =CS_biasing_factor; 
  //G4cout<<bias_factor<<std::endl;
  if (UseMatrixPerElement ) {
  	return DiffCrossSectionPerAtomPrimToSecond(kinEnergyProj,kinEnergyProdForIntegration,ZSelectedNucleus,ASelectedNucleus)*bias_factor;
  }
  else {
  	return DiffCrossSectionPerVolumePrimToSecond(SelectedMaterial,kinEnergyProj,kinEnergyProdForIntegration)*bias_factor;
  }	
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionMoller(G4double kinEnergyProj,G4double kinEnergyProd){
  G4double electron_mass_c2=0.51099906*MeV;
  G4double energy = kinEnergyProj + electron_mass_c2;
  G4double x   = kinEnergyProd/kinEnergyProj;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  G4double beta2  = 1.0 - 1.0/gamma2;
  
  G4double g = (2.0*gam - 1.0)/gamma2;
  G4double y = 1.0 - x;
  G4double fac=twopi_mc2_rcl2/electron_mass_c2;
  G4double dCS = fac*( 1.-g + ((1.0 - g*x)/(x*x)) + ((1.0 - g*y)/(y*y)))/(beta2*(gam-1));
  return dCS/kinEnergyProj;
  
 

}  
////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionFunction2(G4double kinEnergyProj){
  //return kinEnergyProj*kinEnergyProj;	
  G4double bias_factor =  CS_biasing_factor*kinEnergyScatProjForIntegration/kinEnergyProj;
  //ApplyBiasing=false;
 if (!ApplyBiasing) bias_factor = CS_biasing_factor;
 //G4cout<<bias_factor<<std::endl; 
 if (UseMatrixPerElement ) {
  	return DiffCrossSectionPerAtomPrimToScatPrim(kinEnergyProj,kinEnergyScatProjForIntegration,ZSelectedNucleus,ASelectedNucleus)*bias_factor;
 }	
 else {
 	return DiffCrossSectionPerVolumePrimToScatPrim(SelectedMaterial,kinEnergyProj,kinEnergyScatProjForIntegration)*bias_factor;
 
 }	
 		
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerAtomForSecond(      
				G4double kinEnergyProd,
				G4double Z, 
                                G4double A ,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, G4double(G4VEmAdjointModel::*)(G4double)> integral;
  ASelectedNucleus= G4int(A);
  ZSelectedNucleus=G4int(Z);
  kinEnergyProdForIntegration = kinEnergyProd;
  
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
  G4double maxEProj= GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
  G4double E1=minEProj;
  std::vector< G4double >*  log_ESec_vector = new  std::vector< G4double >();
  std::vector< G4double >*  log_Prob_vector = new  std::vector< G4double >();
  log_ESec_vector->clear();
  log_Prob_vector->clear();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);
  
  G4double E2=std::pow(10.,G4double( G4int(std::log10(minEProj)*nbin_pro_decade)+1)/nbin_pro_decade);
  G4double fE=std::pow(10.,1./nbin_pro_decade);
  G4double int_cross_section=0.;
  
  if (std::pow(fE,5.)>(maxEProj/minEProj)) fE = std::pow(maxEProj/minEProj,0.2);
  
  while (E1 <maxEProj*0.9999999){
  	//G4cout<<E1<<'\t'<<E2<<std::endl;
	
  	int_cross_section +=integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction1,E1,std::min(E2,maxEProj*0.99999999), 10);
	//G4cout<<"int_cross_section 1 "<<'\t'<<int_cross_section<<std::endl;
	log_ESec_vector->push_back(std::log(std::min(E2,maxEProj)));
	log_Prob_vector->push_back(std::log(int_cross_section));	
	E1=E2;
	E2*=fE;
  
  }
  std::vector< std::vector<G4double>* > res_mat;
  res_mat.clear();
  if (int_cross_section >0.) {
  	res_mat.push_back(log_ESec_vector);
  	res_mat.push_back(log_Prob_vector);
  }	
  
  return res_mat;
}

/////////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerAtomForScatProj(
      				G4double kinEnergyScatProj,
				G4double Z, 
                                G4double A ,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, G4double(G4VEmAdjointModel::*)(G4double)> integral;
  ASelectedNucleus=G4int(A);
  ZSelectedNucleus=G4int(Z);
  kinEnergyScatProjForIntegration = kinEnergyScatProj;
  
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForScatProjToProjCase(kinEnergyScatProj); 
  G4double maxEProj= GetSecondAdjEnergyMaxForScatProjToProjCase(kinEnergyScatProj);
  G4double dEmax=maxEProj-kinEnergyScatProj;
  G4double dEmin=GetLowEnergyLimit();
  G4double dE1=dEmin;
  G4double dE2=dEmin;
  
  
  std::vector< G4double >*  log_ESec_vector = new std::vector< G4double >();
  std::vector< G4double >*  log_Prob_vector = new std::vector< G4double >();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins=std::max( G4int(std::log10(dEmax/dEmin))*nbin_pro_decade,5);
  G4double fE=std::pow(dEmax/dEmin,1./nbins);
  
  G4double int_cross_section=0.;
  
  while (dE1 <dEmax*0.9999999999999){
  	dE2=dE1*fE;
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction2,minEProj+dE1,std::min(minEProj+dE2,maxEProj), 20);
	//G4cout<<"int_cross_section "<<minEProj+dE1<<'\t'<<int_cross_section<<std::endl;
	log_ESec_vector->push_back(std::log(std::min(dE2,maxEProj)));
	log_Prob_vector->push_back(std::log(int_cross_section));	
	dE1=dE2;

  }
  /*G4cout<<"total int_cross_section"<<'\t'<<int_cross_section<<std::endl;
  G4cout<<"energy "<<kinEnergyScatProj<<std::endl;*/
  
  
  
  
  std::vector< std::vector<G4double> *> res_mat; 
  res_mat.clear();
  if (int_cross_section >0.) {
  	res_mat.push_back(log_ESec_vector);
  	res_mat.push_back(log_Prob_vector);
  }	
  
  return res_mat;
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForSecond(      
				G4Material* aMaterial,
				G4double kinEnergyProd,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, G4double(G4VEmAdjointModel::*)(G4double)> integral;
  SelectedMaterial= aMaterial;
  kinEnergyProdForIntegration = kinEnergyProd;
  //G4cout<<aMaterial->GetName()<<std::endl;
  //G4cout<<kinEnergyProd/MeV<<std::endl;
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
  G4double maxEProj= GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
  G4double E1=minEProj;
  std::vector< G4double >*  log_ESec_vector = new  std::vector< G4double >();
  std::vector< G4double >*  log_Prob_vector = new  std::vector< G4double >();
  log_ESec_vector->clear();
  log_Prob_vector->clear();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);
  
  G4double E2=std::pow(10.,G4double( G4int(std::log10(minEProj)*nbin_pro_decade)+1)/nbin_pro_decade);
  G4double fE=std::pow(10.,1./nbin_pro_decade);
  G4double int_cross_section=0.;
  
  if (std::pow(fE,5.)>(maxEProj/minEProj)) fE = std::pow(maxEProj/minEProj,0.2);
  
  while (E1 <maxEProj*0.9999999){
  	//G4cout<<E1<<'\t'<<E2<<std::endl;
	
  	int_cross_section +=integral.Simpson(this, &G4VEmAdjointModel::DiffCrossSectionFunction1,E1,std::min(E2,maxEProj*0.99999999), 10);
	//G4cout<<"int_cross_section 1 "<<E1<<'\t'<<int_cross_section<<std::endl;
	log_ESec_vector->push_back(std::log(std::min(E2,maxEProj)));
	log_Prob_vector->push_back(std::log(int_cross_section));
	E1=E2;
	E2*=fE;
  
  }
  std::vector< std::vector<G4double>* > res_mat;
  res_mat.clear();
  
  //if (int_cross_section >0.) {
  	res_mat.push_back(log_ESec_vector);
  	res_mat.push_back(log_Prob_vector);
  //}	
  
  return res_mat;
}

/////////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForScatProj(
      				G4Material* aMaterial,
				G4double kinEnergyScatProj,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, G4double(G4VEmAdjointModel::*)(G4double)> integral;
  SelectedMaterial= aMaterial;
  kinEnergyScatProjForIntegration = kinEnergyScatProj;
  /*G4cout<<name<<std::endl;
  G4cout<<aMaterial->GetName()<<std::endl;
  G4cout<<kinEnergyScatProj/MeV<<std::endl;*/
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForScatProjToProjCase(kinEnergyScatProj);
  G4double maxEProj= GetSecondAdjEnergyMaxForScatProjToProjCase(kinEnergyScatProj);
 
  
  G4double dEmax=maxEProj-kinEnergyScatProj;
  G4double dEmin=GetLowEnergyLimit();
  G4double dE1=dEmin;
  G4double dE2=dEmin;
  
  
  std::vector< G4double >*  log_ESec_vector = new std::vector< G4double >();
  std::vector< G4double >*  log_Prob_vector = new std::vector< G4double >();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins=std::max( G4int(std::log10(dEmax/dEmin))*nbin_pro_decade,5);
  G4double fE=std::pow(dEmax/dEmin,1./nbins);
  
  G4double int_cross_section=0.;
  
  while (dE1 <dEmax*0.9999999999999){
  	dE2=dE1*fE;
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction2,minEProj+dE1,std::min(minEProj+dE2,maxEProj), 20);
	//G4cout<<"int_cross_section "<<minEProj+dE1<<'\t'<<int_cross_section<<std::endl;
	log_ESec_vector->push_back(std::log(std::min(dE2,maxEProj)));
	log_Prob_vector->push_back(std::log(int_cross_section));
	dE1=dE2;

  }
  
  
  
  
  
  std::vector< std::vector<G4double> *> res_mat;
  res_mat.clear();
  if (int_cross_section >0.) {
  	res_mat.push_back(log_ESec_vector);
  	res_mat.push_back(log_Prob_vector);
  }	
  
  return res_mat;
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromCSMatrix(size_t MatrixIndex,G4double aPrimEnergy,G4bool IsScatProjToProjCase)
{ 
  
  
  G4AdjointCSMatrix* theMatrix= (*pOnCSMatrixForProdToProjBackwardScattering)[MatrixIndex];
  if (IsScatProjToProjCase) theMatrix= (*pOnCSMatrixForScatProjToProjBackwardScattering)[MatrixIndex];
  std::vector< G4double >* theLogPrimEnergyVector = theMatrix->GetLogPrimEnergyVector();
  //G4double dLog = theMatrix->GetDlog();
  
  
  
  if (theLogPrimEnergyVector->size() ==0){
 	G4cout<<"No data are contained in the given AdjointCSMatrix!"<<std::endl;
	G4cout<<"The sampling procedure will be stopped."<<std::endl;
	return 0.;
	
  }
  
  G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
  G4double aLogPrimEnergy = std::log(aPrimEnergy);
  size_t ind =theInterpolator->FindPositionForLogVector(aLogPrimEnergy,*theLogPrimEnergyVector);
  
  
  G4double aLogPrimEnergy1,aLogPrimEnergy2;
  G4double aLogCS1,aLogCS2;
   G4double log01,log02;
  std::vector< G4double>* aLogSecondEnergyVector1 =0;
  std::vector< G4double>* aLogSecondEnergyVector2  =0;
  std::vector< G4double>* aLogProbVector1=0;
  std::vector< G4double>* aLogProbVector2=0; 
  std::vector< size_t>* aLogProbVectorIndex1=0;
  std::vector< size_t>* aLogProbVectorIndex2=0;
	 							     
  theMatrix->GetData(ind, aLogPrimEnergy1,aLogCS1,log01, aLogSecondEnergyVector1,aLogProbVector1,aLogProbVectorIndex1);
  theMatrix->GetData(ind+1, aLogPrimEnergy2,aLogCS2,log02, aLogSecondEnergyVector2,aLogProbVector2,aLogProbVectorIndex2);
  
  G4double rand_var = G4UniformRand();
  G4double log_rand_var= std::log(rand_var);
  G4double log_Tcut =std::log(currentTcutForDirectSecond);
  G4double Esec=0;
  G4double log_dE1,log_dE2;
  G4double log_rand_var1,log_rand_var2;
  G4double log_E1,log_E2;
  log_rand_var1=log_rand_var;
  log_rand_var2=log_rand_var;
  
  G4double Emin=0.;
  G4double Emax=0.;
  if (theMatrix->IsScatProjToProjCase()){ //case where Tcut plays a role
 	//G4cout<<"Here "<<std::endl;
	if (ApplyCutInRange) {
		if (second_part_of_same_type && currentTcutForDirectSecond>aPrimEnergy) return aPrimEnergy;
		/*if (IsIonisation){
			G4double inv_Tcut= 1./currentTcutForDirectSecond;
			G4double inv_dE=inv_Tcut-rand_var*(inv_Tcut-1./aPrimEnergy);
			Esec= aPrimEnergy+1./inv_dE;
			//return Esec;
			G4double dE1=currentTcutForDirectSecond;
			G4double dE2=currentTcutForDirectSecond*1.00001;
			G4double dCS1=DiffCrossSectionMoller(aPrimEnergy+dE1,dE1);
			G4double dCS2=DiffCrossSectionMoller(aPrimEnergy+dE2,dE2);
			G4double alpha1=std::log(dCS1/dCS2)/std::log(dE1/dE2);
			G4double a1=dCS1/std::pow(dE1,alpha1);
			dCS1=DiffCrossSectionMoller(aPrimEnergy+dE1,dE1);
			dCS2=DiffCrossSectionMoller(aPrimEnergy+dE2,dE2);
			
			return Esec;
			
			
			
			dE1=aPrimEnergy/1.00001;
			dE2=aPrimEnergy;
			dCS1=DiffCrossSectionMoller(aPrimEnergy+dE1,dE1);
			dCS2=DiffCrossSectionMoller(aPrimEnergy+dE2,dE2);
			G4double alpha2=std::log(dCS1/dCS2)/std::log(dE1/dE2);
			G4double a2=dCS1/std::pow(dE1,alpha1);
			return Esec;
			
			
			
			
		}*/	
		log_rand_var1=log_rand_var+theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector1,*aLogProbVector1);
		log_rand_var2=log_rand_var+theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector2,*aLogProbVector2);
		
	}	
	log_dE1 = theInterpolator->Interpolate(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,"Lin");
 	log_dE2 = theInterpolator->Interpolate(log_rand_var2,*aLogProbVector2,*aLogSecondEnergyVector2,"Lin");
	
	/*log_dE1 = theInterpolator->InterpolateWithIndexVector(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,*aLogProbVectorIndex1,log01,dLog);
	log_dE2 = theInterpolator->InterpolateWithIndexVector(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,*aLogProbVectorIndex1,log02,dLog);
	*/				    
					    
	
	
	
	
	Esec = aPrimEnergy +
	std::exp(theInterpolator->LinearInterpolation(aLogPrimEnergy,aLogPrimEnergy1,aLogPrimEnergy2,log_dE1,log_dE2));
	
	Emin=GetSecondAdjEnergyMinForScatProjToProjCase(aPrimEnergy);
	Emax=GetSecondAdjEnergyMaxForScatProjToProjCase(aPrimEnergy);
	Esec=std::max(Esec,Emin);
	Esec=std::min(Esec,Emax);
	
	
	//G4cout<<"Esec "<<Esec<<std::endl;
	//if (Esec > 2.*aPrimEnergy && second_part_of_same_type) Esec = 2.*aPrimEnergy; 
  }
  else { //Tcut condition is already full-filled
        /*G4cout<<"Start "<<std::endl;
  	G4cout<<std::exp((*aLogProbVector1)[0])<<std::endl;
	G4cout<<std::exp((*aLogProbVector2)[0])<<std::endl;*/
	/*G4double inv_E1= .5/aPrimEnergy;
		
	G4double inv_E=inv_E1-rand_var*(inv_E1-0.00001);
	Esec= 1./inv_E;
	return Esec;*/
 	log_E1 = theInterpolator->Interpolate(log_rand_var,*aLogProbVector1,*aLogSecondEnergyVector1,"Lin");
 	log_E2 = theInterpolator->Interpolate(log_rand_var,*aLogProbVector2,*aLogSecondEnergyVector2,"Lin");
	/*log_E1 = theInterpolator->InterpolateWithIndexVector(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,*aLogProbVectorIndex1,log01,dLog);
	log_E2 = theInterpolator->InterpolateWithIndexVector(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,*aLogProbVectorIndex1,log02,dLog);
	*/
	
	
	/*G4cout<<std::exp(log_E1)<<std::endl;
	G4cout<<std::exp(log_E2)<<std::endl;*/
	
	Esec = std::exp(theInterpolator->LinearInterpolation(aLogPrimEnergy,aLogPrimEnergy1,aLogPrimEnergy2,log_E1,log_E2));
	Emin=GetSecondAdjEnergyMinForProdToProjCase(aPrimEnergy);
	Emax=GetSecondAdjEnergyMaxForProdToProjCase(aPrimEnergy);
	Esec=std::max(Esec,Emin);
	Esec=std::min(Esec,Emax);
	
  }
  	
  return Esec;
  
  
  
 
										   
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromDiffCrossSectionPerAtom(G4double prim_energy,G4bool IsScatProjToProjCase)
{  
  // here we try to use the rejection method 
  //-----------------------------------------
  
  G4double E=0;
  G4double x,xmin,greject,q;
  if ( IsScatProjToProjCase){
  	G4double Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(prim_energy);
        G4double Emin= prim_energy+currentTcutForDirectSecond;
	xmin=Emin/Emax;
	G4double grejmax = DiffCrossSectionPerAtomPrimToScatPrim(Emin,prim_energy,1)*prim_energy; 

	do {
        	q = G4UniformRand();
        	x = 1./(q*(1./xmin -1.) +1.);
		E=x*Emax;
		greject = DiffCrossSectionPerAtomPrimToScatPrim( E,prim_energy ,1)*prim_energy; 
		
	}	
        
	while( greject < G4UniformRand()*grejmax );	
	
  }
  else {
  	G4double Emax = GetSecondAdjEnergyMaxForProdToProjCase(prim_energy);
        G4double Emin=  GetSecondAdjEnergyMinForProdToProjCase(prim_energy);;
	xmin=Emin/Emax;
	G4double grejmax = DiffCrossSectionPerAtomPrimToSecond(Emin,prim_energy,1); 
	do {
        	q = G4UniformRand();
        	x = std::pow(xmin, q);
		E=x*Emax;
		greject = DiffCrossSectionPerAtomPrimToSecond( E,prim_energy ,1); 
		
	}	
        
	while( greject < G4UniformRand()*grejmax );
  
  
  
  }
  
  return E;
  
  
  
 
										   
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::GetSecondAdjEnergyMaxForScatProjToProjCase(G4double kinEnergyScatProj)
{ G4double maxEProj= HighEnergyLimit;
  if (second_part_of_same_type)  maxEProj=std::min(kinEnergyScatProj*2.,HighEnergyLimit);
  return maxEProj;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::GetSecondAdjEnergyMinForScatProjToProjCase(G4double PrimAdjEnergy,G4double Tcut)
{ return PrimAdjEnergy+Tcut;
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::GetSecondAdjEnergyMaxForProdToProjCase(G4double )
{ return HighEnergyLimit;
}
//////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy)
{ G4double minEProj=PrimAdjEnergy;
  if (second_part_of_same_type)  minEProj=PrimAdjEnergy*2.;
  return minEProj;
}
////////////////////////////////////////////////////////////////////////////////////////////
//
void  G4VEmAdjointModel::DefineCurrentMaterial(const G4MaterialCutsCouple* couple)
{ if(couple != currentCouple) {
    	currentCouple   = const_cast<G4MaterialCutsCouple*> (couple);
    	currentMaterial = const_cast<G4Material*> (couple->GetMaterial());
    	currentCoupleIndex = couple->GetIndex();
    	currentMaterialIndex = currentMaterial->GetIndex();
   	size_t idx=56;
    
   	if (theAdjEquivOfDirectPrimPartDef) {
    		if (theAdjEquivOfDirectPrimPartDef->GetParticleName() == "adj_gamma") idx = 0;
   		else if (theAdjEquivOfDirectPrimPartDef->GetParticleName() == "adj_e-") idx = 1;
    		else if (theAdjEquivOfDirectPrimPartDef->GetParticleName() == "adj_e+") idx = 2;
    		const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
    		currentTcutForDirectPrim=(*aVec)[currentCoupleIndex];
    	}	
    
    
    	if (theAdjEquivOfDirectPrimPartDef == theAdjEquivOfDirectSecondPartDef) {
    		currentTcutForDirectSecond = currentTcutForDirectPrim;
    	}
   	else { 
    		if (theAdjEquivOfDirectSecondPartDef){
			if (theAdjEquivOfDirectSecondPartDef->GetParticleName() == "adj_gamma") idx = 0;
   			else if (theAdjEquivOfDirectSecondPartDef->GetParticleName() == "adj_e-") idx = 1;
    			else if (theAdjEquivOfDirectSecondPartDef->GetParticleName() == "adj_e+") idx = 2;
    			const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
    			currentTcutForDirectSecond=(*aVec)[currentCoupleIndex];
		}
    	}
  }
}
