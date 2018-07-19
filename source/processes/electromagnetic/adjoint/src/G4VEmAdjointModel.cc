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
// $Id: G4VEmAdjointModel.cc 100341 2016-10-18 08:02:25Z gcosmo $
//
#include "G4VEmAdjointModel.hh"
#include "G4AdjointCSManager.hh"
#include "G4Integrator.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleChange.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4AdjointPositron.hh"
#include "G4AdjointInterpolator.hh"
#include "G4PhysicsTable.hh"

////////////////////////////////////////////////////////////////////////////////
//
G4VEmAdjointModel::G4VEmAdjointModel(const G4String& nam):
name(nam)
// lowLimit(0.1*keV), highLimit(100.0*TeV), fluc(0), name(nam), pParticleChange(0)
{ 
  model_index = G4AdjointCSManager::GetAdjointCSManager()->RegisterEmAdjointModel(this);
  second_part_of_same_type =false;
  theDirectEMModel=0;
  mass_ratio_product=1.;
  mass_ratio_projectile=1.;
  currentCouple=0;
  additional_weight_correction_factor_for_post_step_outside_model=1.;
}
////////////////////////////////////////////////////////////////////////////////
//
G4VEmAdjointModel::~G4VEmAdjointModel()
{;}
////////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double primEnergy,
				G4bool IsScatProjToProjCase)
{ 
  DefineCurrentMaterial(aCouple);
  preStepEnergy=primEnergy;
  
  std::vector<G4double>* CS_Vs_Element = &CS_Vs_ElementForProdToProjCase;
  if (IsScatProjToProjCase)   CS_Vs_Element = &CS_Vs_ElementForScatProjToProjCase;
  lastCS = G4AdjointCSManager::GetAdjointCSManager()->ComputeAdjointCS(currentMaterial,
					    		 		this, 
					    		 		primEnergy,
					    		 		currentTcutForDirectSecond,
					    		 		IsScatProjToProjCase,
							 		*CS_Vs_Element);
  if (IsScatProjToProjCase) lastAdjointCSForScatProjToProjCase = lastCS;
  else lastAdjointCSForProdToProjCase =lastCS;								
	
 
 
  return lastCS;
  									
}
////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase)
{
  return AdjointCrossSection(aCouple, primEnergy,
				IsScatProjToProjCase);
  
  /*
  //To continue
  DefineCurrentMaterial(aCouple);
  preStepEnergy=primEnergy;
  if (IsScatProjToProjCase){
  	G4double ekin=primEnergy*mass_ratio_projectile;
  	lastCS = G4AdjointCSManager::GetAdjointCSManager()->GetAdjointSigma(ekin, model_index,true, aCouple);
	lastAdjointCSForScatProjToProjCase = lastCS;
	//G4cout<<ekin<<std::endl;
  }
  else {
  	G4double ekin=primEnergy*mass_ratio_product;
	lastCS = G4AdjointCSManager::GetAdjointCSManager()->GetAdjointSigma(ekin, model_index,false, aCouple);
	lastAdjointCSForProdToProjCase = lastCS;
	//G4cout<<ekin<<std::endl;
  }
  return lastCS;
  */
}					     				
////////////////////////////////////////////////////////////////////////////////
//
//General implementation correct for energy loss process, for the photoelectric and compton scattering the method should be redefine  
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
	
	/*G4double Tmax=kinEnergyProj;
	if (second_part_of_same_type) Tmax = kinEnergyProj/2.;*/

	G4double E1=kinEnergyProd;
 	G4double E2=kinEnergyProd*1.000001;
 	G4double dE=(E2-E1);
 	G4double sigma1=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E1,1.e20);
 	G4double sigma2=theDirectEMModel->ComputeCrossSectionPerAtom(theDirectPrimaryPartDef,kinEnergyProj,Z,A ,E2,1.e20);
 	
	dSigmadEprod=(sigma1-sigma2)/dE;
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
	/*G4double Tmax=kinEnergyProj;
	if (second_part_of_same_type) Tmax = kinEnergyProj/2.;*/
	G4double E1=kinEnergyProd;
	G4double E2=kinEnergyProd*1.0001;
 	G4double dE=(E2-E1);
 	G4double sigma1=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E1,1.e20);
    G4double sigma2=theDirectEMModel->CrossSectionPerVolume(aMaterial,theDirectPrimaryPartDef,kinEnergyProj,E2,1.e20);
 	dSigmadEprod=(sigma1-sigma2)/dE;
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
  
  
  G4double bias_factor = CS_biasing_factor*kinEnergyProdForIntegration/kinEnergyProj;	


  if (UseMatrixPerElement ) {
  	return DiffCrossSectionPerAtomPrimToSecond(kinEnergyProj,kinEnergyProdForIntegration,ZSelectedNucleus,ASelectedNucleus)*bias_factor;
  }
  else  {
  	return DiffCrossSectionPerVolumePrimToSecond(SelectedMaterial,kinEnergyProj,kinEnergyProdForIntegration)*bias_factor;
  }	
}

////////////////////////////////////////////////////////////////////////////////
//
G4double G4VEmAdjointModel::DiffCrossSectionFunction2(G4double kinEnergyProj){
  
 G4double bias_factor =  CS_biasing_factor*kinEnergyScatProjForIntegration/kinEnergyProj;
 if (UseMatrixPerElement ) {
  	return DiffCrossSectionPerAtomPrimToScatPrim(kinEnergyProj,kinEnergyScatProjForIntegration,ZSelectedNucleus,ASelectedNucleus)*bias_factor;
 }	
 else {  
 	return DiffCrossSectionPerVolumePrimToScatPrim(SelectedMaterial,kinEnergyProj,kinEnergyScatProjForIntegration)*bias_factor;
 
 }	
 		
}
////////////////////////////////////////////////////////////////////////////////
//

G4double G4VEmAdjointModel::DiffCrossSectionPerVolumeFunctionForIntegrationOverEkinProj(G4double kinEnergyProd)
{
  return DiffCrossSectionPerVolumePrimToSecond(SelectedMaterial,kinEnergyProjForIntegration,kinEnergyProd);
}
////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerAtomForSecond(      
				G4double kinEnergyProd,
				G4double Z, 
                                G4double A ,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ 
  G4Integrator<G4VEmAdjointModel, double(G4VEmAdjointModel::*)(double)> integral;
  ASelectedNucleus= int(A);
  ZSelectedNucleus=int(Z);
  kinEnergyProdForIntegration = kinEnergyProd;
  
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
  G4double maxEProj= GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
  G4double E1=minEProj;
  std::vector< double>*  log_ESec_vector = new  std::vector< double>();
  std::vector< double>*  log_Prob_vector = new  std::vector< double>();
  log_ESec_vector->clear();
  log_Prob_vector->clear();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);
  
  G4double E2=std::pow(10.,double( int(std::log10(minEProj)*nbin_pro_decade)+1)/nbin_pro_decade);
  G4double fE=std::pow(10.,1./nbin_pro_decade);
  G4double int_cross_section=0.;
  
  if (std::pow(fE,5.)>(maxEProj/minEProj)) fE = std::pow(maxEProj/minEProj,0.2);
  
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while (E1 <maxEProj*0.9999999){
  	//G4cout<<E1<<'\t'<<E2<<G4endl;
	
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction1,E1,std::min(E2,maxEProj*0.99999999), 5);
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
{ G4Integrator<G4VEmAdjointModel, double(G4VEmAdjointModel::*)(double)> integral;
  ASelectedNucleus=int(A);
  ZSelectedNucleus=int(Z);
  kinEnergyScatProjForIntegration = kinEnergyScatProj;
  
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForScatProjToProjCase(kinEnergyScatProj); 
  G4double maxEProj= GetSecondAdjEnergyMaxForScatProjToProjCase(kinEnergyScatProj);
  G4double dEmax=maxEProj-kinEnergyScatProj;
  G4double dEmin=GetLowEnergyLimit();
  G4double dE1=dEmin;
  G4double dE2=dEmin;
  
  
  std::vector< double>*  log_ESec_vector = new std::vector< double>();
  std::vector< double>*  log_Prob_vector = new std::vector< double>();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins=std::max( int(std::log10(dEmax/dEmin))*nbin_pro_decade,5);
  G4double fE=std::pow(dEmax/dEmin,1./nbins);
  
  
  
  
  
  G4double int_cross_section=0.;
  
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while (dE1 <dEmax*0.9999999999999){
  	dE2=dE1*fE;
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction2,minEProj+dE1,std::min(minEProj+dE2,maxEProj), 5);
	//G4cout<<"int_cross_section "<<minEProj+dE1<<'\t'<<int_cross_section<<G4endl;
	log_ESec_vector->push_back(std::log(std::min(dE2,maxEProj-minEProj)));
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
////////////////////////////////////////////////////////////////////////////////
//
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForSecond(      
				G4Material* aMaterial,
				G4double kinEnergyProd,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, double(G4VEmAdjointModel::*)(double)> integral;
  SelectedMaterial= aMaterial;
  kinEnergyProdForIntegration = kinEnergyProd;
   //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForProdToProjCase(kinEnergyProd);
  G4double maxEProj= GetSecondAdjEnergyMaxForProdToProjCase(kinEnergyProd);
  G4double E1=minEProj;
  std::vector< double>*  log_ESec_vector = new  std::vector< double>();
  std::vector< double>*  log_Prob_vector = new  std::vector< double>();
  log_ESec_vector->clear();
  log_Prob_vector->clear();
  log_ESec_vector->push_back(std::log(E1));
  log_Prob_vector->push_back(-50.);
  
  G4double E2=std::pow(10.,double( int(std::log10(minEProj)*nbin_pro_decade)+1)/nbin_pro_decade);
  G4double fE=std::pow(10.,1./nbin_pro_decade);
  G4double int_cross_section=0.;
  
  if (std::pow(fE,5.)>(maxEProj/minEProj)) fE = std::pow(maxEProj/minEProj,0.2);
  
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while (E1 <maxEProj*0.9999999){
  	
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction1,E1,std::min(E2,maxEProj*0.99999999), 5);
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
std::vector< std::vector<G4double>* > G4VEmAdjointModel::ComputeAdjointCrossSectionVectorPerVolumeForScatProj(
      				G4Material* aMaterial,
				G4double kinEnergyScatProj,
				G4int nbin_pro_decade) //nb bins pro order of magnitude of energy
{ G4Integrator<G4VEmAdjointModel, double(G4VEmAdjointModel::*)(double)> integral;
  SelectedMaterial= aMaterial;
  kinEnergyScatProjForIntegration = kinEnergyScatProj;
 
  //compute the vector of integrated cross sections
  //-------------------
  
  G4double minEProj= GetSecondAdjEnergyMinForScatProjToProjCase(kinEnergyScatProj);
  G4double maxEProj= GetSecondAdjEnergyMaxForScatProjToProjCase(kinEnergyScatProj);
 
  
  G4double dEmax=maxEProj-kinEnergyScatProj;
  G4double dEmin=GetLowEnergyLimit();
  G4double dE1=dEmin;
  G4double dE2=dEmin;
  
  
  std::vector< double>*  log_ESec_vector = new std::vector< double>();
  std::vector< double>*  log_Prob_vector = new std::vector< double>();
  log_ESec_vector->push_back(std::log(dEmin));
  log_Prob_vector->push_back(-50.);
  G4int nbins=std::max( int(std::log10(dEmax/dEmin))*nbin_pro_decade,5);
  G4double fE=std::pow(dEmax/dEmin,1./nbins);
  
  G4double int_cross_section=0.;
  
  // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
  while (dE1 <dEmax*0.9999999999999){
  	dE2=dE1*fE;
  	int_cross_section +=integral.Simpson(this,
	&G4VEmAdjointModel::DiffCrossSectionFunction2,minEProj+dE1,std::min(minEProj+dE2,maxEProj), 5);
	log_ESec_vector->push_back(std::log(std::min(dE2,maxEProj-minEProj)));
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
  std::vector< double>* theLogPrimEnergyVector = theMatrix->GetLogPrimEnergyVector();
  
  if (theLogPrimEnergyVector->size() ==0){
 	G4cout<<"No data are contained in the given AdjointCSMatrix!"<<G4endl;
	G4cout<<"The sampling procedure will be stopped."<<G4endl;
	return 0.;
	
  }
  
  G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
  G4double aLogPrimEnergy = std::log(aPrimEnergy);
  size_t ind =theInterpolator->FindPositionForLogVector(aLogPrimEnergy,*theLogPrimEnergyVector);
  
  
  G4double aLogPrimEnergy1,aLogPrimEnergy2;
  G4double aLogCS1,aLogCS2;
  G4double log01,log02;
  std::vector< double>* aLogSecondEnergyVector1 =0;
  std::vector< double>* aLogSecondEnergyVector2  =0;
  std::vector< double>* aLogProbVector1=0;
  std::vector< double>* aLogProbVector2=0; 
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
 	Emin=GetSecondAdjEnergyMinForScatProjToProjCase(aPrimEnergy,currentTcutForDirectSecond);
	Emax=GetSecondAdjEnergyMaxForScatProjToProjCase(aPrimEnergy);
	G4double dE=0;
	if (Emin < Emax ){
	    if (ApplyCutInRange) {
		if (second_part_of_same_type && currentTcutForDirectSecond>aPrimEnergy) return aPrimEnergy;
		
		log_rand_var1=log_rand_var+theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector1,*aLogProbVector1);
		log_rand_var2=log_rand_var+theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector2,*aLogProbVector2);
		
	    }	
	    log_dE1 = theInterpolator->Interpolate(log_rand_var1,*aLogProbVector1,*aLogSecondEnergyVector1,"Lin");
	    log_dE2 = theInterpolator->Interpolate(log_rand_var2,*aLogProbVector2,*aLogSecondEnergyVector2,"Lin");
	     dE=std::exp(theInterpolator->LinearInterpolation(aLogPrimEnergy,aLogPrimEnergy1,aLogPrimEnergy2,log_dE1,log_dE2));
	}
	
	Esec = aPrimEnergy +dE;
	Esec=std::max(Esec,Emin);
	Esec=std::min(Esec,Emax);
	
  }
  else { //Tcut condition is already full-filled
        
 	log_E1 = theInterpolator->Interpolate(log_rand_var,*aLogProbVector1,*aLogSecondEnergyVector1,"Lin");
 	log_E2 = theInterpolator->Interpolate(log_rand_var,*aLogProbVector2,*aLogSecondEnergyVector2,"Lin");
	
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
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromCSMatrix(G4double aPrimEnergy,G4bool IsScatProjToProjCase)
{ SelectCSMatrix(IsScatProjToProjCase);
  return SampleAdjSecEnergyFromCSMatrix(indexOfUsedCrossSectionMatrix, aPrimEnergy, IsScatProjToProjCase);
}
//////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::SelectCSMatrix(G4bool IsScatProjToProjCase)
{ 
  indexOfUsedCrossSectionMatrix=0;
  if (!UseMatrixPerElement) indexOfUsedCrossSectionMatrix = currentMaterialIndex;
  else if (!UseOnlyOneMatrixForAllElements) { //Select Material
   	std::vector<G4double>* CS_Vs_Element = &CS_Vs_ElementForScatProjToProjCase;
	lastCS=lastAdjointCSForScatProjToProjCase;
  	if ( !IsScatProjToProjCase) {
		CS_Vs_Element = &CS_Vs_ElementForProdToProjCase;
		lastCS=lastAdjointCSForProdToProjCase;
	}	
  	G4double rand_var= G4UniformRand();
  	G4double SumCS=0.;
	size_t ind=0;
  	for (size_t i=0;i<CS_Vs_Element->size();i++){
 		SumCS+=(*CS_Vs_Element)[i];
		if (rand_var<=SumCS/lastCS){
			ind=i;
			break;
		}
  	}
	indexOfUsedCrossSectionMatrix = currentMaterial->GetElement(ind)->GetIndex();
  }
}
//////////////////////////////////////////////////////////////////////////////
//				
G4double G4VEmAdjointModel::SampleAdjSecEnergyFromDiffCrossSectionPerAtom(G4double prim_energy,G4bool IsScatProjToProjCase)
{  
  // here we try to use the rejection method 
  //-----------------------------------------
  
  const G4int iimax = 1000;
  G4double E=0;
  G4double x,xmin,greject,q;
  if ( IsScatProjToProjCase){
  	G4double Emax = GetSecondAdjEnergyMaxForScatProjToProjCase(prim_energy);
        G4double Emin= prim_energy+currentTcutForDirectSecond;
	xmin=Emin/Emax;
	G4double grejmax = DiffCrossSectionPerAtomPrimToScatPrim(Emin,prim_energy,1)*prim_energy; 

        G4int ii =0;
	do {
        	q = G4UniformRand();
        	x = 1./(q*(1./xmin -1.) +1.);
		E=x*Emax;
		greject = DiffCrossSectionPerAtomPrimToScatPrim( E,prim_energy ,1)*prim_energy; 
		++ii;
                if(ii >= iimax) { break; }
	}	
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	while( greject < G4UniformRand()*grejmax );	
	
  }
  else {
  	G4double Emax = GetSecondAdjEnergyMaxForProdToProjCase(prim_energy);
        G4double Emin=  GetSecondAdjEnergyMinForProdToProjCase(prim_energy);;
	xmin=Emin/Emax;
	G4double grejmax = DiffCrossSectionPerAtomPrimToSecond(Emin,prim_energy,1); 
        G4int ii =0;
	do {
        	q = G4UniformRand();
        	x = std::pow(xmin, q);
		E=x*Emax;
		greject = DiffCrossSectionPerAtomPrimToSecond( E,prim_energy ,1); 
		++ii;
                if(ii >= iimax) { break; }		
	}	
	// Loop checking, 07-Aug-2015, Vladimir Ivanchenko
	while( greject < G4UniformRand()*grejmax );
  
  
  
  }
  
  return E;
}

////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::CorrectPostStepWeight(G4ParticleChange* fParticleChange, 
						  G4double old_weight,  
						  G4double adjointPrimKinEnergy, 
						  G4double projectileKinEnergy,
						  G4bool IsScatProjToProjCase) 
{
 G4double new_weight=old_weight;
 G4double w_corr =1./CS_biasing_factor;
 w_corr*=G4AdjointCSManager::GetAdjointCSManager()->GetPostStepWeightCorrection();
 
 
 lastCS=lastAdjointCSForScatProjToProjCase;
 if ( !IsScatProjToProjCase) lastCS=lastAdjointCSForProdToProjCase;
 if ((adjointPrimKinEnergy-preStepEnergy)/preStepEnergy>0.001){ //Is that in all cases needed???
 	G4double post_stepCS=AdjointCrossSection(currentCouple, adjointPrimKinEnergy
						 ,IsScatProjToProjCase );
	if (post_stepCS>0 && lastCS>0) w_corr*=post_stepCS/lastCS;
 }
	
 new_weight*=w_corr;

 //G4cout<<"Post step "<<new_weight<<'\t'<<w_corr<<'\t'<<old_weight<<G4endl;
 new_weight*=projectileKinEnergy/adjointPrimKinEnergy;//This is needed due to the biasing of diff CS
 							//by the factor adjointPrimKinEnergy/projectileKinEnergy
   


 fParticleChange->SetParentWeightByProcess(false);
 fParticleChange->SetSecondaryWeightByProcess(false);
 fParticleChange->ProposeParentWeight(new_weight);	
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
{ G4double Emin=PrimAdjEnergy;
  if (ApplyCutInRange) Emin=PrimAdjEnergy+Tcut;
  return Emin;
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
    	currentTcutForDirectSecond =0.00000000001;
   	if (theAdjEquivOfDirectSecondPartDef) {
    		if (theAdjEquivOfDirectSecondPartDef == G4AdjointGamma::AdjointGamma()) idx = 0;
   		else if (theAdjEquivOfDirectSecondPartDef == G4AdjointElectron::AdjointElectron()) idx = 1;
    		else if (theAdjEquivOfDirectSecondPartDef == G4AdjointPositron::AdjointPositron()) idx = 2;
    		if (idx <56){
			const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
    			currentTcutForDirectSecond=(*aVec)[currentCoupleIndex];
		}	
    	}	
    
    	
  }
}
////////////////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::SetHighEnergyLimit(G4double aVal)
{ HighEnergyLimit=aVal;
  if (theDirectEMModel) theDirectEMModel->SetHighEnergyLimit( aVal);
}
////////////////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::SetLowEnergyLimit(G4double aVal)
{
  LowEnergyLimit=aVal;
  if (theDirectEMModel) theDirectEMModel->SetLowEnergyLimit( aVal);
}
////////////////////////////////////////////////////////////////////////////////////////////
//
void G4VEmAdjointModel::SetAdjointEquivalentOfDirectPrimaryParticleDefinition(G4ParticleDefinition* aPart)
{
  theAdjEquivOfDirectPrimPartDef=aPart;
  if (theAdjEquivOfDirectPrimPartDef->GetParticleName() =="adj_e-")
					theDirectPrimaryPartDef=G4Electron::Electron();
  if (theAdjEquivOfDirectPrimPartDef->GetParticleName() =="adj_gamma")
					theDirectPrimaryPartDef=G4Gamma::Gamma();
}
