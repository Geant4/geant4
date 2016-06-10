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
// $Id: G4AdjointCSManager.cc 93569 2015-10-26 14:53:21Z gcosmo $
//

#include <fstream>
#include <iomanip>

#include "G4AdjointCSManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4AdjointInterpolator.hh"
#include "G4AdjointCSMatrix.hh"
#include "G4VEmAdjointModel.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4ParticleDefinition.hh"
#include "G4Element.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4PhysicsTable.hh" 
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4AdjointProton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCutsTable.hh"

G4ThreadLocal G4AdjointCSManager* G4AdjointCSManager::theInstance = nullptr;
///////////////////////////////////////////////////////
//
G4AdjointCSManager* G4AdjointCSManager::GetAdjointCSManager()
{ 
  if(theInstance == nullptr) {
    static G4ThreadLocalSingleton<G4AdjointCSManager> inst;
    theInstance = inst.Instance();
  }
  return theInstance; 
}

///////////////////////////////////////////////////////
//
G4AdjointCSManager::G4AdjointCSManager()
{ CrossSectionMatrixesAreBuilt=false;
  TotalSigmaTableAreBuilt=false;
  theTotalForwardSigmaTableVector.clear();
  theTotalAdjointSigmaTableVector.clear();
  listOfForwardEmProcess.clear();
  listOfForwardEnergyLossProcess.clear();
  theListOfAdjointParticlesInAction.clear(); 
  EminForFwdSigmaTables.clear();
  EminForAdjSigmaTables.clear();
  EkinofFwdSigmaMax.clear();
  EkinofAdjSigmaMax.clear();
  listSigmaTableForAdjointModelScatProjToProj.clear();
  listSigmaTableForAdjointModelProdToProj.clear();
  Tmin=0.1*keV;
  Tmax=100.*TeV;
  nbins=320; //probably this should be decrease, that was choosen to avoid error in the CS value closed to CS jump.(For example at Tcut)
  
  RegisterAdjointParticle(G4AdjointElectron::AdjointElectron());
  RegisterAdjointParticle(G4AdjointGamma::AdjointGamma());
  RegisterAdjointParticle(G4AdjointProton::AdjointProton());
  
  verbose  = 1;
  currentParticleIndex = 0;
  currentMatIndex = 0;
  eindex = 0;
 
  lastPartDefForCS = nullptr;
  LastEkinForCS = lastPrimaryEnergy = lastTcut = 0.;
  LastCSCorrectionFactor = massRatio = 1.;
  
  forward_CS_is_used = true;
  forward_CS_mode = true;
  
  currentParticleDef = nullptr;
  currentCouple =nullptr;
  currentMaterial=nullptr;
  lastMaterial=nullptr;

  theAdjIon = nullptr;
  theFwdIon = nullptr;
 
  PreadjCS = PostadjCS = PrefwdCS = PostfwdCS = 0.0;
}
///////////////////////////////////////////////////////
//
G4AdjointCSManager::~G4AdjointCSManager()
{;
}
///////////////////////////////////////////////////////
//
size_t G4AdjointCSManager::RegisterEmAdjointModel(G4VEmAdjointModel* aModel)
{listOfAdjointEMModel.push_back(aModel);
 listSigmaTableForAdjointModelScatProjToProj.push_back(new G4PhysicsTable);
 listSigmaTableForAdjointModelProdToProj.push_back(new G4PhysicsTable);
 return listOfAdjointEMModel.size() -1;
 
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterEmProcess(G4VEmProcess* aProcess, G4ParticleDefinition* aFwdPartDef)
{ 
  G4ParticleDefinition* anAdjPartDef = GetAdjointParticleEquivalent(aFwdPartDef);
  if (anAdjPartDef && aProcess){
  	RegisterAdjointParticle(anAdjPartDef);
	G4int index=-1;
	
  	for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  		if (anAdjPartDef->GetParticleName() == theListOfAdjointParticlesInAction[i]->GetParticleName()) index=i;
  	}
  	listOfForwardEmProcess[index]->push_back(aProcess);
  }
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterEnergyLossProcess(G4VEnergyLossProcess* aProcess, G4ParticleDefinition* aFwdPartDef)
{
  G4ParticleDefinition* anAdjPartDef = GetAdjointParticleEquivalent(aFwdPartDef);
  if (anAdjPartDef && aProcess){
        RegisterAdjointParticle(anAdjPartDef);
  	G4int index=-1;
  	for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  		if (anAdjPartDef->GetParticleName() == theListOfAdjointParticlesInAction[i]->GetParticleName()) index=i;
  	}
  	listOfForwardEnergyLossProcess[index]->push_back(aProcess);
   }
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterAdjointParticle(G4ParticleDefinition* aPartDef)
{  G4int index=-1;
   for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	if (aPartDef->GetParticleName() == theListOfAdjointParticlesInAction[i]->GetParticleName()) index=i;
   }
  
   if (index ==-1){
  	listOfForwardEnergyLossProcess.push_back(new std::vector<G4VEnergyLossProcess*>());
	theTotalForwardSigmaTableVector.push_back(new G4PhysicsTable);
	theTotalAdjointSigmaTableVector.push_back(new G4PhysicsTable);
	listOfForwardEmProcess.push_back(new std::vector<G4VEmProcess*>());
	theListOfAdjointParticlesInAction.push_back(aPartDef);
	EminForFwdSigmaTables.push_back(std::vector<G4double> ());
	EminForAdjSigmaTables.push_back(std::vector<G4double> ());
	EkinofFwdSigmaMax.push_back(std::vector<G4double> ());
	EkinofAdjSigmaMax.push_back(std::vector<G4double> ());
	
   }
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::BuildCrossSectionMatrices()
{	
	if (CrossSectionMatrixesAreBuilt) return;
		//Tcut, Tmax 
			//The matrices will be computed probably just once
			  //When Tcut will change some PhysicsTable will be recomputed  
			  // for each MaterialCutCouple but not all the matrices	
			  //The Tcut defines a lower limit in the energy of the Projectile before the scattering
			  //In the Projectile to Scattered Projectile case we have
			  // 			E_ScatProj<E_Proj-Tcut
			  //Therefore in the adjoint case we have
			  //			Eproj> E_ScatProj+Tcut
			  //This implies that when computing the adjoint CS we should integrate over Epro
			  // from E_ScatProj+Tcut to Emax
			  //In the Projectile to Secondary case Tcut plays a role only in the fact that  
			  // Esecond should be greater than Tcut to have the possibility to have any adjoint
			  //process 		 
			  //To avoid to recompute the matrices for all changes of MaterialCutCouple
			  //We propose to compute the matrices only once for the minimum possible Tcut and then
			  //to interpolate the probability for a new Tcut (implemented in G4VAdjointEmModel)
	
	
	theAdjointCSMatricesForScatProjToProj.clear();
	theAdjointCSMatricesForProdToProj.clear();
	const G4ElementTable* theElementTable = G4Element::GetElementTable();
	const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
	
	G4cout<<"========== Computation of cross section matrices for adjoint models =========="<<G4endl;
	for (size_t i=0; i<listOfAdjointEMModel.size();i++){
		G4VEmAdjointModel* aModel =listOfAdjointEMModel[i];
		G4cout<<"Build adjoint cross section matrices for "<<aModel->GetName()<<G4endl;
		if (aModel->GetUseMatrix()){
			std::vector<G4AdjointCSMatrix*>* aListOfMat1 = new std::vector<G4AdjointCSMatrix*>();
			std::vector<G4AdjointCSMatrix*>* aListOfMat2 = new std::vector<G4AdjointCSMatrix*>();
			aListOfMat1->clear();
			aListOfMat2->clear();
			if (aModel->GetUseMatrixPerElement()){
				if (aModel->GetUseOnlyOneMatrixForAllElements()){
						std::vector<G4AdjointCSMatrix*>
						two_matrices=BuildCrossSectionsMatricesForAGivenModelAndElement(aModel,1, 1, 80);
						aListOfMat1->push_back(two_matrices[0]);
						aListOfMat2->push_back(two_matrices[1]);
				}
				else {		
					for (size_t j=0; j<theElementTable->size();j++){
						G4Element* anElement=(*theElementTable)[j];
						G4int Z = G4lrint(anElement->GetZ());
						G4int A = G4lrint(anElement->GetN());
						std::vector<G4AdjointCSMatrix*>
							two_matrices=BuildCrossSectionsMatricesForAGivenModelAndElement(aModel,Z, A, 40);
						aListOfMat1->push_back(two_matrices[0]);
						aListOfMat2->push_back(two_matrices[1]);
					}
				}	
			}
			else { //Per material case
				for (size_t j=0; j<theMaterialTable->size();j++){
					G4Material* aMaterial=(*theMaterialTable)[j];
					std::vector<G4AdjointCSMatrix*>
						two_matrices=BuildCrossSectionsMatricesForAGivenModelAndMaterial(aModel,aMaterial, 40);
					aListOfMat1->push_back(two_matrices[0]);
					aListOfMat2->push_back(two_matrices[1]);
				}
			
			}
			theAdjointCSMatricesForProdToProj.push_back(*aListOfMat1);
			theAdjointCSMatricesForScatProjToProj.push_back(*aListOfMat2);	
			aModel->SetCSMatrices(aListOfMat1, aListOfMat2); 	
		}
		else {  G4cout<<"The model "<<aModel->GetName()<<" does not use cross section matrices"<<G4endl;
			std::vector<G4AdjointCSMatrix*> two_empty_matrices;
			theAdjointCSMatricesForProdToProj.push_back(two_empty_matrices);
			theAdjointCSMatricesForScatProjToProj.push_back(two_empty_matrices);
			
		}		
	}
	G4cout<<"              All adjoint cross section matrices are computed!"<<G4endl;
	G4cout<<"======================================================================"<<G4endl;
	
	CrossSectionMatrixesAreBuilt = true;


}


///////////////////////////////////////////////////////
//
void G4AdjointCSManager::BuildTotalSigmaTables()
{ if (TotalSigmaTableAreBuilt) return;

  
  const G4ProductionCutsTable* theCoupleTable= G4ProductionCutsTable::GetProductionCutsTable();
 
 
 //Prepare the Sigma table for all AdjointEMModel, will be filled later on 
  for (size_t i=0; i<listOfAdjointEMModel.size();i++){
  	listSigmaTableForAdjointModelScatProjToProj[i]->clearAndDestroy();
	listSigmaTableForAdjointModelProdToProj[i]->clearAndDestroy();
	for (size_t j=0;j<theCoupleTable->GetTableSize();j++){
		listSigmaTableForAdjointModelScatProjToProj[i]->push_back(new G4PhysicsLogVector(Tmin, Tmax, nbins));
		listSigmaTableForAdjointModelProdToProj[i]->push_back(new G4PhysicsLogVector(Tmin, Tmax, nbins));
	}
  }
  


  for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	G4ParticleDefinition* thePartDef = theListOfAdjointParticlesInAction[i];
	DefineCurrentParticle(thePartDef);
  	theTotalForwardSigmaTableVector[i]->clearAndDestroy();
	theTotalAdjointSigmaTableVector[i]->clearAndDestroy();
	EminForFwdSigmaTables[i].clear();
	EminForAdjSigmaTables[i].clear();
	EkinofFwdSigmaMax[i].clear();
	EkinofAdjSigmaMax[i].clear();
	//G4cout<<thePartDef->GetParticleName();
	
  	for (size_t j=0;j<theCoupleTable->GetTableSize();j++){
		const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
		
		/*
		G4String file_name1=couple->GetMaterial()->GetName()+"_"+thePartDef->GetParticleName()+"_adj_totCS.txt";
		G4String file_name2=couple->GetMaterial()->GetName()+"_"+thePartDef->GetParticleName()+"_fwd_totCS.txt";
		
		std::fstream FileOutputAdjCS(file_name1, std::ios::out);
		std::fstream FileOutputFwdCS(file_name2, std::ios::out);
		
		
		
		FileOutputAdjCS<<std::setiosflags(std::ios::scientific);
 		FileOutputAdjCS<<std::setprecision(6);
		FileOutputFwdCS<<std::setiosflags(std::ios::scientific);
 		FileOutputFwdCS<<std::setprecision(6);	
		*/
			  

		//make first the total fwd CS table for FwdProcess
		G4PhysicsVector* aVector =  new G4PhysicsLogVector(Tmin, Tmax, nbins);
		G4bool Emin_found=false;
		G4double sigma_max =0.;
		G4double e_sigma_max =0.;
		for(size_t l=0; l<aVector->GetVectorLength(); l++) { 
			G4double totCS=0.;
			G4double e=aVector->GetLowEdgeEnergy(l);
			for (size_t k=0; k<listOfForwardEmProcess[i]->size(); k++){
				totCS+=(*listOfForwardEmProcess[i])[k]->GetLambda(e, couple);
			}
			for (size_t k=0; k<listOfForwardEnergyLossProcess[i]->size(); k++){
				if (thePartDef == theAdjIon) { // e is considered already as the scaled energy
				        size_t mat_index = couple->GetIndex();
					G4VEmModel* currentModel = (*listOfForwardEnergyLossProcess[i])[k]->SelectModelForMaterial(e,mat_index);
  					G4double chargeSqRatio =  currentModel->GetChargeSquareRatio(theFwdIon,couple->GetMaterial(),e/massRatio);
					(*listOfForwardEnergyLossProcess[i])[k]->SetDynamicMassCharge(massRatio,chargeSqRatio);
				}
				G4double e1=e/massRatio;
				totCS+=(*listOfForwardEnergyLossProcess[i])[k]->GetLambda(e1, couple);
			}
			aVector->PutValue(l,totCS);
			if (totCS>sigma_max){
				sigma_max=totCS;
				e_sigma_max = e;
				
			}
			//FileOutputFwdCS<<e<<'\t'<<totCS<<G4endl;
			
			if (totCS>0 && !Emin_found) {
				EminForFwdSigmaTables[i].push_back(e);
				Emin_found=true;
			}
			
		
		}
		//FileOutputFwdCS.close();
		
		EkinofFwdSigmaMax[i].push_back(e_sigma_max);
		
		
		if(!Emin_found)	EminForFwdSigmaTables[i].push_back(Tmax);
		
		theTotalForwardSigmaTableVector[i]->push_back(aVector);
		
		
		Emin_found=false;
		sigma_max=0;
		e_sigma_max =0.;
		G4PhysicsVector* aVector1 =  new G4PhysicsLogVector(Tmin, Tmax, nbins);
		for(eindex=0; eindex<aVector->GetVectorLength(); eindex++) { 
			G4double e=aVector->GetLowEdgeEnergy(eindex);
			G4double totCS =ComputeTotalAdjointCS(couple,thePartDef,e*0.9999999/massRatio); //massRatio needed for ions
			aVector1->PutValue(eindex,totCS);
			if (totCS>sigma_max){
				sigma_max=totCS;
				e_sigma_max = e;
				
			}
			//FileOutputAdjCS<<e<<'\t'<<totCS<<G4endl;
			if (totCS>0 && !Emin_found) {
				EminForAdjSigmaTables[i].push_back(e);
				Emin_found=true;
			}
		
		}
		//FileOutputAdjCS.close();
		EkinofAdjSigmaMax[i].push_back(e_sigma_max);
		if(!Emin_found)	EminForAdjSigmaTables[i].push_back(Tmax);
			
		theTotalAdjointSigmaTableVector[i]->push_back(aVector1);
		
	}
  }
  TotalSigmaTableAreBuilt =true;
   
}
///////////////////////////////////////////////////////
//
G4double G4AdjointCSManager::GetTotalAdjointCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple)
{ DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);	
  G4bool b;
  return (((*theTotalAdjointSigmaTableVector[currentParticleIndex])[currentMatIndex])->GetValue(Ekin*massRatio, b));
  
  
  
}				     
///////////////////////////////////////////////////////
//				     
G4double G4AdjointCSManager::GetTotalForwardCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple)
{ DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  G4bool b;
  return (((*theTotalForwardSigmaTableVector[currentParticleIndex])[currentMatIndex])->GetValue(Ekin*massRatio, b));
  
  
}
///////////////////////////////////////////////////////
//
G4double G4AdjointCSManager::GetAdjointSigma(G4double Ekin_nuc, size_t index_model,G4bool is_scat_proj_to_proj,
	 		   	     const G4MaterialCutsCouple* aCouple)
{ DefineCurrentMaterial(aCouple);
  G4bool b;
  if (is_scat_proj_to_proj) return (((*listSigmaTableForAdjointModelScatProjToProj[index_model])[currentMatIndex])->GetValue(Ekin_nuc, b));
  else return (((*listSigmaTableForAdjointModelProdToProj[index_model])[currentMatIndex])->GetValue(Ekin_nuc, b));
}				     				     
///////////////////////////////////////////////////////
//				     
void G4AdjointCSManager::GetEminForTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double& emin_adj, G4double& emin_fwd)
{ DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  emin_adj = EminForAdjSigmaTables[currentParticleIndex][currentMatIndex]/massRatio;
  emin_fwd = EminForFwdSigmaTables[currentParticleIndex][currentMatIndex]/massRatio;
  
  
  
}
///////////////////////////////////////////////////////
//				     
void G4AdjointCSManager::GetMaxFwdTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double& e_sigma_max, G4double& sigma_max)
{ DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  e_sigma_max = EkinofFwdSigmaMax[currentParticleIndex][currentMatIndex];
  G4bool b;
  sigma_max =((*theTotalForwardSigmaTableVector[currentParticleIndex])[currentMatIndex])->GetValue(e_sigma_max, b);
  e_sigma_max/=massRatio;
  
  
}
///////////////////////////////////////////////////////
//				     
void G4AdjointCSManager::GetMaxAdjTotalCS(G4ParticleDefinition* aPartDef,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double& e_sigma_max, G4double& sigma_max)
{ DefineCurrentMaterial(aCouple);
  DefineCurrentParticle(aPartDef);
  e_sigma_max = EkinofAdjSigmaMax[currentParticleIndex][currentMatIndex];
  G4bool b;
  sigma_max =((*theTotalAdjointSigmaTableVector[currentParticleIndex])[currentMatIndex])->GetValue(e_sigma_max, b);
  e_sigma_max/=massRatio;
  
  
}				     
///////////////////////////////////////////////////////
//				     
G4double G4AdjointCSManager::GetCrossSectionCorrection(G4ParticleDefinition* aPartDef,G4double PreStepEkin,const G4MaterialCutsCouple* aCouple, G4bool& fwd_is_used,
										G4double& fwd_TotCS)
{ G4double corr_fac = 1.;
  if (forward_CS_mode && aPartDef ) {
  	fwd_TotCS=PrefwdCS;
  	if (LastEkinForCS != PreStepEkin || aPartDef != lastPartDefForCS || aCouple!=currentCouple) {
		DefineCurrentMaterial(aCouple);
		PreadjCS = GetTotalAdjointCS(aPartDef, PreStepEkin,aCouple);
		PrefwdCS = GetTotalForwardCS(aPartDef, PreStepEkin,aCouple);
		LastEkinForCS = PreStepEkin;
		lastPartDefForCS = aPartDef;
		if (PrefwdCS >0. &&  PreadjCS >0.) {
			forward_CS_is_used = true;
			LastCSCorrectionFactor = PrefwdCS/PreadjCS;
		}
		else {
			forward_CS_is_used = false;
			LastCSCorrectionFactor = 1.;
			
		}
		
	}
	corr_fac =LastCSCorrectionFactor;
	
	
	
  }  
  else {
  	forward_CS_is_used = false;
  	LastCSCorrectionFactor = 1.;
  }
  fwd_TotCS=PrefwdCS;	
  fwd_is_used = forward_CS_is_used;
  return  corr_fac;
}				     
///////////////////////////////////////////////////////
//				     			     
G4double G4AdjointCSManager::GetContinuousWeightCorrection(G4ParticleDefinition* aPartDef, G4double PreStepEkin,G4double AfterStepEkin,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double step_length)
{  G4double corr_fac = 1.;
  //return corr_fac;
  //G4double after_adjCS = GetTotalAdjointCS(aPartDef, AfterStepEkin,aCouple);
  G4double after_fwdCS = GetTotalForwardCS(aPartDef, AfterStepEkin,aCouple);
  G4double pre_adjCS = GetTotalAdjointCS(aPartDef, PreStepEkin,aCouple);
  if (!forward_CS_is_used || pre_adjCS ==0. ||  after_fwdCS==0.) {
	forward_CS_is_used=false;
	G4double pre_fwdCS = GetTotalForwardCS(aPartDef, PreStepEkin,aCouple);
	corr_fac *=std::exp((pre_adjCS-pre_fwdCS)*step_length);
	LastCSCorrectionFactor = 1.;
  }
  else {
	LastCSCorrectionFactor = after_fwdCS/pre_adjCS;
  }	
	

 
  return corr_fac; 
}				     
///////////////////////////////////////////////////////
//				     			     
G4double G4AdjointCSManager::GetPostStepWeightCorrection( )
{//return 1.; 
 return  1./LastCSCorrectionFactor;
	
}							  
///////////////////////////////////////////////////////
//
G4double  G4AdjointCSManager::ComputeAdjointCS(G4Material* aMaterial,
					    		 G4VEmAdjointModel* aModel, 
					    		 G4double PrimEnergy,
					    		 G4double Tcut,
					    		 G4bool IsScatProjToProjCase,
							 std::vector<G4double>& CS_Vs_Element)
{ 
  
  G4double EminSec=0;
  G4double EmaxSec=0;
  
  if (IsScatProjToProjCase){
	EminSec= aModel->GetSecondAdjEnergyMinForScatProjToProjCase(PrimEnergy,Tcut);
	EmaxSec= aModel->GetSecondAdjEnergyMaxForScatProjToProjCase(PrimEnergy);
  }
  else if (PrimEnergy > Tcut || !aModel->GetApplyCutInRange()) {
	EminSec= aModel->GetSecondAdjEnergyMinForProdToProjCase(PrimEnergy);
	EmaxSec= aModel->GetSecondAdjEnergyMaxForProdToProjCase(PrimEnergy);
  }
  if (EminSec >= EmaxSec) return 0.;


  G4bool need_to_compute=false;
  if ( aMaterial!= lastMaterial || PrimEnergy != lastPrimaryEnergy || Tcut != lastTcut){
  	lastMaterial =aMaterial;
  	lastPrimaryEnergy = PrimEnergy;
  	lastTcut=Tcut;
	listOfIndexOfAdjointEMModelInAction.clear();
	listOfIsScatProjToProjCase.clear();
	lastAdjointCSVsModelsAndElements.clear();
	need_to_compute=true;
	
  }
  size_t ind=0;
  if (!need_to_compute){
  	need_to_compute=true;
  	for (size_t i=0;i<listOfIndexOfAdjointEMModelInAction.size();i++){
		size_t ind1=listOfIndexOfAdjointEMModelInAction[i];
		if (aModel == listOfAdjointEMModel[ind1] && IsScatProjToProjCase == listOfIsScatProjToProjCase[i]){
			need_to_compute=false;
			CS_Vs_Element = lastAdjointCSVsModelsAndElements[ind];
		}
		ind++;
	}
  }
  
  if (need_to_compute){
  	size_t ind_model=0;
	for (size_t i=0;i<listOfAdjointEMModel.size();i++){
		if (aModel == listOfAdjointEMModel[i]){
			ind_model=i;
			break;
		}
	}
	G4double Tlow=Tcut;
	if (!listOfAdjointEMModel[ind_model]->GetApplyCutInRange()) Tlow =listOfAdjointEMModel[ind_model]->GetLowEnergyLimit();
	listOfIndexOfAdjointEMModelInAction.push_back(ind_model);	
	listOfIsScatProjToProjCase.push_back(IsScatProjToProjCase);
	CS_Vs_Element.clear();
	if (!aModel->GetUseMatrix()){
		CS_Vs_Element.push_back(aModel->AdjointCrossSection(currentCouple,PrimEnergy,IsScatProjToProjCase));
				         
	
	}
	else if (aModel->GetUseMatrixPerElement()){
			size_t n_el = aMaterial->GetNumberOfElements();
		if (aModel->GetUseOnlyOneMatrixForAllElements()){
			G4AdjointCSMatrix* theCSMatrix;
			if (IsScatProjToProjCase){
				theCSMatrix=theAdjointCSMatricesForScatProjToProj[ind_model][0];
			}
			else  	theCSMatrix=theAdjointCSMatricesForProdToProj[ind_model][0];
			G4double CS =0.;
			if (PrimEnergy > Tlow)
					CS = ComputeAdjointCS(PrimEnergy,theCSMatrix,Tlow);
			G4double factor=0.;
			for (size_t i=0;i<n_el;i++){ //this could be computed only once
				//size_t ind_el = aMaterial->GetElement(i)->GetIndex();
				factor+=aMaterial->GetElement(i)->GetZ()*aMaterial->GetVecNbOfAtomsPerVolume()[i];
			}
			CS *=factor;
			CS_Vs_Element.push_back(CS);
								
		}
		else {
			for (size_t i=0;i<n_el;i++){
				size_t ind_el = aMaterial->GetElement(i)->GetIndex();
				//G4cout<<aMaterial->GetName()<<G4endl;
				G4AdjointCSMatrix* theCSMatrix;
				if (IsScatProjToProjCase){
					theCSMatrix=theAdjointCSMatricesForScatProjToProj[ind_model][ind_el];
				}
				else  	theCSMatrix=theAdjointCSMatricesForProdToProj[ind_model][ind_el];
				G4double CS =0.;
				if (PrimEnergy > Tlow)
					CS = ComputeAdjointCS(PrimEnergy,theCSMatrix,Tlow);
				//G4cout<<CS<<G4endl;			
				CS_Vs_Element.push_back(CS*(aMaterial->GetVecNbOfAtomsPerVolume()[i]));	
			}
		}		
		
	}
	else {
		size_t ind_mat = aMaterial->GetIndex();
		G4AdjointCSMatrix* theCSMatrix;
		if (IsScatProjToProjCase){
			theCSMatrix=theAdjointCSMatricesForScatProjToProj[ind_model][ind_mat];
		}
		else  	theCSMatrix=theAdjointCSMatricesForProdToProj[ind_model][ind_mat];
		G4double CS =0.;
		if (PrimEnergy > Tlow)
			CS = ComputeAdjointCS(PrimEnergy,theCSMatrix,Tlow);
		CS_Vs_Element.push_back(CS);							
			
		
	}
	lastAdjointCSVsModelsAndElements.push_back(CS_Vs_Element);
	
  }
  
  
  G4double CS=0;
  for (size_t i=0;i<CS_Vs_Element.size();i++){
  	CS+=CS_Vs_Element[i]; //We could put the progressive sum of the CS instead of the CS of an element itself
  	
  }
  return CS;
}							
///////////////////////////////////////////////////////
//	
G4Element* G4AdjointCSManager::SampleElementFromCSMatrices(G4Material* aMaterial,
						       		   G4VEmAdjointModel* aModel,
					 	       		   G4double PrimEnergy,
						        	   G4double Tcut,
						       		   G4bool IsScatProjToProjCase)
{ std::vector<G4double> CS_Vs_Element;
  G4double CS = ComputeAdjointCS(aMaterial,aModel,PrimEnergy,Tcut,IsScatProjToProjCase,CS_Vs_Element);
  G4double rand_var= G4UniformRand();
  G4double SumCS=0.;
  size_t ind=0;
  for (size_t i=0;i<CS_Vs_Element.size();i++){
 	SumCS+=CS_Vs_Element[i];
	if (rand_var<=SumCS/CS){
		ind=i;
		break;
	}
  }
 
  return const_cast<G4Element*>(aMaterial->GetElement(ind));
 
 
					    
}								
///////////////////////////////////////////////////////
//
G4double G4AdjointCSManager::ComputeTotalAdjointCS(const G4MaterialCutsCouple* aCouple,
							     G4ParticleDefinition* aPartDef,
							     G4double Ekin)
{
 G4double TotalCS=0.;

 DefineCurrentMaterial(aCouple);

  
 std::vector<G4double> CS_Vs_Element;
 G4double CS;
 for (size_t i=0; i<listOfAdjointEMModel.size();i++){
 	
	G4double Tlow=0;
	if (!listOfAdjointEMModel[i]->GetApplyCutInRange()) Tlow =listOfAdjointEMModel[i]->GetLowEnergyLimit();
	else {
		G4ParticleDefinition* theDirSecondPartDef = 
			GetForwardParticleEquivalent(listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectSecondaryParticleDefinition());
		size_t idx=56;
 		if (theDirSecondPartDef->GetParticleName() == "gamma") idx = 0;
		else if (theDirSecondPartDef->GetParticleName() == "e-") idx = 1;
 		else if (theDirSecondPartDef->GetParticleName() == "e+") idx = 2;
		if (idx <56) {
			const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
 			Tlow =(*aVec)[aCouple->GetIndex()];
		}	
		
	
	}
 	if ( Ekin<=listOfAdjointEMModel[i]->GetHighEnergyLimit() && Ekin>=listOfAdjointEMModel[i]->GetLowEnergyLimit()){
		if (aPartDef == listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectPrimaryParticleDefinition()){
			CS=ComputeAdjointCS(currentMaterial,
					    	       listOfAdjointEMModel[i], 
					    	       Ekin, Tlow,true,CS_Vs_Element);
			TotalCS += CS;
			(*listSigmaTableForAdjointModelScatProjToProj[i])[currentMatIndex]->PutValue(eindex,CS);			       
		}
		if (aPartDef == listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectSecondaryParticleDefinition()){
			CS = ComputeAdjointCS(currentMaterial,
					    	       listOfAdjointEMModel[i], 
					    	       Ekin, Tlow,false, CS_Vs_Element);
			TotalCS += CS;
			(*listSigmaTableForAdjointModelProdToProj[i])[currentMatIndex]->PutValue(eindex,CS);
		}
		
	}
	else {
		(*listSigmaTableForAdjointModelScatProjToProj[i])[currentMatIndex]->PutValue(eindex,0.);
		(*listSigmaTableForAdjointModelProdToProj[i])[currentMatIndex]->PutValue(eindex,0.);
		
	}
 }
 return TotalCS;
    
 
}	
///////////////////////////////////////////////////////
//
std::vector<G4AdjointCSMatrix*>
G4AdjointCSManager::BuildCrossSectionsMatricesForAGivenModelAndElement(G4VEmAdjointModel* aModel,G4int Z,G4int A,
								        G4int nbin_pro_decade)
{ 
  G4AdjointCSMatrix* theCSMatForProdToProjBackwardScattering = new G4AdjointCSMatrix(false);
  G4AdjointCSMatrix* theCSMatForScatProjToProjBackwardScattering = new G4AdjointCSMatrix(true);
 
 
  //make the vector of primary energy of the adjoint particle, could try to make this just once ? 
  
   G4double EkinMin =aModel->GetLowEnergyLimit();
   G4double EkinMaxForScat =aModel->GetHighEnergyLimit()*0.999;
   G4double EkinMaxForProd =aModel->GetHighEnergyLimit()*0.999;
   if (aModel->GetSecondPartOfSameType() )EkinMaxForProd =EkinMaxForProd/2.;
   
    
   //Product to projectile backward scattering
   //-----------------------------------------
   G4double fE=std::pow(10.,1./nbin_pro_decade);
   G4double E2=std::pow(10.,double( int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   G4double E1=EkinMin;
   // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
   while (E1 <EkinMaxForProd){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForProd,E1);
	std::vector< std::vector< double>* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerAtomForSecond(E1,Z,A,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< double>* log_ESecVec=aMat[0];
		std::vector< double>* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
		//normalise CSVec such that it becomes a probability vector
 		for (size_t j=0;j<log_CSVec->size();j++) {
	        	if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS) +1e-50);
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-std::log(1000.);
		theCSMatForProdToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
	}	
   	E1=E2;
   	E2*=fE;
   }
   
   //Scattered projectile to projectile backward scattering
   //-----------------------------------------
   
   E2=std::pow(10.,double( int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   E1=EkinMin;
   // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
   while (E1 <EkinMaxForScat){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForScat,E1);
	std::vector< std::vector< double>* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerAtomForScatProj(E1,Z,A,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< double>* log_ESecVec=aMat[0];
		std::vector< double>* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
		//normalise CSVec such that it becomes a probability vector
 		for (size_t j=0;j<log_CSVec->size();j++) {
	        	if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS)+1e-50);
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-std::log(1000.);
		theCSMatForScatProjToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
   	}
	E1=E2;
   	E2*=fE;
   }
   
   
  std::vector<G4AdjointCSMatrix*> res;
  res.clear();
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering);
  

/*
  G4String file_name;
  std::stringstream astream;
  G4String str_Z;
  astream<<Z;
  astream>>str_Z;  
  theCSMatForProdToProjBackwardScattering->Write(aModel->GetName()+G4String("_CSMat_Z")+str_Z+"_ProdToProj.txt"); 
  theCSMatForScatProjToProjBackwardScattering->Write(aModel->GetName()+G4String("_CSMat_Z")+str_Z+"_ScatProjToProj.txt");
 
*/

  
  return res;
  
  
}
///////////////////////////////////////////////////////
//
std::vector<G4AdjointCSMatrix*>
G4AdjointCSManager::BuildCrossSectionsMatricesForAGivenModelAndMaterial(G4VEmAdjointModel* aModel,
								        G4Material* aMaterial,
								        G4int nbin_pro_decade)
{ 
  G4AdjointCSMatrix* theCSMatForProdToProjBackwardScattering = new G4AdjointCSMatrix(false);
  G4AdjointCSMatrix* theCSMatForScatProjToProjBackwardScattering = new G4AdjointCSMatrix(true);
 
 
  //make the vector of primary energy of the adjoint particle, could try to make this just once ? 
  
   G4double EkinMin =aModel->GetLowEnergyLimit();
   G4double EkinMaxForScat =aModel->GetHighEnergyLimit()*0.999;
   G4double EkinMaxForProd =aModel->GetHighEnergyLimit()*0.999;
   if (aModel->GetSecondPartOfSameType() )EkinMaxForProd =EkinMaxForProd/2.;
   
    
   
   
   
   
   
   //Product to projectile backward scattering
   //-----------------------------------------
   G4double fE=std::pow(10.,1./nbin_pro_decade);
   G4double E2=std::pow(10.,double( int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   G4double E1=EkinMin;
   // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
   while (E1 <EkinMaxForProd){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForProd,E1);
	std::vector< std::vector< double>* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerVolumeForSecond(aMaterial,E1,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< double>* log_ESecVec=aMat[0];
		std::vector< double>* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
	
		//normalise CSVec such that it becomes a probability vector
		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<G4endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
			//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<G4endl;
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-std::log(1000.);
		theCSMatForProdToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
	}	
	
 	
	
   	E1=E2;
   	E2*=fE;
   }
   
   //Scattered projectile to projectile backward scattering
   //-----------------------------------------
   
   E2=std::pow(10.,double( int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   E1=EkinMin;
   // Loop checking, 07-Aug-2015, Vladimir Ivanchenko
   while (E1 <EkinMaxForScat){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForScat,E1);
	std::vector< std::vector< double>* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerVolumeForScatProj(aMaterial,E1,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< double>* log_ESecVec=aMat[0];
		std::vector< double>* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
	
		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<G4endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
		//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<G4endl;if (theAdjPartDef->GetParticleName() == "adj_gamma") return G4Gamma::Gamma();
 
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-std::log(1000.);
	
		theCSMatForScatProjToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
   	}
	E1=E2;
   	E2*=fE;	
   }
   
   
   
   
   
   
   
  std::vector<G4AdjointCSMatrix*> res;
  res.clear();
  
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering); 
  
 /*
  theCSMatForProdToProjBackwardScattering->Write(aModel->GetName()+"_CSMat_"+aMaterial->GetName()+"_ProdToProj.txt");
  theCSMatForScatProjToProjBackwardScattering->Write(aModel->GetName()+"_CSMat_"+aMaterial->GetName()+"_ScatProjToProj.txt");
*/


  return res;
  
  
}

///////////////////////////////////////////////////////
//
G4ParticleDefinition* G4AdjointCSManager::GetAdjointParticleEquivalent(G4ParticleDefinition* theFwdPartDef)
{
 if (theFwdPartDef->GetParticleName() == "e-") return G4AdjointElectron::AdjointElectron();
 else if (theFwdPartDef->GetParticleName() == "gamma") return G4AdjointGamma::AdjointGamma();
 else if (theFwdPartDef->GetParticleName() == "proton") return G4AdjointProton::AdjointProton();
 else if (theFwdPartDef ==theFwdIon) return theAdjIon;
 
 return 0; 	
}
///////////////////////////////////////////////////////
//
G4ParticleDefinition* G4AdjointCSManager::GetForwardParticleEquivalent(G4ParticleDefinition* theAdjPartDef)
{
 if (theAdjPartDef->GetParticleName() == "adj_e-") return G4Electron::Electron();
 else if (theAdjPartDef->GetParticleName() == "adj_gamma") return G4Gamma::Gamma();
 else if (theAdjPartDef->GetParticleName() == "adj_proton") return G4Proton::Proton();
 else if (theAdjPartDef == theAdjIon) return theFwdIon;
 return 0; 	
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::DefineCurrentMaterial(const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = const_cast<G4MaterialCutsCouple*> (couple);
    currentMaterial = const_cast<G4Material*> (couple->GetMaterial());
    currentMatIndex = couple->GetIndex();
    lastPartDefForCS =0;
    LastEkinForCS =0;
    LastCSCorrectionFactor =1.;
  }  
}

///////////////////////////////////////////////////////
//
void G4AdjointCSManager::DefineCurrentParticle(const G4ParticleDefinition* aPartDef)
{
  if(aPartDef != currentParticleDef) {
  
  	currentParticleDef= const_cast< G4ParticleDefinition* > (aPartDef);
	massRatio=1;
	if (aPartDef == theAdjIon) massRatio = proton_mass_c2/aPartDef->GetPDGMass();
  	currentParticleIndex=1000000;
 	for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  		if (aPartDef == theListOfAdjointParticlesInAction[i]) currentParticleIndex=i;
  	}	
	 
  }  
}



/////////////////////////////////////////////////////////////////////////////////////////////////
//
G4double G4AdjointCSManager::ComputeAdjointCS(G4double aPrimEnergy,G4AdjointCSMatrix*
				anAdjointCSMatrix,G4double Tcut)
{ 
  std::vector< double> *theLogPrimEnergyVector = anAdjointCSMatrix->GetLogPrimEnergyVector();
  if (theLogPrimEnergyVector->size() ==0){
 	G4cout<<"No data are contained in the given AdjointCSMatrix!"<<G4endl;
	G4cout<<"The s"<<G4endl;
	return 0.;
	
  }
  G4double log_Tcut = std::log(Tcut);
  G4double log_E =std::log(aPrimEnergy);
  
  if (aPrimEnergy <= Tcut || log_E > theLogPrimEnergyVector->back()) return 0.;
  
  

  G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
 
  size_t ind =theInterpolator->FindPositionForLogVector(log_E,*theLogPrimEnergyVector);
  G4double aLogPrimEnergy1,aLogPrimEnergy2;
  G4double aLogCS1,aLogCS2;
  G4double log01,log02;
  std::vector< double>* aLogSecondEnergyVector1 =0;
  std::vector< double>* aLogSecondEnergyVector2  =0;
  std::vector< double>* aLogProbVector1=0;
  std::vector< double>* aLogProbVector2=0;
  std::vector< size_t>* aLogProbVectorIndex1=0;
  std::vector< size_t>* aLogProbVectorIndex2=0;
  
	 							     
  anAdjointCSMatrix->GetData(ind, aLogPrimEnergy1,aLogCS1,log01, aLogSecondEnergyVector1,aLogProbVector1,aLogProbVectorIndex1);
  anAdjointCSMatrix->GetData(ind+1, aLogPrimEnergy2,aLogCS2,log02, aLogSecondEnergyVector2,aLogProbVector2,aLogProbVectorIndex2);
  if (anAdjointCSMatrix->IsScatProjToProjCase()){ //case where the Tcut plays a role
	G4double log_minimum_prob1, log_minimum_prob2;
	log_minimum_prob1=theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector1,*aLogProbVector1);
	log_minimum_prob2=theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector2,*aLogProbVector2);
	aLogCS1+= log_minimum_prob1;
	aLogCS2+= log_minimum_prob2;
  }
 
  G4double log_adjointCS = theInterpolator->LinearInterpolation(log_E,aLogPrimEnergy1,aLogPrimEnergy2,aLogCS1,aLogCS2);
  return std::exp(log_adjointCS); 
  
  
}	 
