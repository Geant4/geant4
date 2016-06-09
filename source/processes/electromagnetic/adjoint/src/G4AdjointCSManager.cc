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
#include "G4AdjointCSManager.hh"
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
#include "G4AdjointElectron.hh"
#include "G4AdjointGamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProductionCutsTable.hh"


G4AdjointCSManager* G4AdjointCSManager::theInstance = 0;
///////////////////////////////////////////////////////
//
G4AdjointCSManager* G4AdjointCSManager::GetAdjointCSManager()
{ if(theInstance == 0) {
    static G4AdjointCSManager ins;
     theInstance = &ins;
  }
 return theInstance; 
}

///////////////////////////////////////////////////////
//
G4AdjointCSManager::G4AdjointCSManager()
{ CrossSectionMatrixesAreBuilt=false;
  theTotalForwardSigmaTableVector.clear();
  theTotalAdjointSigmaTableVector.clear();
  listOfForwardEmProcess.clear();
  listOfForwardEnergyLossProcess.clear();
  theListOfAdjointParticlesInAction.clear();
  Tmin=0.1*keV;
  Tmax=100.*TeV;
  nbins=240;
  
  RegisterAdjointParticle(G4AdjointElectron::AdjointElectron());
  RegisterAdjointParticle(G4AdjointGamma::AdjointGamma());
  
  verbose  = 1;
  
  consider_continuous_weight_correction =true;
  consider_poststep_weight_correction =false;
 
}
///////////////////////////////////////////////////////
//
G4AdjointCSManager::~G4AdjointCSManager()
{;
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterEmAdjointModel(G4VEmAdjointModel* aModel)
{listOfAdjointEMModel.push_back(aModel);
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterEmProcess(G4VEmProcess* aProcess, G4ParticleDefinition* aFwdPartDef)
{ 
  G4ParticleDefinition* anAdjPartDef = GetAdjointParticleEquivalent(aFwdPartDef);
  if (anAdjPartDef && aProcess){
  	RegisterAdjointParticle(anAdjPartDef);
	int index=-1;
	
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
  	int index=-1;
  	for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  		if (anAdjPartDef->GetParticleName() == theListOfAdjointParticlesInAction[i]->GetParticleName()) index=i;
  	}
  	listOfForwardEnergyLossProcess[index]->push_back(aProcess);
   }
}
///////////////////////////////////////////////////////
//
void G4AdjointCSManager::RegisterAdjointParticle(G4ParticleDefinition* aPartDef)
{  int index=-1;
   for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	if (aPartDef->GetParticleName() == theListOfAdjointParticlesInAction[i]->GetParticleName()) index=i;
   }
  
   if (index ==-1){
  	listOfForwardEnergyLossProcess.push_back(new std::vector<G4VEnergyLossProcess*>());
	theTotalForwardSigmaTableVector.push_back(new G4PhysicsTable);
	theTotalAdjointSigmaTableVector.push_back(new G4PhysicsTable);
	listOfForwardEmProcess.push_back(new std::vector<G4VEmProcess*>());
	theListOfAdjointParticlesInAction.push_back(aPartDef);
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
	for (size_t i=0; i<listOfAdjointEMModel.size();i++){
		G4VEmAdjointModel* aModel =listOfAdjointEMModel[i];
		G4cout<<"Build adjoint cross section matrices for "<<aModel->GetName()<<std::endl;
		if (aModel->GetUseMatrix()){
			std::vector<G4AdjointCSMatrix*>* aListOfMat1 = new std::vector<G4AdjointCSMatrix*>();
			std::vector<G4AdjointCSMatrix*>* aListOfMat2 = new std::vector<G4AdjointCSMatrix*>();
			aListOfMat1->clear();
			aListOfMat2->clear();
			if (aModel->GetUseMatrixPerElement()){
				if (aModel->GetUseOnlyOneMatrixForAllElements()){
						std::vector<G4AdjointCSMatrix*>
						two_matrices=BuildCrossSectionsMatricesForAGivenModelAndElement(aModel,1, 1, 10);
						aListOfMat1->push_back(two_matrices[0]);
						aListOfMat2->push_back(two_matrices[1]);
				}
				else {		
					for (size_t j=0; j<theElementTable->size();j++){
						G4Element* anElement=(*theElementTable)[j];
						G4int Z = G4int(anElement->GetZ());
						G4int A = G4int(anElement->GetA());
						std::vector<G4AdjointCSMatrix*>
							two_matrices=BuildCrossSectionsMatricesForAGivenModelAndElement(aModel,Z, A, 10);
						aListOfMat1->push_back(two_matrices[0]);
						aListOfMat2->push_back(two_matrices[1]);
					}
				}	
			}
			else { //Per material case
				for (size_t j=0; j<theMaterialTable->size();j++){
					G4Material* aMaterial=(*theMaterialTable)[j];
					std::vector<G4AdjointCSMatrix*>
						two_matrices=BuildCrossSectionsMatricesForAGivenModelAndMaterial(aModel,aMaterial, 10);
					aListOfMat1->push_back(two_matrices[0]);
					aListOfMat2->push_back(two_matrices[1]);
				}
			
			}
			theAdjointCSMatricesForProdToProj.push_back(*aListOfMat1);
			theAdjointCSMatricesForScatProjToProj.push_back(*aListOfMat2);	
			aModel->SetCSMatrices(aListOfMat1, aListOfMat2); 	
		}
		else {  std::vector<G4AdjointCSMatrix*> two_empty_matrices;
			theAdjointCSMatricesForProdToProj.push_back(two_empty_matrices);
			theAdjointCSMatricesForScatProjToProj.push_back(two_empty_matrices);
			
		}		
	}
	G4cout<<"All adjoint cross section matrices are built "<<std::endl;
	CrossSectionMatrixesAreBuilt = true;
}


///////////////////////////////////////////////////////
//
void G4AdjointCSManager::BuildTotalSigmaTables()
{ 
  const G4ProductionCutsTable* theCoupleTable= G4ProductionCutsTable::GetProductionCutsTable();
  for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	G4ParticleDefinition* thePartDef = theListOfAdjointParticlesInAction[i];
  	theTotalForwardSigmaTableVector[i]->clearAndDestroy();
	theTotalAdjointSigmaTableVector[i]->clearAndDestroy();
  	for (size_t j=0;j<theCoupleTable->GetTableSize();j++){
		const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
		
		//make first the total fwd CS table for FwdProcess
		G4PhysicsVector* aVector =  new G4PhysicsLogVector(Tmin, Tmax, nbins);
		for(size_t l=0; l<aVector->GetVectorLength(); l++) { 
			G4double totCS=0;
			G4double e=aVector->GetLowEdgeEnergy(l);
			for (size_t k=0; k<listOfForwardEmProcess[i]->size(); k++){
				totCS+=(*listOfForwardEmProcess[i])[k]->GetLambda(e, couple);
			}
			for (size_t k=0; k<listOfForwardEnergyLossProcess[i]->size(); k++){
				totCS+=(*listOfForwardEnergyLossProcess[i])[k]->GetLambda(e, couple);
			}
			//G4cout<<totCS<<std::endl;
			aVector->PutValue(l,totCS);
		
		}
		theTotalForwardSigmaTableVector[i]->push_back(aVector);
		
		G4PhysicsVector* aVector1 =  new G4PhysicsLogVector(Tmin, Tmax, nbins);
		for(size_t l=0; l<aVector->GetVectorLength(); l++) { 
			G4double e=aVector->GetLowEdgeEnergy(l);
			G4double totCS =ComputeTotalAdjointCS(couple,thePartDef,e);
			//G4cout<<totCS<<std::endl;
			aVector1->PutValue(l,totCS);
		
		}	
		theTotalAdjointSigmaTableVector[i]->push_back(aVector1);
		
	}
  }
   
}
///////////////////////////////////////////////////////
//
G4double G4AdjointCSManager::GetTotalAdjointCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple)
{ DefineCurrentMaterial(aCouple);
  int index=-1;
  for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	if (aPartDef == theListOfAdjointParticlesInAction[i]) index=i;
  }	
  if (index == -1) return 0.;
  
  G4bool b;
  return (((*theTotalAdjointSigmaTableVector[index])[currentMatIndex])->GetValue(Ekin, b));
  
  
  
}				     
///////////////////////////////////////////////////////
//				     
G4double G4AdjointCSManager::GetTotalForwardCS(G4ParticleDefinition* aPartDef, G4double Ekin,
	 		   	     const G4MaterialCutsCouple* aCouple)
{ DefineCurrentMaterial(aCouple);
  int index=-1;
  for (size_t i=0;i<theListOfAdjointParticlesInAction.size();i++){
  	if (aPartDef == theListOfAdjointParticlesInAction[i]) index=i;
  }	
  if (index == -1) return 0.;
  G4bool b;
  return (((*theTotalForwardSigmaTableVector[index])[currentMatIndex])->GetValue(Ekin, b));
  
  
}				     
///////////////////////////////////////////////////////
//				     			     
G4double G4AdjointCSManager::GetContinuousWeightCorrection(G4ParticleDefinition* aPartDef, G4double PreStepEkin,G4double AfterStepEkin,
	 		   	     const G4MaterialCutsCouple* aCouple, G4double step_length)
{ //G4double fwdCS = GetTotalForwardCS(aPartDef, AfterStepEkin,aCouple);
  
  G4double corr_fac = 1.;
  if (consider_continuous_weight_correction) {
  	
	G4double adjCS = GetTotalAdjointCS(aPartDef, PreStepEkin,aCouple);
	G4double PrefwdCS;
	PrefwdCS = GetTotalForwardCS(aPartDef, PreStepEkin,aCouple);
  	G4double fwdCS = GetTotalForwardCS(aPartDef, (AfterStepEkin+PreStepEkin)/2.,aCouple);
	G4cout<<adjCS<<'\t'<<fwdCS<<std::endl;
	//if (aPartDef ==G4AdjointGamma::AdjointGamma()) G4cout<<adjCS<<'\t'<<fwdCS<<std::endl;
  	/*if (adjCS >0 ) corr_fac = std::exp((PrefwdCS-fwdCS)*step_length);
	else corr_fac = std::exp(-fwdCS*step_length);*/
	corr_fac *=std::exp((adjCS-fwdCS)*step_length);
	corr_fac=std::max(corr_fac,1.e-6);
	corr_fac *=PreStepEkin/AfterStepEkin;
   	
  }
  G4cout<<"Cont "<<corr_fac<<std::endl;
  G4cout<<"Ekin0 "<<PreStepEkin<<std::endl;
  G4cout<<"Ekin1 "<<AfterStepEkin<<std::endl;
  G4cout<<"step_length "<<step_length<<std::endl;
  return corr_fac; 
}				     
///////////////////////////////////////////////////////
//				     			     
G4double G4AdjointCSManager::GetPostStepWeightCorrection(G4ParticleDefinition* , G4ParticleDefinition* ,
					    		  G4double ,G4double ,
	 		   	     			  const G4MaterialCutsCouple* )
{ G4double corr_fac = 1.;
  if (consider_poststep_weight_correction) {
	/*G4double fwdCS = GetTotalForwardCS(aSecondPartDef, EkinPrim,aCouple);
  	G4double adjCS = GetTotalAdjointCS(aPrimPartDef, EkinPrim,aCouple);*/
  	//G4double fwd1CS = GetTotalForwardCS(aPrimPartDef, EkinPrim,aCouple);
  	//if (adjCS>0 && fwd1CS>0) adjCS = fwd1CS;
  	//corr_fac =fwdCS*EkinSecond/adjCS/EkinPrim;
  	//corr_fac = adjCS/fwdCS;
  }
  return corr_fac;
}							  
///////////////////////////////////////////////////////
//
double  G4AdjointCSManager::ComputeAdjointCS(G4Material* aMaterial,
					    		 G4VEmAdjointModel* aModel, 
					    		 G4double PrimEnergy,
					    		 G4double Tcut,
					    		 G4bool IsScatProjToProjCase,
							 std::vector<double>& CS_Vs_Element)
{ 
  
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
		return aModel->AdjointCrossSection(currentCouple,PrimEnergy,IsScatProjToProjCase);
				         
	
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
			for (size_t i=0;i<n_el;i++){
				size_t ind_el = aMaterial->GetElement(i)->GetIndex();
				factor+=aMaterial->GetElement(i)->GetZ()*aMaterial->GetVecNbOfAtomsPerVolume()[i];
				G4AdjointCSMatrix* theCSMatrix;
				if (IsScatProjToProjCase){
					theCSMatrix=theAdjointCSMatricesForScatProjToProj[ind_model][ind_el];
				}
				else  	theCSMatrix=theAdjointCSMatricesForProdToProj[ind_model][ind_el];
				//G4double CS =0.;
				
				//G4cout<<CS<<std::endl;			
				
			}
			CS *=factor;
			CS_Vs_Element.push_back(CS);
								
		}
		else {
			for (size_t i=0;i<n_el;i++){
				size_t ind_el = aMaterial->GetElement(i)->GetIndex();
				//G4cout<<aMaterial->GetName()<<std::endl;
				G4AdjointCSMatrix* theCSMatrix;
				if (IsScatProjToProjCase){
					theCSMatrix=theAdjointCSMatricesForScatProjToProj[ind_model][ind_el];
				}
				else  	theCSMatrix=theAdjointCSMatricesForProdToProj[ind_model][ind_el];
				G4double CS =0.;
				if (PrimEnergy > Tlow)
					CS = ComputeAdjointCS(PrimEnergy,theCSMatrix,Tlow);
				//G4cout<<CS<<std::endl;			
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
  	CS+=CS_Vs_Element[i];
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
{ std::vector<double> CS_Vs_Element;
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
// G4ParticleDefinition* theDirPartDef = GetForwardParticleEquivalent(aPartDef);
 DefineCurrentMaterial(aCouple);
/* size_t idx=-1;
 if (theDirPartDef->GetParticleName() == "gamma") idx = 0;
 else if (theDirPartDef->GetParticleName() == "e-") idx = 1;
 else if (theDirPartDef->GetParticleName() == "e+") idx = 2;
 
 //THe tCut computation is wrong this should be on Tcut per model the secondary determioming the Tcut
 const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
 //G4cout<<aVec<<std::endl;
 G4double Tcut =(*aVec)[aCouple->GetIndex()];*/
 //G4cout<<"Tcut "<<Tcut<<std::endl;
 //G4cout<<(*aVec)[0]<<std::endl;
// G4double Tcut =converters[idx]->Convert(Rcut,aCouple->GetMaterial());
 
  
 std::vector<double> CS_Vs_Element;
 for (size_t i=0; i<listOfAdjointEMModel.size();i++){
 	/*G4ParticleDefinition* theDirSecondPartDef = 
			GetForwardParticleEquivalent(listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectSecondaryParticleDefinition());
 	
	*/
	
	
	G4double Tlow=0;
	if (!listOfAdjointEMModel[i]->GetApplyCutInRange()) Tlow =listOfAdjointEMModel[i]->GetLowEnergyLimit();
	else {
		G4ParticleDefinition* theDirSecondPartDef = 
			GetForwardParticleEquivalent(listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectSecondaryParticleDefinition());
		G4int idx=-1;
 		if (theDirSecondPartDef->GetParticleName() == "gamma") idx = 0;
		else if (theDirSecondPartDef->GetParticleName() == "e-") idx = 1;
 		else if (theDirSecondPartDef->GetParticleName() == "e+") idx = 2;
		const std::vector<G4double>* aVec = G4ProductionCutsTable::GetProductionCutsTable()->GetEnergyCutsVector(idx);
 		Tlow =(*aVec)[aCouple->GetIndex()];
		
	
	}
 	
 	if ( Ekin<=listOfAdjointEMModel[i]->GetHighEnergyLimit() && Ekin>=listOfAdjointEMModel[i]->GetLowEnergyLimit()){
		if (aPartDef == listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectPrimaryParticleDefinition()){
			//G4cout<<"Yes1 before "<<std::endl;
			TotalCS += ComputeAdjointCS(currentMaterial,
					    	       listOfAdjointEMModel[i], 
					    	       Ekin, Tlow,true,CS_Vs_Element);
			//G4cout<<"Yes1 "<<Ekin<<'\t'<<TotalCS<<std::endl;			       
		}
		if (aPartDef == listOfAdjointEMModel[i]->GetAdjointEquivalentOfDirectSecondaryParticleDefinition()){
			TotalCS += ComputeAdjointCS(currentMaterial,
					    	       listOfAdjointEMModel[i], 
					    	       Ekin, Tlow,false, CS_Vs_Element);
						       
			//G4cout<<"Yes2 "<<TotalCS<<std::endl;
		}
		
	}
 }
 return TotalCS;
    
 
}	
///////////////////////////////////////////////////////
//
std::vector<G4AdjointCSMatrix*>
G4AdjointCSManager::BuildCrossSectionsMatricesForAGivenModelAndElement(G4VEmAdjointModel* aModel,G4int Z,G4int A,
								        int nbin_pro_decade)
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
   G4double E2=std::pow(10.,G4double( G4int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   G4double E1=EkinMin;
   while (E1 <EkinMaxForProd){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForProd,E1);
	std::vector< std::vector< G4double >* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerAtomForSecond(E1,Z,A,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< G4double >* log_ESecVec=aMat[0];
		std::vector< G4double >* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
		//normalise CSVec such that it becomes a probability vector
 		/*for (size_t j=0;j<log_CSVec->size();j++) (*log_CSVec)[j]=(*log_CSVec)[j]-log_adjointCS;
		(*log_CSVec)[0]=-90.;*/
	
	
		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<std::endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
			//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<std::endl;
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-1.;
		theCSMatForProdToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
	}	
   	E1=E2;
   	E2*=fE;
   }
   
   //Scattered projectile to projectile backward scattering
   //-----------------------------------------
   
   E2=std::pow(10.,G4double( G4int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   E1=EkinMin;
   while (E1 <EkinMaxForScat){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForScat,E1);
	std::vector< std::vector< G4double >* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerAtomForScatProj(E1,Z,A,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< G4double >* log_ESecVec=aMat[0];
		std::vector< G4double >* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
		//normalise CSVec such that it becomes a probability vector
 		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<std::endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
		//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<std::endl;
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-1.;
		theCSMatForScatProjToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
   	}
	E1=E2;
   	E2*=fE;
   }
   
   
   
   
   
   
   
  std::vector<G4AdjointCSMatrix*> res;
  res.clear();
  
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering);
  

#ifdef TEST_MODE
  G4String file_name;
  std::stringstream astream;
  G4String str_Z;
  astream<<Z;
  astream>>str_Z;  
  theCSMatForProdToProjBackwardScattering->Write(aModel->GetName()+G4String("_CSMat_Z")+str_Z+"_ProdToProj.txt"); 
  theCSMatForScatProjToProjBackwardScattering->Write(aModel->GetName()+G4String("_CSMat_Z")+str_Z+"_ScatProjToProj.txt");
 
  /*G4AdjointCSMatrix* aMat1 = new G4AdjointCSMatrix(false);
  G4AdjointCSMatrix* aMat2 = new G4AdjointCSMatrix(true);
  
  aMat1->Read(G4String("test_Z")+str_Z+"_1.txt");
  aMat2->Read(G4String("test_Z")+str_Z+"_2.txt");
  aMat1->Write(G4String("test_Z")+str_Z+"_11.txt");
  aMat2->Write(G4String("test_Z")+str_Z+"_22.txt"); */
#endif 
  
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
   G4double E2=std::pow(10.,G4double( G4int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   G4double E1=EkinMin;
   while (E1 <EkinMaxForProd){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForProd,E1);
	std::vector< std::vector< G4double >* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerVolumeForSecond(aMaterial,E1,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< G4double >* log_ESecVec=aMat[0];
		std::vector< G4double >* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
	
		//normalise CSVec such that it becomes a probability vector
		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<std::endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
			//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<std::endl;
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-1.;
		theCSMatForProdToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
	}	
	
 	
	
   	E1=E2;
   	E2*=fE;
   }
   
   //Scattered projectile to projectile backward scattering
   //-----------------------------------------
   
   E2=std::pow(10.,G4double( G4int(std::log10(EkinMin)*nbin_pro_decade)+1)/nbin_pro_decade)/fE;
   E1=EkinMin;
   while (E1 <EkinMaxForScat){
   	E1=std::max(EkinMin,E2);
	E1=std::min(EkinMaxForScat,E1);
	std::vector< std::vector< G4double >* >  aMat= aModel->ComputeAdjointCrossSectionVectorPerVolumeForScatProj(aMaterial,E1,nbin_pro_decade);
	if (aMat.size()>=2) {
		std::vector< G4double >* log_ESecVec=aMat[0];
		std::vector< G4double >* log_CSVec=aMat[1];
		G4double log_adjointCS=log_CSVec->back();
	
		for (size_t j=0;j<log_CSVec->size();j++) {
	        	//G4cout<<"CSMan1 "<<(*log_CSVec)[j]<<std::endl;
			if (j==0) (*log_CSVec)[j] = 0.; 
			else (*log_CSVec)[j]=std::log(1.-std::exp((*log_CSVec)[j]-log_adjointCS));
		//G4cout<<"CSMan2 "<<(*log_CSVec)[j]<<std::endl;
		}	
		(*log_CSVec)[log_CSVec->size()-1]=(*log_CSVec)[log_CSVec->size()-2]-1.;
	
		theCSMatForScatProjToProjBackwardScattering->AddData(std::log(E1),log_adjointCS,log_ESecVec,log_CSVec,0);
   	}
	E1=E2;
   	E2*=fE;	
   }
   
   
   
   
   
   
   
  std::vector<G4AdjointCSMatrix*> res;
  res.clear();
  
  res.push_back(theCSMatForProdToProjBackwardScattering);
  res.push_back(theCSMatForScatProjToProjBackwardScattering); 
  
#ifdef TEST_MODE
  theCSMatForProdToProjBackwardScattering->Write(aModel->GetName()+"_CSMat_"+aMaterial->GetName()+"_ProdToProj.txt");
  theCSMatForScatProjToProjBackwardScattering->Write(aModel->GetName()+"_CSMat_"+aMaterial->GetName()+"_ScatProjToProj.txt");
#endif


  return res;
  
  
}

///////////////////////////////////////////////////////
//
G4ParticleDefinition* G4AdjointCSManager::GetAdjointParticleEquivalent(G4ParticleDefinition* theFwdPartDef)
{
 if (theFwdPartDef->GetParticleName() == "e-") return G4AdjointElectron::AdjointElectron();
 if (theFwdPartDef->GetParticleName() == "gamma") return G4AdjointGamma::AdjointGamma();
 return 0; 	
}
///////////////////////////////////////////////////////
//
G4ParticleDefinition* G4AdjointCSManager::GetForwardParticleEquivalent(G4ParticleDefinition* theAdjPartDef)
{
 if (theAdjPartDef->GetParticleName() == "adj_e-") return G4Electron::Electron();
 if (theAdjPartDef->GetParticleName() == "adj_gamma") return G4Gamma::Gamma();
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
    //G4cout<<"Index material "<<currentMatIndex<<std::endl;
  }  
}



///////////////////////////////////////////////////////
//
double G4AdjointCSManager::ComputeAdjointCS(G4double aPrimEnergy,G4AdjointCSMatrix*
				anAdjointCSMatrix,G4double Tcut)
{ 
  std::vector< G4double > *theLogPrimEnergyVector = anAdjointCSMatrix->GetLogPrimEnergyVector();
  if (theLogPrimEnergyVector->size() ==0){
 	G4cout<<"No data are contained in the given AdjointCSMatrix!"<<std::endl;
	G4cout<<"The sampling procedure will be stopped."<<std::endl;
	return 0.;
	
  }
  //G4cout<<"A prim/Tcut "<<aPrimEnergy<<'\t'<<Tcut<<std::endl; 
  G4double log_Tcut = std::log(Tcut);
  G4double log_E =std::log(aPrimEnergy);
  
  if (aPrimEnergy <= Tcut || log_E > theLogPrimEnergyVector->back()) return 0.;
  
  

  G4AdjointInterpolator* theInterpolator=G4AdjointInterpolator::GetInstance();
 
  size_t ind =theInterpolator->FindPositionForLogVector(log_E,*theLogPrimEnergyVector);
  //G4cout<<"Prim energy "<<(*thePrimEnergyVector)[0]<<std::endl;
  //G4cout<<"Prim energy[ind]"<<(*thePrimEnergyVector)[ind]<<std::endl;
  //G4cout<<"Prim energy ind"<<ind<<std::endl;
  
  G4double aLogPrimEnergy1,aLogPrimEnergy2;
  G4double aLogCS1,aLogCS2;
  G4double log01,log02;
  std::vector< G4double>* aLogSecondEnergyVector1 =0;
  std::vector< G4double>* aLogSecondEnergyVector2  =0;
  std::vector< G4double>* aLogProbVector1=0;
  std::vector< G4double>* aLogProbVector2=0;
  std::vector< size_t>* aLogProbVectorIndex1=0;
  std::vector< size_t>* aLogProbVectorIndex2=0;
  
	 							     
  anAdjointCSMatrix->GetData(ind, aLogPrimEnergy1,aLogCS1,log01, aLogSecondEnergyVector1,aLogProbVector1,aLogProbVectorIndex1);
  anAdjointCSMatrix->GetData(ind+1, aLogPrimEnergy2,aLogCS2,log02, aLogSecondEnergyVector2,aLogProbVector2,aLogProbVectorIndex2);
  //G4cout<<"aSecondEnergyVector1.size() "<<aSecondEnergyVector1->size()<<std::endl;
  //G4cout<<aSecondEnergyVector1<<std::endl;
  //G4cout<<"aSecondEnergyVector2.size() "<<aSecondEnergyVector2->size()<<std::endl;
  if (anAdjointCSMatrix->IsScatProjToProjCase()){ //case where the Tcut plays a role
	G4double log_minimum_prob1, log_minimum_prob2;
	
	//G4cout<<aSecondEnergyVector1->size()<<std::endl;
	log_minimum_prob1=theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector1,*aLogProbVector1);
	log_minimum_prob2=theInterpolator->InterpolateForLogVector(log_Tcut,*aLogSecondEnergyVector2,*aLogProbVector2);
	//G4cout<<"minimum_prob1 "<< std::exp(log_minimum_prob1)<<std::endl;
	//G4cout<<"minimum_prob2 "<< std::exp(log_minimum_prob2)<<std::endl;
	//G4cout<<"Tcut "<<std::endl;
	aLogCS1+= log_minimum_prob1;
	aLogCS2+= log_minimum_prob2;
  }
 
  G4double log_adjointCS = theInterpolator->LinearInterpolation(log_E,aLogPrimEnergy1,aLogPrimEnergy2,aLogCS1,aLogCS2);
  return std::exp(log_adjointCS); 
  
  
}	 
