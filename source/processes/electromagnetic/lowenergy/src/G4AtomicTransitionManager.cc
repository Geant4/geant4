//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AtomicTransitionManager.cc,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Authors: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//          Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// ?????      Created
//
// -------------------------------------------------------------------

#include "G4AtomicTransitionManager.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "G4ShellData.hh"
#include "G4FluoData.hh"

// constructor

G4AtomicTransitionManager::G4AtomicTransitionManager(G4int minZ, G4int maxZ, G4int limitInfTable,G4int limitSupTable)
  :zMin(minZ), zMax(maxZ),infTableLimit(limitInfTable),supTableLimit(limitSupTable)
{
 
  G4ShellData* shellManager = new G4ShellData;
  
  shellManager->LoadData("binding");
 
  for (G4int Z = zMin; Z<= zMax; Z++)
    {
      G4std::vector<G4AtomicShell*> vectorOfShells;  
    
      size_t numberOfShells=shellManager->NumberOfShells(Z);
      for (G4int shellIndex = 0; shellIndex<numberOfShells; shellIndex++)
	{ 
	  G4int shellId = shellManager->ShellId(Z,shellIndex);
	  G4double bindingEnergy = shellManager->BindingEnergy(Z,shellIndex);
	 
	  G4AtomicShell * shell = new G4AtomicShell(shellId,bindingEnergy);
	
	  vectorOfShells.push_back(shell);
	}
    
      shellTable.insert(std::make_pair(Z, vectorOfShells));
    }
  for (G4int Znum= infTableLimit; Znum<=supTableLimit; Znum++)
    {  G4FluoData* fluoManager = new G4FluoData;
      G4std::vector<G4AtomicTransition*> vectorOfTransitions;
      fluoManager->LoadData(Znum);
     
      size_t numberOfVacancies = fluoManager-> NumberOfVacancies();
     
      for (size_t vacancyIndex = 0; vacancyIndex<numberOfVacancies;  vacancyIndex++)
	
	{
	  G4std::vector<G4int>  vectorOfIds;
	  G4DataVector vectorOfEnergies;
	  G4DataVector vectorOfProbabilities;
	
	  G4int finalShell = fluoManager->VacancyId(vacancyIndex);
	  size_t numberOfTransitions = fluoManager->NumberOfTransitions(vacancyIndex);
	  for (size_t origShellIndex = 0; origShellIndex <= numberOfTransitions;origShellIndex++)

	    {

	      G4int originatingShellId = fluoManager->StartShellId(origShellIndex,vacancyIndex);

	      vectorOfIds.push_back(originatingShellId);

	      G4double transitionEnergy = fluoManager->StartShellEnergy(origShellIndex,vacancyIndex);
	      vectorOfEnergies.push_back(transitionEnergy);
	      G4double transitionProbability = fluoManager->StartShellProb(origShellIndex,vacancyIndex);
	      vectorOfProbabilities.push_back(transitionProbability);
	    }
	  G4AtomicTransition * transition = new G4AtomicTransition (finalShell,vectorOfIds,
								    vectorOfEnergies,vectorOfProbabilities);
	  vectorOfTransitions.push_back(transition); 
	}
      transitionTable.insert(std::make_pair(Znum, vectorOfTransitions));
    
  delete fluoManager;
    }
 delete shellManager;
}
/*
  
  G4int ShellIndex=0;
  
  for (G4int index=0; index <(ZVector.size()); index++){  
  
  G4std::vector<G4AtomicShell*> vectorOfShells;  
  
  G4int Z = ZVector[index];
  
    G4FirstLevel* theBindEnVec = (*theBindingEnergyTable)[Z-1];
    
    G4int size = ((*theBindEnVec)[0])->size();
   
    for (ShellIndex=0; ShellIndex<size; ShellIndex++){
      
      G4int  id = (G4int) (*(*theBindEnVec)[0])[ShellIndex];
      
       G4DataVector finalId = 0;
    
       G4DataVector prob = 0;
    
       G4DataVector energies = 0;
      
      G4double energy = ((*(*theBindEnVec)[1])[ShellIndex])*MeV;
      
      G4SecondLevel* oneAtomTrans = (*theTransitionTable)[index];
      
      G4int FinalShellIndex = 0;
      
      G4int ShellCol = 0,   ProbCol = 1, EnergyCol = 2;
      
      G4int SizeOneAtomTrans = oneAtomTrans->size();
      
      for (FinalShellIndex=0; FinalShellIndex<SizeOneAtomTrans; FinalShellIndex++){
	
	G4int sizeOneFinalShellVec =(*(*oneAtomTrans)[FinalShellIndex])[ShellCol]->size();

	for(G4int k = 0; k<sizeOneFinalShellVec; k++){

	  G4int tmpId = (G4int) (*(*(*oneAtomTrans)[FinalShellIndex])[ShellCol])[0];

	  if (id == tmpId){

	  prob.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[ProbCol])[k]);

	  energies.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[EnergyCol])[k]);

	  finalId.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[ShellCol])[k]);

	  }
	  
	}
	
	
      }  
      
      //no radiative transition allowed  

      if (prob.size()==0)
	
	{prob.push_back(0.);
	
	energies.push_back(0.);
	
	finalId.push_back(0.);
	}
      
      G4AtomicShell *shell = new G4AtomicShell (id,energy,prob,energies,finalId);
      
      vectorOfShells.push_back(shell);
      
    }
    shellTable.insert(std::make_pair(Z, vectorOfShells));
    
  }
  
  delete theBindingEnergyTable;
  delete theTransitionTable;
  
}
*/
G4AtomicTransitionManager::~G4AtomicTransitionManager()

{ G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;
 
 for (pos = shellTable.begin(); pos != shellTable.end(); pos++){
   
   G4std::vector< G4AtomicShell*>vec = pos->second;
   
   G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
     
     delete vec[i]; 
   }
   
 }
 
 G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> >::iterator ppos;
 
 for (ppos = transitionTable.begin(); ppos != transitionTable.end(); ppos++){
   
   G4std::vector< G4AtomicTransition*>vec = ppos->second;
   
   G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
     
     delete vec[i]; 
   }
   
 }   
 //delete fluoManager;
 //delete shellManager;
}

 G4AtomicTransitionManager* G4AtomicTransitionManager::instance = 0;
 
G4AtomicTransitionManager* G4AtomicTransitionManager::Instance()
{
  if (instance==0)
    {
      instance = new G4AtomicTransitionManager;
     
    }
  return instance;
}

  //build the table for the binding energies.
//theBindingEnergyTable has three degrees of fredom (consequently three indexes):
//the first ranges from 0 to the size of ZVector and is used to choose the element.
//its value os the index of the position the Z of the choosen element occupies in ZVector
//the second index is used to discriminate between shells[0] and binding energies[0]
//the third index is used to choose a shell and its binding energy

/*
void G4AtomicTransitionManager::BuildZVec(){
  
  const G4MaterialTable* theMaterialTable=G4Material::GetMaterialTable();
  G4int numOfMaterials = theMaterialTable->length();
   
   for (G4int J=0 ; J < numOfMaterials; J++){ 
    
    const G4Material* material= (*theMaterialTable)[J];        
    const G4ElementVector* theElementVector = material->GetElementVector();
    const G4int NumberOfElements = material->GetNumberOfElements() ;
    
    for (G4int iel=0; iel<NumberOfElements; iel++ ){
      
      G4double Zel = (*theElementVector)(iel)->GetZ();
      
      G4int ZVectorSize = ZVector.size();
      
      G4int k = 0;
      
      for( G4int i = 0; i<ZVectorSize; i++){
	
	if (Zel !=ZVector[i]){k=k+1;}
	
	else{}
	
      }
      
      if ( k==ZVectorSize){
      
	ZVector.push_back(Zel);
	
      } 
      else
	{
	  continue;
	  
	}
    }
    
   }
   
}
*/
const G4AtomicShell* G4AtomicTransitionManager::Shell(G4int Z, G4int shellIndex)
{ 
  G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

  pos = shellTable.find(Z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = pos->second;

    if (shellIndex<v.size()){
      
      return(v[shellIndex]);
      
    }
    else {
      G4Exception("G4AtomicTransitionManager:shell not found");
      
    }
  }
  else{
    G4Exception("G4AtomicTransitionManager:Z not found");
    
  } 
}

const G4AtomicTransition* G4AtomicTransitionManager:: ReachableShell(G4int Z,G4int shellIndex)
{
 G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = pos->second;

    if (shellIndex<v.size()){
      
      return(v[shellIndex]);
      
    }
    else {
      G4Exception("G4AtomicTransitionManager:reachable shell not found");
      
    }
  }
  else{
    G4Exception("G4AtomicTransitionManager:Z not found");
    
  } 
}

G4int G4AtomicTransitionManager::NumberOfShells (G4int Z)
{

G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

  pos = shellTable.find(Z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = pos->second;

    return v.size();
  }

  else{
    G4Exception( "G4AtomicTransitionManager: Z not found" );
  } 
}

G4int G4AtomicTransitionManager::NumberOfReachableShells(G4int Z)
{
G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = pos->second;

    return v.size();
  }

  else{
    G4Exception( "G4AtomicTransitionManager: Z not found" );
  } 
}
G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(G4int Z, G4int shellIndex)

{

G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = pos->second;
  
    /* G4int index = v.size();

    for (G4int i=0; i<v.size(); i++){

      G4int tmpId = (v[i])->ShellId();
      
      if (tmpId==shellId)
	
	{index = i;}
      
      else{}
    }  
  
    if (index< v.size()){
  
    G4AtomicShell* shell=v[index];*/ 

    if (shellIndex<v.size()){

      G4AtomicTransition* transition=v[shellIndex];

   G4DataVector transProb = transition->TransitionProbabilities();

   G4double totalRadTransProb = 0;

   for(G4int j = 1; j<transProb.size(); j++){

     totalRadTransProb = totalRadTransProb + transProb[j];

   }

   return totalRadTransProb;   
      
    }
 else {
   G4Exception( "G4AtomicTransitionManager: shell not found" );
	
 }
  }
  else{
    G4Exception( "G4AtomicTransitionManager: Z not found");
  } 
}

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(G4int Z, G4int shellIndex)

{

  G4std::map<G4int,G4std::vector<G4AtomicTransition*>,std::less<G4int> >::iterator pos;
  
  pos = transitionTable.find(Z);
  
  if (pos!= transitionTable.end()){
    
    G4std::vector<G4AtomicTransition*> v = pos->second;
    
    /*G4int index = v.size();
    
    for (G4int i=0; i<v.size(); i++){
      
      G4int tmpId=(v[i])->ShellId();
      if (tmpId==shellId)
	{index = i;}
      else{}
    }  
    
    if (index< v.size()){
      
    G4AtomicShell* shell=v[index];*/ 
    
    if (shellIndex<v.size()){

      G4AtomicTransition* transition=v[shellIndex];
      
      G4DataVector transProb = transition->TransitionProbabilities();
      
      G4double totalRadTransProb = 0;
      
      for(G4int j = 1; j<transProb.size(); j++){

	totalRadTransProb = totalRadTransProb + transProb[j];
      }
      
      G4double totalNonRadTransProb= (1 - totalRadTransProb);
      
      return totalNonRadTransProb;    }
    
    else {
      G4Exception( "shell not found");
    }
  }
  else{
    G4Exception("Z not found");
  } 
}




