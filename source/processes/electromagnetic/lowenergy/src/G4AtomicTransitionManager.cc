// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4AtomicTransitionManager.cc
//
//      Authors:        Alfonso Mantero (alfonso.mantero@ge.infn.it)
//                      Elena Guardincerri (elena.guardincerri@ge.infn.it)
// 
//      Creation date: 4 May 2001
// -------------------------------------------------------------------

#include "G4AtomicTransitionManager.hh"
#include "G4LowEnergyUtilities.hh"
#include "G4EnergyLossTables.hh"
#include "G4Element.hh"
#include "G4ios.hh"
#include "g4std/fstream"

//----------..................-------------------....................

// constructor

G4AtomicTransitionManager::G4AtomicTransitionManager( )
  :ZVector(0),
   tableLimit(0)
{
  G4int tableLimit = 5;

  G4SecondLevel* theBindingEnergyTable=BuildBindingEnergiesTable();
  
  G4ThirdLevel* theTransitionTable= BuildTransitionTable();
  
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
	
	G4int sizeOneFinalShellVec = (*(*oneAtomTrans)[FinalShellIndex]).size();
	
	
	for(G4int k = 0; k<sizeOneFinalShellVec; k++){
	  
	  G4int tmpId = (G4int) (*(*(*oneAtomTrans)[FinalShellIndex])[ShellCol])[k];  
	  
	  if (id == tmpId){
	    
	    prob.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[ProbCol])[k]);
	    
	    energies.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[EnergyCol])[k]);
	    
	    finalId.push_back((*(*(*oneAtomTrans)[FinalShellIndex])[ShellCol])[0]);
	    
	  }
	  else{ 
	    
	  }
	  
	} 
	
      }  
      //no radiative transition allowed
     
      while (prob.size()<=ShellIndex)
	{
	  prob.push_back(0.);
	  
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

G4AtomicTransitionManager::~G4AtomicTransitionManager()

{ G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

 for (pos = shellTable.begin(); pos != shellTable.end(); pos++){
   
   G4std::vector< G4AtomicShell*>vec = pos->second;

     G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
     
     delete vec[i]; 
  }

 }
   
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

G4SecondLevel* G4AtomicTransitionManager::BuildBindingEnergiesTable()
{
  
  G4int dataNum = 2;
  G4SecondLevel*theBindingEnergyTable = util.BuildSecondLevelTables(0,dataNum,"fluor/binding");
  
return theBindingEnergyTable;
}
//this function builds the table for the radiative transitions
//they have four degrees of fredoom (four indexes):
//the first ranges from 0 to the size of ZVector and is used to choose the element.
//the second is used to select the values relative to all the radiative transition
//towards a given final shell
//the third is used to disctiminate between identities of the shells [0], transition 
//probabilities[1] and transition energies[2]
//the fourth is used to select a given starting shell for the transition
G4ThirdLevel* G4AtomicTransitionManager::BuildTransitionTable()
{
 
  BuildZVec();
 
   G4ThirdLevel* theTransitionTable = new G4ThirdLevel();
  
  G4int dataNumb = 3;
  
 for(G4int tableInd = 0; tableInd <ZVector.size()  ; tableInd++){

    G4int atomInd =  ZVector[tableInd];

    if(atomInd > tableLimit){
  
      G4SecondLevel* oneAtomTransProb = util.BuildSecondLevelTables(atomInd, dataNumb, "fluor/fl-tr-pr-");
     
      theTransitionTable->push_back(oneAtomTransProb);
    }
    else{
      ZVector.erase(remove( ZVector.begin(), ZVector.end(),atomInd),ZVector.end()) ;
     
       tableInd=tableInd-1;
    } 
   
  }

   return theTransitionTable;
}

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

const G4AtomicShell* G4AtomicTransitionManager::Shell(G4int z, G4int shellIdentifier)
{ 
  G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

  pos = shellTable.find(z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = pos->second;

G4int index =v.size();
 
    for (G4int i=0; i<v.size(); i++){

      G4int tmpId = (v[i])->ShellId();
      if (tmpId==shellIdentifier)
	{index = i;}
      else{}
    }  
  
 if (index< v.size()){
  
  
 return (v[index]);
      
    }
    else {
      G4Exception("shell not found");
	
    }
  }
  else{
    G4Exception("Z not found");
	
  } 
}

G4int G4AtomicTransitionManager::NumberOfShells (G4int z)

{

G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

  pos = shellTable.find(z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = pos->second;

    return v.size();
  }

  else{
    G4Exception( "Z not found" );
  } 
}

G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(G4int z, G4int shellId)

{

G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;

  pos = shellTable.find(z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = pos->second;
  
    G4int index = v.size();

    for (G4int i=0; i<v.size(); i++){

      G4int tmpId = (v[i])->ShellId();
      
      if (tmpId==shellId)
	
	{index = i;}
      
      else{}
    }  
  
 if (index< v.size()){
  
   G4AtomicShell* shell=v[index];

   G4DataVector transProb = shell->TransitionProbabilities();

   G4double totalRadTransProb = 0;

   for(G4int j = 0; j<transProb.size(); j++){

     totalRadTransProb = totalRadTransProb + transProb[j];

   }

   return totalRadTransProb;   
      
    }
 else {
   G4Exception( "shell not found" );
	
 }
  }
  else{
    G4Exception( "Z not found");
  } 
}

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(G4int z, G4int shellId)

{

  G4std::map<G4int,G4std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;
  
  pos = shellTable.find(z);
  
  if (pos!= shellTable.end()){
    
    G4std::vector<G4AtomicShell*> v = pos->second;
    
    G4int index = v.size();
    
    for (G4int i=0; i<v.size(); i++){
      
      G4int tmpId=(v[i])->ShellId();
      if (tmpId==shellId)
	{index = i;}
      else{}
    }  
    
    if (index< v.size()){
      
      G4AtomicShell* shell=v[index];
      
      G4DataVector transProb = shell->TransitionProbabilities();
      
      G4double totalRadTransProb = 0;
      
      for(G4int j = 0; j<transProb.size(); j++){

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

