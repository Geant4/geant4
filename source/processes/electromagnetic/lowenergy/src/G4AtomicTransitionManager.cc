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
// 16 Sep 2001 E. Guardincerri  First Committed to cvs
//
// -------------------------------------------------------------------

#include "G4AtomicTransitionManager.hh"

G4AtomicTransitionManager::G4AtomicTransitionManager(G4int minZ, G4int maxZ, G4int limitInfTable,G4int limitSupTable)
  :zMin(minZ), zMax(maxZ),infTableLimit(limitInfTable),supTableLimit(limitSupTable)
{
  // infTableLimit is initialized to 6 because EADL lacks data for Z<=5
  G4ShellData* shellManager = new G4ShellData;
  
  shellManager->LoadData("/fluor/binding");
  
  // Fills shellTable with the data from EADL, identities and binding 
  // energies of shells
  for (G4int Z = zMin; Z<= zMax; Z++)
    {
      G4std::vector<G4AtomicShell*> vectorOfShells;  
    
      size_t numberOfShells=shellManager->NumberOfShells(Z);
      for (size_t shellIndex = 0; shellIndex<numberOfShells; shellIndex++)
	{ 
	  G4int shellId = shellManager->ShellId(Z,shellIndex);
	  G4double bindingEnergy = shellManager->BindingEnergy(Z,shellIndex);
	 
	  G4AtomicShell * shell = new G4AtomicShell(shellId,bindingEnergy);
	
	  vectorOfShells.push_back(shell);
	}
    
      //     shellTable.insert(G4std::make_pair(Z, vectorOfShells));
      shellTable[Z] = vectorOfShells;
    }
  
  // Fills transitionTable with the data from EADL, identities, transition 
  // energies and transition probabilities
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
    //      transitionTable.insert(G4std::make_pair(Znum, vectorOfTransitions));
    transitionTable[Znum] = vectorOfTransitions;
      
      delete fluoManager;
    }
  delete shellManager;
}

G4AtomicTransitionManager::~G4AtomicTransitionManager()
  
{ G4std::map<G4int,G4std::vector<G4AtomicShell*>,G4std::less<G4int> >::iterator pos;
 
 for (pos = shellTable.begin(); pos != shellTable.end(); pos++){
   
   G4std::vector< G4AtomicShell*>vec = (*pos).second;
   
   G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
     
     delete vec[i]; 
   }
   
 }
 
 G4std::map<G4int,G4std::vector<G4AtomicTransition*>,G4std::less<G4int> >::iterator ppos;
 
 for (ppos = transitionTable.begin(); ppos != transitionTable.end(); ppos++){
   
   G4std::vector< G4AtomicTransition*>vec = (*ppos).second;
   
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


const G4AtomicShell* G4AtomicTransitionManager::Shell(G4int Z, size_t shellIndex)
{ 
  G4std::map<G4int,G4std::vector<G4AtomicShell*>,G4std::less<G4int> >::iterator pos;
  
  pos = shellTable.find(Z);
  
  if (pos!= shellTable.end()){
    
    G4std::vector<G4AtomicShell*> v = (*pos).second;
    
    if (shellIndex<v.size()){
      
      return(v[shellIndex]);
      
    }
    else {
      G4Exception("G4AtomicTransitionManager:shell not found");
      return 0;
    }
  }
  else{
    G4Exception("G4AtomicTransitionManager:Z not found");
    return 0;
  } 
}

const G4AtomicTransition* G4AtomicTransitionManager:: ReachableShell(G4int Z,size_t shellIndex)
{
 G4std::map<G4int,G4std::vector<G4AtomicTransition*>,G4std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = (*pos).second;

    if (shellIndex<v.size()){
      
      return(v[shellIndex]);
      
    }
    else {
      G4Exception("G4AtomicTransitionManager:reachable shell not found");
      return 0;
    }
  }
  else{
    G4Exception("G4AtomicTransitionManager:Z not found");
    return 0;
  } 
}

G4int G4AtomicTransitionManager::NumberOfShells (G4int Z)
{

G4std::map<G4int,G4std::vector<G4AtomicShell*>,G4std::less<G4int> >::iterator pos;

  pos = shellTable.find(Z);

  if (pos!= shellTable.end()){

    G4std::vector<G4AtomicShell*> v = (*pos).second;

    return v.size();
  }

  else{
    G4Exception( "G4AtomicTransitionManager: Z not found" );
    return 0;
  } 
}

G4int G4AtomicTransitionManager::NumberOfReachableShells(G4int Z)
{
G4std::map<G4int,G4std::vector<G4AtomicTransition*>,G4std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = (*pos).second;

    return v.size();
  }

  else{
    G4Exception( "G4AtomicTransitionManager: Z not found" );
    return 0;
  } 
}
G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(G4int Z, size_t shellIndex)

{

G4std::map<G4int,G4std::vector<G4AtomicTransition*>,G4std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end()){

    G4std::vector<G4AtomicTransition*> v = (*pos).second;
  
    if (shellIndex<v.size()){

      G4AtomicTransition* transition=v[shellIndex];

   G4DataVector transProb = transition->TransitionProbabilities();

   G4double totalRadTransProb = 0;

   for(size_t j = 1; j<transProb.size(); j++){

     totalRadTransProb = totalRadTransProb + transProb[j];

   }

   return totalRadTransProb;   
      
    }
 else {
   G4Exception( "G4AtomicTransitionManager: shell not found" );
   return 0;
	
 }
  }
  else{
    G4Exception( "G4AtomicTransitionManager: Z not found");
    return 0;
  } 
}

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(G4int Z, size_t shellIndex)

{

  G4std::map<G4int,G4std::vector<G4AtomicTransition*>,G4std::less<G4int> >::iterator pos;
  
  pos = transitionTable.find(Z);
  
  if (pos!= transitionTable.end()){
    
    G4std::vector<G4AtomicTransition*> v = (*pos).second;
  
    
    if (shellIndex<v.size()){

      G4AtomicTransition* transition=v[shellIndex];
      
      G4DataVector transProb = transition->TransitionProbabilities();
      
      G4double totalRadTransProb = 0;
      
      for(size_t j = 1; j<transProb.size(); j++){

	totalRadTransProb = totalRadTransProb + transProb[j];
      }
      
      G4double totalNonRadTransProb= (1 - totalRadTransProb);
      
      return totalNonRadTransProb;    }
    
    else {
      G4Exception( "shell not found");
      return 0;
    }
  }
  else{
    G4Exception("Z not found");
    return 0;
  } 
}










