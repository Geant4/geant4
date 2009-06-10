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

G4AtomicTransitionManager::G4AtomicTransitionManager(G4int minZ, G4int maxZ, 
  G4int limitInfTable,G4int limitSupTable)
  :zMin(minZ), 
  zMax(maxZ),
  infTableLimit(limitInfTable),
  supTableLimit(limitSupTable)
{
  // infTableLimit is initialized to 6 because EADL lacks data for Z<=5
  G4ShellData* shellManager = new G4ShellData;

  // initialization of the data for auger effect
  
  augerData = new G4AugerData;

  shellManager->LoadData("/fluor/binding");
  
  // Fills shellTable with the data from EADL, identities and binding 
  // energies of shells
  for (G4int Z = zMin; Z<= zMax; Z++)
    {
      std::vector<G4AtomicShell*> vectorOfShells;  
      size_t shellIndex = 0; 

      size_t numberOfShells=shellManager->NumberOfShells(Z);
      for (shellIndex = 0; shellIndex<numberOfShells; shellIndex++) 
	{ 
	  G4int shellId = shellManager->ShellId(Z,shellIndex);
	  G4double bindingEnergy = shellManager->BindingEnergy(Z,shellIndex);
	 
	  G4AtomicShell * shell = new G4AtomicShell(shellId,bindingEnergy);
	
	  vectorOfShells.push_back(shell);
	}
    
      //     shellTable.insert(std::make_pair(Z, vectorOfShells));
      shellTable[Z] = vectorOfShells;
    }
  
  // Fills transitionTable with the data from EADL, identities, transition 
  // energies and transition probabilities
  for (G4int Znum= infTableLimit; Znum<=supTableLimit; Znum++)
    {  G4FluoData* fluoManager = new G4FluoData;
    std::vector<G4FluoTransition*> vectorOfTransitions;
    fluoManager->LoadData(Znum);
    
    size_t numberOfVacancies = fluoManager-> NumberOfVacancies();
    
    for (size_t vacancyIndex = 0; vacancyIndex<numberOfVacancies;  vacancyIndex++)
      
      {
	std::vector<G4int>  vectorOfIds;
	G4DataVector vectorOfEnergies;
	G4DataVector vectorOfProbabilities;
	
	G4int finalShell = fluoManager->VacancyId(vacancyIndex);
	size_t numberOfTransitions = fluoManager->NumberOfTransitions(vacancyIndex);
	for (size_t origShellIndex = 0; origShellIndex < numberOfTransitions;
	     origShellIndex++)
	    
	  {
	    
	    G4int originatingShellId = fluoManager->StartShellId(origShellIndex,vacancyIndex);
	    
	    vectorOfIds.push_back(originatingShellId);
	    
	    G4double transitionEnergy = fluoManager->StartShellEnergy(origShellIndex,vacancyIndex);
	    vectorOfEnergies.push_back(transitionEnergy);
	    G4double transitionProbability = fluoManager->StartShellProb(origShellIndex,vacancyIndex);
	    vectorOfProbabilities.push_back(transitionProbability);
	  }
	  G4FluoTransition * transition = new G4FluoTransition (finalShell,vectorOfIds,
								    vectorOfEnergies,vectorOfProbabilities);
	  vectorOfTransitions.push_back(transition); 
      }
    //      transitionTable.insert(std::make_pair(Znum, vectorOfTransitions));
    transitionTable[Znum] = vectorOfTransitions;
      
      delete fluoManager;
    }
  delete shellManager;
}

G4AtomicTransitionManager::~G4AtomicTransitionManager()
  
{ 

  delete augerData;

std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;
 
 for (pos = shellTable.begin(); pos != shellTable.end(); pos++){
   
   std::vector< G4AtomicShell*>vec = (*pos).second;
   
   G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
     G4AtomicShell* shell = vec[i];
     delete shell;     
   }
   
 }
 
 std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::iterator ppos;
 
 for (ppos = transitionTable.begin(); ppos != transitionTable.end(); ppos++){
   
   std::vector<G4FluoTransition*>vec = (*ppos).second;
   
   G4int vecSize=vec.size();
   
   for (G4int i=0; i< vecSize; i++){
	 G4FluoTransition* transition = vec[i];
	 delete transition;      
   }
   
 }   
 
}

G4AtomicTransitionManager* G4AtomicTransitionManager::instance = 0;

G4AtomicTransitionManager* G4AtomicTransitionManager::Instance()
{
  if (instance == 0)
    {
      instance = new G4AtomicTransitionManager;
     
    }
  return instance;
}


G4AtomicShell* G4AtomicTransitionManager::Shell(G4int Z, size_t shellIndex) const
{ 
  std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::const_iterator pos;
  
  pos = shellTable.find(Z);
  
  if (pos!= shellTable.end())
    {
      std::vector<G4AtomicShell*> v = (*pos).second;
      if (shellIndex<v.size())
	{
	  return(v[shellIndex]);
	}
      else 
	{
	  size_t lastShell = v.size();
	  G4cout << "G4AtomicTransitionManager::Shell - Z = " 
		 << Z << ", shellIndex = " << shellIndex 
		 << " not found; number of shells = " << lastShell << G4endl;
	  //  G4Exception("G4AtomicTransitionManager:shell not found");
	  if (lastShell > 0)
	    {
	      return v[lastShell - 1];
	    }
	  else
	    {	    
	      return 0;
	    }
	}
    }
  else
    {
      G4Exception("G4AtomicTransitionManager:Z not found");
      return 0;
    } 
}

// This function gives, upon Z and the Index of the initial shell where te vacancy is, 
// the radiative transition that can happen (originating shell, energy, probability)

const G4FluoTransition* G4AtomicTransitionManager::ReachableShell(G4int Z,size_t shellIndex) const
{
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::const_iterator pos;
  pos = transitionTable.find(Z);
  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;      
      if (shellIndex < v.size()) return(v[shellIndex]);
      else {
	G4Exception("G4AtomicTransitionManager:reachable shell not found");
	return 0;
      }
  }
  else{
    G4cout << "G4AtomicTransitionMagare warning: No fluorescence or Auger for Z=" << Z << G4endl;
    G4cout << "Absorbed enrgy deposited locally" << G4endl;

    //    G4Exception("G4AtomicTransitionManager:Z not found");
    return 0;
  } 
}

const G4AugerTransition* G4AtomicTransitionManager::ReachableAugerShell(G4int Z, G4int vacancyShellIndex) const
{
  
  G4AugerTransition* augerTransition = augerData->GetAugerTransition(Z,vacancyShellIndex);
  return augerTransition;
}



G4int G4AtomicTransitionManager::NumberOfShells (G4int Z) const
{

std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::const_iterator pos;

  pos = shellTable.find(Z);

  if (pos!= shellTable.end()){

    std::vector<G4AtomicShell*> v = (*pos).second;

    return v.size();
  }

  else{
    G4cout << "G4AtomicTransitionMagare warning: No fluorescence or Auger for Z=" << Z << G4endl;
    G4cout << "Absorbed enrgy deposited locally" << G4endl;

    //    G4Exception("G4AtomicTransitionManager:Z not found");
    return 0;
  } 
}

// This function returns the number of possible radiative transitions for the atom with atomic number Z
// i.e. the number of shell in wich a vacancy can be filled with a radiative transition 

G4int G4AtomicTransitionManager::NumberOfReachableShells(G4int Z) const
{
std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::const_iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;
      return v.size();
    }
  else
    {
      G4cout << "G4AtomicTransitionMagare warning: No fluorescence or Auger for Z=" << Z << G4endl;
      G4cout << "Absorbed enrgy deposited locally" << G4endl;

      //    G4Exception("G4AtomicTransitionManager:Z not found");
      return 0;
    } 
}

// This function returns the number of possible NON-radiative transitions for the atom with atomic number Z
// i.e. the number of shell in wich a vacancy can be filled with a NON-radiative transition 

G4int G4AtomicTransitionManager::NumberOfReachableAugerShells(G4int Z)const 
{
  G4int n = augerData->NumberOfVacancies(Z);
  return n;
}



G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(G4int Z, 
									size_t shellIndex)

{
std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::iterator pos;

  pos = transitionTable.find(Z);

  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;
      
    if (shellIndex < v.size())
      {
	G4FluoTransition* transition = v[shellIndex];
	G4DataVector transProb = transition->TransitionProbabilities();
	G4double totalRadTransProb = 0;
	
	for (size_t j = 0; j<transProb.size(); j++) // AM -- corrected, it was 1
	{
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
    G4cout << "G4AtomicTransitionMagare warning: No fluorescence or Auger for Z=" << Z << G4endl;
    G4cout << "Absorbed enrgy deposited locally" << G4endl;

    //    G4Exception("G4AtomicTransitionManager:Z not found");

    return 0;
  } 
}

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(G4int Z, size_t shellIndex)

{

  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::iterator pos;
  
  pos = transitionTable.find(Z);
  
  if (pos!= transitionTable.end()){
    
    std::vector<G4FluoTransition*> v = (*pos).second;
  
    
    if (shellIndex<v.size()){

      G4FluoTransition* transition=v[shellIndex];
      G4DataVector transProb = transition->TransitionProbabilities();
      G4double totalRadTransProb = 0;
      
      for(size_t j = 0; j<transProb.size(); j++) // AM -- Corrected, was 1
	{
	  totalRadTransProb = totalRadTransProb + transProb[j];
	}
      
      if (totalRadTransProb > 1) {
      G4Exception( "Wrong Total Probability");
      return 0;
}
      G4double totalNonRadTransProb= (1 - totalRadTransProb);
      
      return totalNonRadTransProb;    }
    
    else {
      G4Exception( "shell not found");
      return 0;
    }
  }
  else{
    G4cout << "G4AtomicTransitionMagare warning: No fluorescence or Auger for Z=" << Z << G4endl;
    G4cout << "Absorbed enrgy deposited locally" << G4endl;

    //    G4Exception("G4AtomicTransitionManager:Z not found");
    return 0;
  } 
}










