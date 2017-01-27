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
#include "G4EmParameters.hh"
#include "G4FluoData.hh"
#include "G4AugerData.hh"

G4AtomicTransitionManager* G4AtomicTransitionManager::instance = 0;

G4AtomicTransitionManager* G4AtomicTransitionManager::Instance()
{
  if (instance == 0) {
    instance = new G4AtomicTransitionManager();
  }
  return instance;
}

G4AtomicTransitionManager::G4AtomicTransitionManager()
 : augerData(0),
   zMin(1), 
   zMax(100),
   infTableLimit(6),
   supTableLimit(100),
   isInitialized(false),
   verboseLevel(0)
{}

G4AtomicTransitionManager::~G4AtomicTransitionManager()
{ 
  delete augerData;

  std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::iterator pos;
 
  for(pos = shellTable.begin(); pos != shellTable.end(); ++pos){
   
    std::vector<G4AtomicShell*>vec = (*pos).second;
    G4int vecSize = vec.size();
   
    for (G4int i=0; i< vecSize; ++i){
      G4AtomicShell* shell = vec[i];
      delete shell;     
    }
  }
 
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::iterator ppos;
  for (ppos = transitionTable.begin(); ppos != transitionTable.end(); ++ppos){
   
    std::vector<G4FluoTransition*>vec = (*ppos).second;
   
    G4int vecSize=vec.size();
   
    for (G4int i=0; i< vecSize; i++){
      G4FluoTransition* transition = vec[i];
      delete transition;      
    }
  }   
}

G4AtomicShell* 
G4AtomicTransitionManager::Shell(G4int Z, size_t shellIndex) const
{ 
  std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::const_iterator pos;
  
  pos = shellTable.find(Z);
  
  if (pos!= shellTable.end())
    {
      std::vector<G4AtomicShell*> v = (*pos).second;
      if (shellIndex < v.size()) { return v[shellIndex]; }

      else 
	{
	  size_t lastShell = v.size();
	  G4ExceptionDescription ed;
	  ed << "No de-excitation for Z= " << Z 
	     << "  shellIndex= " << shellIndex
	     << ">=  numberOfShells= " << lastShell;
	  if (verboseLevel > 0) 
              G4Exception("G4AtomicTransitionManager::Shell()","de0001",
		JustWarning,ed," AtomicShell not found");
	  if (lastShell > 0) { return v[lastShell - 1]; }
	}
    }
  else
    {
      G4ExceptionDescription ed;
      ed << "No de-excitation for Z= " << Z 
	 << "  shellIndex= " << shellIndex
	 << ". AtomicShell not found - check if data are uploaded";
      G4Exception("G4AtomicTransitionManager::Shell()","de0001",
		  FatalException,ed,"");
    } 
  return 0;
}

// This function gives, upon Z and the Index of the initial shell where 
// the vacancy is, the radiative transition that can happen (originating 
// shell, energy, probability)

const G4FluoTransition* 
G4AtomicTransitionManager::ReachableShell(G4int Z,size_t shellIndex) const
{
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::const_iterator pos;
  pos = transitionTable.find(Z);
  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;      
      if (shellIndex < v.size()) { return(v[shellIndex]); }

      else {
	G4ExceptionDescription ed;
	ed << "No fluo transition for Z= " << Z 
	   << "  shellIndex= " << shellIndex;
	G4Exception("G4AtomicTransitionManager::ReachebleShell()","de0002", 
		    FatalException,ed,"");
      }
    } 
  else 
    {
      G4ExceptionDescription ed;
      ed << "No transition table for Z= " << Z 
	 << "  shellIndex= " << shellIndex;
      G4Exception("G4AtomicTransitionManager::ReachableShell()","de0001",
		  FatalException,ed,"");
    } 
  return 0;
}

const G4AugerTransition* 
G4AtomicTransitionManager::ReachableAugerShell(G4int Z, 
					       G4int vacancyShellIndex) const
{
  return augerData->GetAugerTransition(Z,vacancyShellIndex);
}

G4int G4AtomicTransitionManager::NumberOfShells (G4int Z) const
{
  std::map<G4int,std::vector<G4AtomicShell*>,std::less<G4int> >::const_iterator pos;
  pos = shellTable.find(Z);

  G4int res = 0;
  if (pos != shellTable.end()){

    res = ((*pos).second).size();

  } else {
    G4ExceptionDescription ed;
    ed << "No deexcitation for Z= " << Z;
    G4Exception("G4AtomicTransitionManager::NumberOfShells()","de0001",
		FatalException, ed, "");
  } 
  return res;
}

// This function returns the number of possible radiative transitions for 
// the atom with atomic number Z i.e. the number of shell in wich a vacancy 
// can be filled with a radiative transition 
G4int G4AtomicTransitionManager::NumberOfReachableShells(G4int Z) const
{
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::const_iterator pos;
  pos = transitionTable.find(Z);
  G4int res = 0;
  if (pos!= transitionTable.end())
    {
      res = ((*pos).second).size();
    }
  else
    {
      G4ExceptionDescription ed;
      ed << "No deexcitation for Z= " << Z
	 << ", so energy deposited locally";
      G4Exception("G4AtomicTransitionManager::NumberOfReachebleShells()",
		  "de0001",FatalException,ed,"");
    } 
  return res;
}

// This function returns the number of possible NON-radiative transitions 
// for the atom with atomic number Z i.e. the number of shell in wich a 
// vacancy can be filled with a NON-radiative transition 

G4int G4AtomicTransitionManager::NumberOfReachableAugerShells(G4int Z)const 
{
  return augerData->NumberOfVacancies(Z);
}

G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(
         G4int Z, size_t shellIndex) const
{
  std::map<G4int,std::vector<G4FluoTransition*>,std::less<G4int> >::const_iterator pos;

  pos = transitionTable.find(Z);
  G4double totalRadTransProb = 0.0;

  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;
      
      if (shellIndex < v.size())
	{
	  G4FluoTransition* transition = v[shellIndex];
	  G4DataVector transProb = transition->TransitionProbabilities();
	
	  for (size_t j=0; j<transProb.size(); ++j) // AM -- corrected, it was 1
	    {
	      totalRadTransProb += transProb[j];
	    }
	}
      else 
	{
	  G4ExceptionDescription ed;
	  ed << "Zero transition probability for Z=" << Z 
	     << "  shellIndex= " << shellIndex;
	  G4Exception(
	  "G4AtomicTransitionManager::TotalRadiativeTransitionProbability()",
	  "de0002",FatalException,"Incorrect de-excitation");
	}
    }
  else
    {
      G4ExceptionDescription ed;
      ed << "No deexcitation for Z=" << Z 
	 << "  shellIndex= " << shellIndex;
      G4Exception(
      "G4AtomicTransitionManager::TotalRadiativeTransitionProbability()",
      "de0001",FatalException,ed,"Cannot compute transition probability");
    } 
  return totalRadTransProb;
}

G4double G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability(
                 G4int Z, size_t shellIndex) const
{
  G4double prob = 1.0 - TotalRadiativeTransitionProbability(Z, shellIndex);
  if(prob > 1.0 || prob < 0.0) {
    G4ExceptionDescription ed;
    ed << "Total probability mismatch Z= " << Z 
       << "  shellIndex= " << shellIndex
       << "  prob= " << prob;
    G4Exception( 
    "G4AtomicTransitionManager::TotalNonRadiativeTransitionProbability()",
    "de0003",FatalException,ed,"Cannot compute non-radiative probability");
    return 0.0;
  }
  return prob;    
}

#include "G4AutoLock.hh"
namespace { G4Mutex AtomicTransitionManagerMutex = G4MUTEX_INITIALIZER; }

void G4AtomicTransitionManager::Initialise()
{
  G4AutoLock l(&AtomicTransitionManagerMutex);

  //G4cout << "!!! G4AtomicTransitionManager::Initialise " << isInitialized 
  // << G4endl;
  if(isInitialized) { return; }
  isInitialized = true;

  // Selection of fluorescence files
  G4String fluoDirectory = (G4EmParameters::Instance()->BeardenFluoDir()?
			    "/fluor_Bearden":"/fluor");

  // infTableLimit is initialized to 6 because EADL lacks data for Z<=5
  G4ShellData* shellManager = new G4ShellData;
  shellManager->LoadData(fluoDirectory+"/binding");

  // initialization of the data for auger effect
  augerData = new G4AugerData;
  
  // Fills shellTable with the data from EADL, identities and binding 
  // energies of shells
  for (G4int Z = zMin; Z<= zMax; ++Z) 
    {
      std::vector<G4AtomicShell*> vectorOfShells;  
      size_t shellIndex = 0; 

      size_t numberOfShells = shellManager->NumberOfShells(Z);
      // G4cout << "For Z= " << Z << "  " << numberOfShells << " shells" << G4endl;
      for (shellIndex = 0; shellIndex<numberOfShells; ++shellIndex) 
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
  for (G4int Znum= infTableLimit; Znum<=supTableLimit; ++Znum)
    {  
      G4FluoData* fluoManager = new G4FluoData(fluoDirectory);
      std::vector<G4FluoTransition*> vectorOfTransitions;
      fluoManager->LoadData(Znum);
    
      size_t numberOfVacancies = fluoManager-> NumberOfVacancies();
    
      for(size_t vacancyIndex = 0; vacancyIndex<numberOfVacancies;  
	  ++vacancyIndex)
	{
	  std::vector<G4int>  vectorOfIds;
	  G4DataVector vectorOfEnergies;
	  G4DataVector vectorOfProbabilities;
	
	  G4int finalShell = fluoManager->VacancyId(vacancyIndex);
	  size_t numberOfTransitions = 
	    fluoManager->NumberOfTransitions(vacancyIndex);
	  for (size_t origShellIndex = 0; origShellIndex < numberOfTransitions;
	       ++origShellIndex)
	    {
	      G4int originatingShellId = 
		fluoManager->StartShellId(origShellIndex,vacancyIndex);
	      vectorOfIds.push_back(originatingShellId);
	    
	      G4double transitionEnergy = 
		fluoManager->StartShellEnergy(origShellIndex,vacancyIndex);
	      vectorOfEnergies.push_back(transitionEnergy);
	      G4double transitionProbability = 
		fluoManager->StartShellProb(origShellIndex,vacancyIndex);
	      vectorOfProbabilities.push_back(transitionProbability);
	    }
	  G4FluoTransition* transition = 
	    new G4FluoTransition (finalShell,vectorOfIds,
				  vectorOfEnergies,vectorOfProbabilities);
	  vectorOfTransitions.push_back(transition); 
	}
      transitionTable[Znum] = vectorOfTransitions;
      delete fluoManager;
    }
  delete shellManager;
}





