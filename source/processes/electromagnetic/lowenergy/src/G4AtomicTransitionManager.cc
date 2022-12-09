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
#include "G4AutoLock.hh"
namespace { G4Mutex AtomicTransitionManagerMutex = G4MUTEX_INITIALIZER; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4AtomicTransitionManager* G4AtomicTransitionManager::instance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4AtomicTransitionManager* G4AtomicTransitionManager::Instance()
{
  if (instance == nullptr) {
    instance = new G4AtomicTransitionManager();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4AtomicTransitionManager::G4AtomicTransitionManager()
 : augerData(nullptr),
   verboseLevel(0),
   isInitialized(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4AtomicTransitionManager::~G4AtomicTransitionManager()
{ 
  delete augerData;
 
  for (auto& pos : shellTable){
    std::vector<G4AtomicShell*>vec = pos.second;
    std::size_t vecSize = vec.size();   
    for (std::size_t i=0; i< vecSize; ++i){
      G4AtomicShell* shell = vec[i];
      delete shell;     
    }
  }
 
  for (auto& ppos : transitionTable)
    {
      std::vector<G4FluoTransition*>vec = ppos.second;
      std::size_t vecSize=vec.size();
   
      for (std::size_t i=0; i< vecSize; ++i){
	G4FluoTransition* transition = vec[i];
	delete transition;      
      }
    }   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4AtomicShell* 
G4AtomicTransitionManager::Shell(G4int Z, size_t shellIndex) const
{   
  auto pos = shellTable.find(Z);
  
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// This function gives, upon Z and the Index of the initial shell where 
// the vacancy is, the radiative transition that can happen (originating 
// shell, energy, probability)
const G4FluoTransition* 
G4AtomicTransitionManager::ReachableShell(G4int Z,size_t shellIndex) const
{
  auto pos = transitionTable.find(Z);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
const G4AugerTransition* 
G4AtomicTransitionManager::ReachableAugerShell(G4int Z, 
					       G4int vacancyShellIndex) const
{
  return augerData->GetAugerTransition(Z,vacancyShellIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4int G4AtomicTransitionManager::NumberOfShells (G4int Z) const
{
  auto pos = shellTable.find(Z);

  std::size_t res = 0;
  if (pos != shellTable.cend()){

    res = ((*pos).second).size();

  } else {
    G4ExceptionDescription ed;
    ed << "No deexcitation for Z= " << Z;
    G4Exception("G4AtomicTransitionManager::NumberOfShells()","de0001",
		FatalException, ed, "");
  } 
  return (G4int)res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function returns the number of possible radiative transitions for 
// the atom with atomic number Z i.e. the number of shell in wich a vacancy 
// can be filled with a radiative transition 
G4int G4AtomicTransitionManager::NumberOfReachableShells(G4int Z) const
{
  auto pos = transitionTable.find(Z);
  std::size_t res = 0;
  if (pos!= transitionTable.cend())
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
  return (G4int)res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
// This function returns the number of possible NON-radiative transitions 
// for the atom with atomic number Z i.e. the number of shell in wich a 
// vacancy can be filled with a NON-radiative transition 
G4int G4AtomicTransitionManager::NumberOfReachableAugerShells(G4int Z)const 
{
  return (G4int)augerData->NumberOfVacancies(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4double G4AtomicTransitionManager::TotalRadiativeTransitionProbability(
         G4int Z, size_t shellIndex) const
{
  auto pos = transitionTable.find(Z);
  G4double totalRadTransProb = 0.0;

  if (pos!= transitionTable.end())
    {
      std::vector<G4FluoTransition*> v = (*pos).second;
      
      if (shellIndex < v.size())
	{
	  G4FluoTransition* transition = v[shellIndex];
	  G4DataVector transProb = transition->TransitionProbabilities();
	
	  for (size_t j=0; j<transProb.size(); ++j) 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4AtomicTransitionManager::Initialise()
{
  if(isInitialized) { return; }
  G4AutoLock l(&AtomicTransitionManagerMutex);

  if(isInitialized) { return; }
  isInitialized = true;

  // Selection of fluorescence files
  
  G4String defaultDirectory = "/fluor";
  G4String fluoDirectory = defaultDirectory;
  G4String bindingDirectory = defaultDirectory;
  G4EmFluoDirectory fdir = G4EmParameters::Instance()->FluoDirectory();
  G4int zLim = zMax + 1;
  if(fdir == fluoBearden) {
    zMax = 100;
    supTableLimit = 100;
    bindingDirectory = fluoDirectory = "/fluor_Bearden";
  } else if(fdir == fluoANSTO) {
    zLim = 93;
    fluoDirectory = "/fluor_ANSTO";
  } else if(fdir == fluoXDB_EADL) {
    zMax = 100;
    supTableLimit = 100;
    bindingDirectory = fluoDirectory = "/fluor_XDB_EADL";
  }
 
  // infTableLimit is initialized to 6 because EADL lacks data for Z<=5
  G4ShellData* shellManager = new G4ShellData(1, zMax, false);
  shellManager->LoadData(bindingDirectory + "/binding");

  // initialization of the data for auger effect
  augerData = new G4AugerData;
  
  // Fills shellTable with the data from EADL, identities and binding 
  // energies of shells
  for (G4int Z = zMin; Z <= zMax; ++Z)
    {    
      std::vector<G4AtomicShell*> vectorOfShells;  
      G4int shellIndex = 0; 

      G4int numberOfShells = (G4int)shellManager->NumberOfShells(Z);
      for (shellIndex = 0; shellIndex<numberOfShells; ++shellIndex) 
	{ 
	  G4int shellId = shellManager->ShellId(Z,shellIndex);
	  G4double bindingEnergy = shellManager->BindingEnergy(Z,shellIndex);
	  G4AtomicShell * shell = new G4AtomicShell(shellId,bindingEnergy);
	  vectorOfShells.push_back(shell);
	}
      shellTable[Z] = vectorOfShells;
    }
  
  // Fills transitionTable with the data on identities, transition 
  // energies and transition probabilities
  G4String dir = fluoDirectory;
  for (G4int Znum= infTableLimit; Znum<=supTableLimit; ++Znum)
    {  
      if (Znum == zLim) { dir = defaultDirectory; }
      G4FluoData* fluoManager = new G4FluoData(dir);
      std::vector<G4FluoTransition*> vectorOfTransitions;
      fluoManager->LoadData(Znum);
    
      G4int numberOfVacancies = (G4int)fluoManager->NumberOfVacancies();
      for(G4int vacancyIndex = 0; vacancyIndex<numberOfVacancies;  
	  ++vacancyIndex)
	{
	  std::vector<G4int> vectorOfIds;
	  G4DataVector vectorOfEnergies;
	  G4DataVector vectorOfProbabilities;
	
	  G4int finalShell = fluoManager->VacancyId(vacancyIndex);
	  G4int numberOfTransitions = (G4int)
                fluoManager->NumberOfTransitions(vacancyIndex);
	  for (G4int origShellIndex = 0; origShellIndex < numberOfTransitions;
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
  l.unlock();
}
