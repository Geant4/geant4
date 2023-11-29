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
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// Based on G4FluoData by Elena Guardincerri
// 
// Modified: 
// 30.07.02 VI Add select active Z + clean up against pedantic compiler
//
// -------------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "G4AugerData.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4AugerData::G4AugerData()
{
  numberOfVacancies.resize(105, 0);
  BuildAugerTransitionTable(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4AugerData::NumberOfVacancies(G4int Z) const
{
  return numberOfVacancies[Z];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4AugerData::VacancyId(G4int Z, G4int vacancyIndex) const
{
  G4int n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {
      G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0;
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    n = (G4int) dataSet[vacancyIndex].FinalShellId();
  }  
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Attention: this method wants the vacancy index, not the Id
std::size_t G4AugerData::NumberOfTransitions(G4int Z, G4int vacancyIndex) const
{
  std::size_t n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {
      G4Exception("G4AugerData::VacancyId()","de0002", JustWarning, "Energy deposited locally");
      return 0;
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) 
      {
	G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
	return 0;
      }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    n = dataSet[vacancyIndex].TransitionOriginatingShellIds()->size();
  }
  return  n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4AugerData::NumberOfAuger(G4int Z, G4int initIndex, G4int vacancyId) const
{
  std::size_t n = 0;
  if (initIndex<0 || initIndex>=numberOfVacancies[Z])
    {
      G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
      return 0;
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0;
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    const std::vector<G4int>* temp =
      dataSet[initIndex].AugerOriginatingShellIds(vacancyId);
    n = temp->size();
  }
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4AugerData::AugerShellId(G4int Z, G4int vacancyIndex, G4int transId, G4int augerIndex) const
{
  std::size_t n = 0;  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {
      G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
      return 0;
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0;
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    n = dataSet[vacancyIndex].AugerOriginatingShellId(augerIndex,transId);
  }
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4AugerData::StartShellId(G4int Z, G4int vacancyIndex, G4int transitionShellIndex) const
{
  G4int n = 0; 

  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z]) {
    G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
    return 0;
  }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0; 
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    n = dataSet[vacancyIndex].TransitionOriginatingShellId(transitionShellIndex);
  } 
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4AugerData::StartShellEnergy(G4int Z, G4int vacancyIndex, G4int transitionId, G4int augerIndex) const
{
  G4double energy = 0;
  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {
      G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
      return 0;
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0;
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    energy = dataSet[vacancyIndex].AugerTransitionEnergy(augerIndex,transitionId);
      
  }
  return energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4AugerData::StartShellProb(G4int Z, G4int vacancyIndex,G4int transitionId,G4int augerIndex) const
{
  G4double prob = 0;
    
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z]) 
    {
      G4Exception("G4AugerData::VacancyId()","de0002", FatalErrorInArgument,"");
      return 0;
    }
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end()) {
      G4Exception("G4AugerData::VacancyId()","de0004", FatalErrorInArgument,  "Check element");
      return 0;
    }
    std::vector<G4AugerTransition> dataSet = (*element).second;
    prob = dataSet[vacancyIndex].AugerTransitionProbability(augerIndex, transitionId);
  }
  return prob;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AugerTransition> G4AugerData::LoadData(G4int Z)
{ 
  // Build the complete string identifying the file with the data set
  std::ostringstream ost;
  if(Z != 0){
    ost << "au-tr-pr-"<< Z << ".dat";
  }
  else{
    ost << "au-tr-pr-"<<".dat"; 
  }
  G4String name(ost.str());
  
  const char* path = G4FindDataDir("G4LEDATA");
  if (nullptr == path)
    { 
      G4String excep = "G4AugerData::LoadData";
      G4Exception(excep,"em0006", FatalException,"" );
      std::vector<G4AugerTransition> a;
      return a;
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + "/auger/" + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4String excep = "G4AugerData::LoadData";
      G4String msg = "Missing" + dirFile;
      G4Exception(excep,"em0003", FatalException, msg);
    }
 
  G4double a = 0;
  G4int k = 1;
  G4int sLocal = 0;
  
  G4int vacId = 0;
  std::vector<G4int>* initIds = new std::vector<G4int>;
  std::vector<G4int>* newIds = new std::vector<G4int>;
  G4DataVector* transEnergies = new G4DataVector;
  G4DataVector* transProbabilities = new G4DataVector;
  std::vector<G4AugerTransition> augerTransitionVector;
  std::map<G4int,std::vector<G4int>,std::less<G4int> >* newIdMap = 
    new std::map<G4int,std::vector<G4int>,std::less<G4int> >;
  std::map<G4int,G4DataVector,std::less<G4int> >* newEnergyMap =
    new std::map<G4int,G4DataVector,std::less<G4int> >;
  std::map<G4int,G4DataVector,std::less<G4int> >* newProbabilityMap = 
    new std::map<G4int,G4DataVector,std::less<G4int> >;
  
  do {
    file >> a;
    G4int nColumns = 4;
    if (a == -1)
      {
	if (sLocal == 0)
	  {
	    // End of a shell data set
	    std::vector<G4int>::iterator vectorIndex = initIds->begin();
	    
	    vacId = *vectorIndex;
	    std::vector<G4int> identifiers;
	    for (vectorIndex = initIds->begin()+1 ; vectorIndex != initIds->end(); ++vectorIndex){
		identifiers.push_back(*vectorIndex);
	    }
	    vectorIndex = (initIds->end())-1;
	    G4int augerShellId = *(vectorIndex);
	    	    
	    (*newIdMap)[augerShellId] = *newIds;
	    (*newEnergyMap)[augerShellId] = *transEnergies;
	    (*newProbabilityMap)[augerShellId] = *transProbabilities;
	    
	    augerTransitionVector.push_back(G4AugerTransition(vacId, identifiers, 
							      newIdMap, newEnergyMap, newProbabilityMap));
	    // Now deleting all the variables I used, and creating new ones for the next shell
	    delete newIdMap;
	    delete newEnergyMap;
	    delete newProbabilityMap;
	    
	    G4int n = (G4int)initIds->size();	    
	    nInitShells.push_back(n);
	    numberOfVacancies[Z]++;
	    delete initIds;
	    delete newIds;
	    delete transEnergies;	    
	    delete transProbabilities;
	    initIds = new std::vector<G4int>;
	    newIds = new std::vector<G4int>;
	    transEnergies = new G4DataVector;
	    transProbabilities = new G4DataVector;
	    newIdMap = 
	      new std::map<G4int,std::vector<G4int>,std::less<G4int> >;
	    newEnergyMap = new std::map<G4int,G4DataVector,std::less<G4int> >;
	    newProbabilityMap =
	      new std::map<G4int,G4DataVector,std::less<G4int> >;	
	  }      
	++sLocal;
	if (sLocal == nColumns)
	  {
	    sLocal = 0;
	  }
      }
    else
      {
	
	if (k%nColumns == 3){
	  // 3rd column is the transition probabilities
	  transProbabilities->push_back(a);	  
	  ++k;
	}
	else if(k%nColumns == 2){	 
	  // 2nd column is new auger vacancy	  
	  G4int l = (G4int)a;
	  newIds->push_back(l);
	  ++k;
	}
	else if (k%nColumns == 1)
	  {
	    // 1st column is shell id
	    if(initIds->size() == 0) {
	      // if this is the first data of the shell, all the colums are equal 
		// to the shell Id; so we skip the next colums ang go to the next row
	      initIds->push_back((G4int)a);
	      // first line of initIds is the original shell of the vacancy
	      file >> a;
	      file >> a;
	      file >> a;
	      k = k+3;
	    }
	    else {	      
	      if((G4int)a != initIds->back()){
		if((initIds->size()) == 1) { 
		  initIds->push_back((G4int)a);
		}  
		else {
		  
		  
		  G4int augerShellId = 0;
		  augerShellId = initIds->back();
		  
		  (*newIdMap)[augerShellId] = *newIds;
		  (*newEnergyMap)[augerShellId] = *transEnergies;
		  (*newProbabilityMap)[augerShellId] = *transProbabilities;
		  delete newIds;
		  delete transEnergies;
		  delete transProbabilities;
		  newIds = new std::vector<G4int>;
		  transEnergies = new G4DataVector;
		  transProbabilities = new G4DataVector;
		  initIds->push_back((G4int)a);
		}
	      }
	    }
	    k++;    
	  }
	else if (k%nColumns == 0)
	  {//fourth column is transition energies
	    G4double e = a * MeV; 
	    transEnergies->push_back(e);
	    k=1;
	  }
      }
  }
  while (a != -2); // end of file
  file.close();
  delete initIds;
  delete newIds;
  delete transEnergies;
  delete transProbabilities;
  delete newIdMap ;
  delete newEnergyMap;
  delete newProbabilityMap;
  return augerTransitionVector; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4AugerData::BuildAugerTransitionTable()
{
  for (G4int element = 6; element < 105; ++element) { 
    augerTransitionTable.insert(trans_Table::value_type(element,LoadData(element)));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4AugerData::PrintData(G4int Z) 
{
  for (G4int i = 0; i < numberOfVacancies[Z]; ++i)
    {
      G4cout << "---- TransitionData for the vacancy nb "
	     <<i
	     <<" of the atomic number elemnt " 
	     << Z
	     <<"----- "
	     <<G4endl;
      
      for (G4int k = 0; k<=(G4int)NumberOfTransitions(Z,i); ++k)
	{ 
	  G4int id = StartShellId(Z,i,k);
	  
	  for (G4int a = 0; a <= (G4int)NumberOfAuger(Z,i,id); ++a) {
	    G4double e = StartShellEnergy(Z,i,id,a)/MeV;
	    G4double p = StartShellProb(Z,i,id,a);
	    std::size_t augerId = AugerShellId(Z, i, id, a);
	    G4cout << k <<") Shell id: " << id <<G4endl;
	    G4cout << "    Auger Originatig Shell Id :"<< augerId <<G4endl;
	    G4cout << " - Transition energy = " << e << " MeV "<<G4endl;
	    G4cout << " - Transition probability = " << p <<G4endl;
	  }
	}
      G4cout << "-------------------------------------------------" 
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4AugerTransition* 
G4AugerData::GetAugerTransition(G4int Z, G4int vacancyShellIndex)
{
  std::vector<G4AugerTransition>* dataSet = &augerTransitionTable[Z];
  std::vector<G4AugerTransition>::iterator vectorIndex =
    dataSet->begin() + vacancyShellIndex;
  
  G4AugerTransition* augerTransition = &(*vectorIndex);
  return augerTransition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::vector<G4AugerTransition>* G4AugerData::GetAugerTransitions(G4int Z)
{
  std::vector<G4AugerTransition>* dataSet = &augerTransitionTable[Z];
  return dataSet;
}





















