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
//
// Author: Alfonso Mmantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
// Based on G4RDFluoData by Elena Guardincerri
// 
// Modified: 30.07.02 VI Add select active Z + clean up against pedantic compiler
//
// -------------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "G4RDAugerData.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"

G4RDAugerData::G4RDAugerData()
{

  G4int n = 0;
  G4int pos = 0;

    for (pos = 0 ; pos < 100; pos++) 
      {
	numberOfVacancies.push_back(n);
      }

  BuildAugerTransitionTable(); 


}

G4RDAugerData::~G4RDAugerData()
{ 
  /*
  std::map<G4int,std::vector<G4RDAugerTransition>,std::less<G4int> >::iterator pos;

  for (pos = augerTransitionTable.begin(); pos != augerTransitionTable.end(); pos++)
    {
      std::vector<G4RDAugerTransition> dataSet = (*pos).second;
      delete dataSet;
    }
  for (pos = energyMap.begin(); pos != energyMap.end(); pos++)
    {
      std::map<G4Int,G4DataVector*,std::less<G4int>>* dataMap = (*pos).second;
      for (pos2 = newIdProbabilityMap.begin(); pos2 != idMap.end(); pos2++)
	{
	  G4DataVector* dataSet = (*pos2).second;
	  delete dataSet;
	}
    }
  for (pos = probabilityMap.begin(); pos != probabilityMap.end(); pos++)
    {
      std::map<G4Int,G4DataVector*,std::less<G4int>>* dataMap = (*pos).second;
      for (pos2 = newIdProbabilityMap.begin(); pos2 != idMap.end(); pos2++)
	{
	  G4DataVector* dataSet = (*pos2).second;
	  delete dataSet;
	}
    }
  for (pos2 = newIdMap.begin(); pos2 != idMap.end(); pos2++)
    {
      G4DataVector* dataSet = (*pos2).second;
      delete dataSet;
    }
  for (pos2 = newIdEnergyMap.begin(); pos2 != idMap.end(); pos2++)
    {
      G4DataVector* dataSet = (*pos2).second;
      delete dataSet;
    }
  for (pos2 = newIdProbabilityMap.begin(); pos2 != idMap.end(); pos2++)
    {
      G4DataVector* dataSet = (*pos2).second;
      delete dataSet;
    }
  */

}

size_t G4RDAugerData::NumberOfVacancies(G4int Z) const
{
  return numberOfVacancies[Z];
}

G4int G4RDAugerData::VacancyId(G4int Z, G4int vacancyIndex) const
{

  G4int n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {G4Exception("G4RDAugerData::VacancyId()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::VacancyId()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    n = (G4int) dataSet[vacancyIndex].FinalShellId();
  }
  
  return n;
}



// Attenzione: questa funzione vuole l'indice della vacanza, non l'Id


size_t G4RDAugerData::NumberOfTransitions(G4int Z, G4int vacancyIndex) const
{
  G4int n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {G4Exception("G4RDAugerData::NumberOfTransitions()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::NumberOfTransitions()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    n = (G4int)dataSet[vacancyIndex].TransitionOriginatingShellIds()->size();
  }
 return  n;
}



size_t G4RDAugerData::NumberOfAuger(G4int Z, G4int initIndex, G4int vacancyId) const
{
  size_t n = 0;
  if (initIndex<0 || initIndex>=numberOfVacancies[Z])
    {G4Exception("G4RDAugerData::NumberOfAuger()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::NumberOfAuger()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    const std::vector<G4int>* temp =  dataSet[initIndex].AugerOriginatingShellIds(vacancyId);
    n = temp->size();
  }
  return n;
}

size_t G4RDAugerData::AugerShellId(G4int Z, G4int vacancyIndex, G4int transId, G4int augerIndex) const
{
  size_t n = 0;  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {G4Exception("G4RDAugerData::AugerShellId()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::AugerShellId()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    n = dataSet[vacancyIndex].AugerOriginatingShellId(augerIndex,transId);
  }
  return n;
}

G4int G4RDAugerData::StartShellId(G4int Z, G4int vacancyIndex, G4int transitionShellIndex) const
{
  G4int n = 0; 

  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z]) 
    {G4Exception("G4RDAugerData::StartShellId()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::StartShellId()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
     n = dataSet[vacancyIndex].TransitionOriginatingShellId(transitionShellIndex);
  }
   
 
 return n;
}

G4double G4RDAugerData::StartShellEnergy(G4int Z, G4int vacancyIndex, G4int transitionId, G4int augerIndex) const
{
  G4double energy = 0;
  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z])
    {G4Exception("G4RDAugerData::StartShellEnergy()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::StartShellEnergy()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    energy = dataSet[vacancyIndex].AugerTransitionEnergy(augerIndex,transitionId);
      
  }
  return energy;
}


G4double G4RDAugerData::StartShellProb(G4int Z, G4int vacancyIndex,G4int transitionId,G4int augerIndex) const
{
  G4double prob = 0;
    
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies[Z]) 
    {G4Exception("G4RDAugerData::StartShellProb()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else {
    trans_Table::const_iterator element = augerTransitionTable.find(Z);
    if (element == augerTransitionTable.end())
      {G4Exception("G4RDAugerData::StartShellProb()", "NoDataFound",
                   FatalException, "Data not loaded!");}
    std::vector<G4RDAugerTransition> dataSet = (*element).second;
    prob = dataSet[vacancyIndex].AugerTransitionProbability(augerIndex, transitionId);



  }
     return prob;
}

std::vector<G4RDAugerTransition> G4RDAugerData::LoadData(G4int Z)
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
  
    char* path = std::getenv("G4LEDATA");
    if (!path)
      { 
	G4String excep = "G4LEDATA environment variable not set";
	G4Exception("G4RDAugerData::LoadData()", "InvalidSetup",
                    FatalException, excep);
      }
  
    G4String pathString(path);
    G4String dirFile = pathString + "/auger/" + name;
    std::ifstream file(dirFile);
    std::filebuf* lsdp = file.rdbuf();
  
    if (! (lsdp->is_open()) )
      {
	G4String excep = "Data file: " + dirFile + " not found!";
	G4Exception("G4RDAugerData::LoadData()", "DataNotFound",
                    FatalException, excep);
      }
 

    G4double a = 0;
    G4int k = 1;
    G4int s = 0;
  
    G4int vacId = 0;
    std::vector<G4int>* initIds = new std::vector<G4int>;
    std::vector<G4int>* newIds = new std::vector<G4int>;
    G4DataVector* transEnergies = new G4DataVector;
    G4DataVector* transProbabilities = new G4DataVector;
    std::vector<G4RDAugerTransition> augerTransitionVector;
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



	  if (s == 0)
	    {
	      // End of a shell data set
	    
	    
	    
	      std::vector<G4int>::iterator vectorIndex = initIds->begin();

	      vacId = *vectorIndex;
	    
	      //initIds->erase(vectorIndex);
	    


	      std::vector<G4int> identifiers;
	      for (vectorIndex = initIds->begin()+1 ; vectorIndex != initIds->end(); ++vectorIndex){

		identifiers.push_back(*vectorIndex);
	      }

	      vectorIndex = (initIds->end())-1;

	      G4int augerShellId = *(vectorIndex);
	    
	    
	      (*newIdMap)[augerShellId] = *newIds;
	      (*newEnergyMap)[augerShellId] = *transEnergies;
	      (*newProbabilityMap)[augerShellId] = *transProbabilities;

	      augerTransitionVector.push_back(G4RDAugerTransition(vacId, identifiers, newIdMap, newEnergyMap, newProbabilityMap));

	      // Now deleting all the variables I used, and creating new ones for the next shell

	      delete newIdMap;
	      delete newEnergyMap;
	      delete newProbabilityMap;
	    
	      G4int n = initIds->size();	    
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
	      newIdMap = new std::map<G4int,std::vector<G4int>,std::less<G4int> >;
	      newEnergyMap = new std::map<G4int,G4DataVector,std::less<G4int> >;
	      newProbabilityMap = new std::map<G4int,G4DataVector,std::less<G4int> >;	
	    


	    }      
	  s++;
	  if (s == nColumns)
	    {
	      s = 0;
	    }
	}
      else if (a == -2)
	{
	  // End of file; delete the empty vectors created 
	  //when encountering the last -1 -1 row
	  delete initIds;
	  delete newIds;
	  delete transEnergies;
	  delete transProbabilities;
	  delete newIdMap ;
	  delete newEnergyMap;
	  delete newProbabilityMap;	
	} 
      else
	{
	
	  if (k%nColumns == 3){
	    // 3rd column is the transition probabilities
	    transProbabilities->push_back(a);

	    k++;}
	  else if(k%nColumns == 2){	 
	    // 2nd column is new auger vacancy

	    // 2nd column is new auger vacancy

	    G4int l = (G4int)a;
	    newIds->push_back(l);


	    k++;
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

		//	      std::vector<G4int>::iterator vectorIndex = (initIds->end())-1;
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
    return augerTransitionVector;

}

void G4RDAugerData::BuildAugerTransitionTable()
{

  //  trans_Table::iterator pos = augerTransitionTable.begin();

  const G4MaterialTable* materialTable = G4Material::GetMaterialTable();

  G4int nMaterials = G4Material::GetNumberOfMaterials();

  G4DataVector activeZ;
  activeZ.clear();
  
  for (G4int m=0; m<nMaterials; m++) {

    const G4Material* material= (*materialTable)[m];        
    const G4ElementVector* elementVector = material->GetElementVector();
    const size_t nElements = material->GetNumberOfElements();
      
    for (size_t iEl=0; iEl<nElements; iEl++) {
      G4Element* element = (*elementVector)[iEl];
      G4double Z = element->GetZ();
      if (!(activeZ.contains(Z))) {
	activeZ.push_back(Z);
      }
    }
  }


  for (G4int element = 6; element < 101; element++)
    { 
      if(nMaterials == 0 || activeZ.contains(element)) {
        augerTransitionTable.insert(trans_Table::value_type(element,LoadData(element)));
      
        G4cout << "G4RDAugerData for Element no. " << element << " are loaded" << G4endl;
      //      PrintData(element);
      }    
    }
  
  G4cout << "AugerTransitionTable complete"<< G4endl;
}

void G4RDAugerData::PrintData(G4int Z) 
{
  
  for (G4int i = 0; i < numberOfVacancies[Z]; i++)
    {
      G4cout << "---- TransitionData for the vacancy nb "
	     <<i
	     <<" of the atomic number elemnt " 
	     << Z
	     <<"----- "
	     <<G4endl;
      
      for (size_t k = 0; k<=NumberOfTransitions(Z,i); k++)
	{ 
	  G4int id = StartShellId(Z,i,k);
	  
	  for (size_t a = 0; a <= NumberOfAuger(Z,i,id); a++) {
	    
	    G4double e = StartShellEnergy(Z,i,id,a) /MeV;
	    G4double p = StartShellProb(Z,i,id,a);
	    G4int augerId = AugerShellId(Z, i, id, a);
	    G4cout << k <<") Shell id: " << id <<G4endl;
	    G4cout << "    Auger Originatig Shell Id :"<< augerId <<G4endl;
	    G4cout << " - Transition energy = " << e << " MeV "<<G4endl;
	    G4cout   << " - Transition probability = " << p <<G4endl;
	  }
	}
      G4cout << "-------------------------------------------------" 
	     << G4endl;
    }
}
G4RDAugerTransition* G4RDAugerData::GetAugerTransition(G4int Z,G4int vacancyShellIndex)
    {
      std::vector<G4RDAugerTransition>* dataSet = &augerTransitionTable[Z];
      std::vector<G4RDAugerTransition>::iterator vectorIndex = dataSet->begin() + vacancyShellIndex;

      G4RDAugerTransition* augerTransition = &(*vectorIndex);
      return augerTransition;
    }
  
std::vector<G4RDAugerTransition>* G4RDAugerData::GetAugerTransitions(G4int Z)
  {
    std::vector<G4RDAugerTransition>* dataSet = &augerTransitionTable[Z];
    return dataSet;
  }
 




















