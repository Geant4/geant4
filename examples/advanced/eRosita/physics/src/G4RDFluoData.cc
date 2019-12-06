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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 16 Sept 2001  First committed to cvs
//
// -------------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "G4RDFluoData.hh"
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"
#include "G4RDFluoTransition.hh"

G4RDFluoData::G4RDFluoData()
{
  numberOfVacancies=0; 
}

G4RDFluoData::~G4RDFluoData()
{ 
 std::map<G4int,G4DataVector*,std::less<G4int> >::iterator pos;

  for (pos = idMap.begin(); pos != idMap.end(); ++pos)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
    }
  for (pos = energyMap.begin(); pos != energyMap.end(); ++pos)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
    }
 for (pos = probabilityMap.begin(); pos != probabilityMap.end(); ++pos)
    {
      G4DataVector* dataSet = (*pos).second;
      delete dataSet;
    }
}

size_t G4RDFluoData::NumberOfVacancies() const
{
  return numberOfVacancies;
}

G4int G4RDFluoData::VacancyId(G4int vacancyIndex) const
{
  G4int n = -1;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4RDFluoData::VacancyId()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else
    {
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
      pos = idMap.find(vacancyIndex);
      if (pos!= idMap.end())
	{ G4DataVector dataSet = (*(*pos).second);
	n = (G4int) dataSet[0];
	
	}
    }
  return n;
}

size_t G4RDFluoData::NumberOfTransitions(G4int vacancyIndex) const
{
  G4int n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4RDFluoData::NumberOfTransitions()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else
    {
      n = nInitShells[vacancyIndex]-1;
      //-1 is necessary because the elements of the vector nInitShells
      //include also the vacancy shell:
      // -1 subtracts this last one
  }
 return n;
}
G4int G4RDFluoData::StartShellId(G4int initIndex, G4int vacancyIndex) const
{
 G4int n = -1;

 if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4RDFluoData::StartShellId()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
 else
   {
     std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
    
     pos = idMap.find(vacancyIndex);
     
     G4DataVector dataSet = *((*pos).second);
   
     G4int nData = dataSet.size();
     //The first Element of idMap's dataSets is the original shell of the vacancy, 
     //so we must start from the first element of dataSet
 if (initIndex >= 0 && initIndex < nData)
	    {
	      n =  (G4int) dataSet[initIndex+1];
	    
	    }
   }
 return n;
}
 
G4double G4RDFluoData::StartShellEnergy(G4int initIndex, G4int vacancyIndex) const
{
  G4double n = -1;
  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4RDFluoData::StartShellEnergy()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
 else
   {
     std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
     
     pos = energyMap.find(vacancyIndex);
     
     G4DataVector dataSet = *((*pos).second);
     
     G4int nData = dataSet.size();
     if (initIndex >= 0 && initIndex < nData)
       {
	 n =  dataSet[initIndex];
	 
       }
   }
  return n;
}

G4double G4RDFluoData::StartShellProb(G4int initIndex, G4int vacancyIndex) const
{
  G4double n = -1;

  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4RDFluoData::StartShellProb()", "OutOfRange",
                 FatalException, "VacancyIndex outside boundaries!");}
  else
    {
     std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator pos;
     
     pos = probabilityMap.find(vacancyIndex);
     
     G4DataVector dataSet = *((*pos).second);
     
     G4int nData = dataSet.size();
     if (initIndex >= 0 && initIndex < nData)
       {
	 n =  dataSet[initIndex];
	 
       }
    }
  return n;
}

void G4RDFluoData::LoadData(G4int Z)
{ 
  // Build the complete string identifying the file with the data set
  
  std::ostringstream ost;
  if(Z != 0){
    ost << "fl-tr-pr-"<< Z << ".dat";
  }
  else{
    ost << "fl-tr-pr-"<<".dat"; 
  }
  G4String name(ost.str());
  
  char* path = std::getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4LEDATA environment variable not set!");
      G4Exception("G4RDEMDataSet::LoadData()", "InvalidSetup",
                  FatalException, excep);
    }
  
  G4String pathString(path);
  G4String fluor("/fluor/");
  G4String dirFile = pathString + fluor + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4String excep = "Data file: " + dirFile + " not found";
      G4Exception("G4RDEMDataSet::LoadData()", "DataNotFound",
                  FatalException, excep);
    }
  
  G4double a = 0;
  G4int k = 1;
  G4int s = 0;
  
  G4int vacIndex = 0;
  G4DataVector* initIds = new G4DataVector;
  G4DataVector* transEnergies = new G4DataVector;
  G4DataVector* transProbabilities = new G4DataVector;
  
  do {
    file >> a;
    G4int nColumns = 3;
    if (a == -1)
      {
	if (s == 0)
	  {
	    // End of a shell data set
	    idMap[vacIndex] = initIds;
            energyMap[vacIndex] = transEnergies;
	    probabilityMap[vacIndex] = transProbabilities;
	    //	    G4double size=transProbabilities->size();
            G4int n = initIds->size();
	    
	    nInitShells.push_back(n);
	    numberOfVacancies++;
	    // Start of new shell data set
	    initIds = new G4DataVector;
            transEnergies = new G4DataVector;
	    transProbabilities = new G4DataVector;
            vacIndex++;	
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
	delete transEnergies;
	delete transProbabilities;
      } 
    else
      {
	
	if(k%nColumns == 2)
	  {	 
	    // 2nd column is transition  probabilities

	   if (a != -1) transProbabilities->push_back(a);
	    
	    k++;
	  }
	else if (k%nColumns == 1)
	  {
	    // 1st column is shell id
	    // if this is the first data of the shell, all the colums are equal 
	    // to the shell Id; so we skip the next colums ang go to the next row
	    if(initIds->size() == 0) {
	      if (a != -1) initIds->push_back((G4int)a);
	      file >> a;
	      file >> a;
	      k=k+2;
	    } 
	    else{ 
	      if (a != -1) initIds->push_back(a);
	    }
	    k++;    
	  }
	else if (k%nColumns == 0)

	  {//third column is transition energies

	    if (a != -1) 
	      {G4double e = a * MeV;
	      transEnergies->push_back(e);}
	   
	    k=1;
	  }
      }
  } 
  while (a != -2); // end of file
  file.close();    
}

void G4RDFluoData::PrintData() 
{
  
  for (G4int i = 0; i <numberOfVacancies; i++)
    {
      G4cout << "---- TransitionData for the vacancy nb "
	     <<i
	     <<" ----- "
	     <<G4endl;
      
      for (size_t k = 0; k<NumberOfTransitions(i); k++)
	{ 
	  G4int id = StartShellId(k,i);
	// let's start from 1 because the first (index = 0) element of the vector
	// is the id of the intial vacancy
	  G4double e = StartShellEnergy(k,i) /MeV;
	  G4double p = StartShellProb(k,i); 
	  G4cout << k <<") Shell id: " << id <<G4endl;
	  G4cout << " - Transition energy = " << e << " MeV "<<G4endl;
	  G4cout   << " - Transition probability = " << p <<G4endl;
	  
	}
      G4cout << "-------------------------------------------------" 
	     << G4endl;
    }
}
