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
// $Id: G4FluoDataData.cc,v 1.2 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 16 Sept 2001  First committed to cvs
//
// -------------------------------------------------------------------

#include "G4FluoData.hh"
#include "G4DataVector.hh"
#include "g4std/fstream"
#include "g4std/strstream"

G4FluoData::G4FluoData()
{
  numberOfVacancies=0; 
}

G4FluoData::~G4FluoData()
{ 
 G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::iterator pos;

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

size_t G4FluoData::NumberOfVacancies() const
{
  return numberOfVacancies;
}

G4int G4FluoData::VacancyId(G4int vacancyIndex) const
{
  G4int n = -1;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4FluoData::vacancyIndex outside boundaries");}
  else
    {
      G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
      pos = idMap.find(vacancyIndex);
      if (pos!= idMap.end())
	{ G4DataVector dataSet = (*(*pos).second);
	n = (G4int) dataSet[0];
	
	}
    }
  return n;
}

size_t G4FluoData::NumberOfTransitions(G4int vacancyIndex) const
{
  G4int n = 0;
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4FluoData::vacancyIndex outside boundaries");}
  else
    {
      n = nInitShells[vacancyIndex]-1;
      //-1 is necessary because the elements of the vector nInitShells
      //include also the vacancy shell:
      // -1 subtracts this last one
  }
 return n;
}
G4int G4FluoData::StartShellId(G4int initIndex,G4int vacancyIndex) const
{
 G4int n = -1;

 if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4FluoData::vacancyIndex outside boundaries");}
 else
   {
     G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
    
     pos = idMap.find(vacancyIndex);
     
     G4DataVector dataSet = *((*pos).second);
   
     G4int nData = dataSet.size();
 if (initIndex >= 0 && initIndex < nData)
	    {
	      n =  (G4int) dataSet[initIndex];
	    
	    }
   }
 return n;
}
 
G4double G4FluoData::StartShellEnergy(G4int initIndex,G4int vacancyIndex) const
{
  G4double n = -1;
  
  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4FluoData::vacancyIndex outside boundaries");}
 else
   {
     G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
     
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

G4double G4FluoData::StartShellProb(G4int initIndex,G4int vacancyIndex) const
{
  G4double n = -1;

  if (vacancyIndex<0 || vacancyIndex>=numberOfVacancies)
    {G4Exception("G4FluoData::vacancyIndex outside boundaries");}
  else
    {
     G4std::map<G4int,G4DataVector*,G4std::less<G4int> >::const_iterator pos;
     
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

void G4FluoData::LoadData(G4int Z)
{ 
  // Build the complete string identifying the file with the data set
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  if(Z != 0){
    ost << "fl-tr-pr-"<< Z << ".dat";
  }
  else{
    ost << "fl-tr-pr-"<<".dat"; 
  }
  G4String name(nameChar);
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4EMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String fluor("/fluor/");
  G4String dirFile = pathString + fluor + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) )
    {
      G4String excep("G4FluoData - data file: " + dirFile + " not found");
      G4Exception(excep);
    }
  
  G4double a = 0;
  G4int k = 1;
  G4int s = 0;
  
  G4int vacId = 0;
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
	    idMap[vacId] = initIds;
            energyMap[vacId] = transEnergies;
	    probabilityMap[vacId] = transProbabilities;
	    //	    G4double size=transProbabilities->size();
            G4int n = initIds->size();
	    
	    nInitShells.push_back(n);
	    numberOfVacancies++;
	    // Start of new shell data set
	    initIds = new G4DataVector;
            transEnergies = new G4DataVector;
	    transProbabilities = new G4DataVector;
            vacId++;	
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
	   
	    transProbabilities->push_back(a);
	    
	    k++;
	  }
	else if (k%nColumns == 1)
	  {
	  // 1st column is shell id
	       
	    initIds->push_back(a);
	    k++;    
	   
	  }
	else if (k%nColumns == 0)

	  {//third column is transition energies
	   
	    G4double e = a * MeV;
	    transEnergies->push_back(e);
	    
	    k=1;
	  }
      }
  } 
  while (a != -2); // end of file
  file.close();    
}

void G4FluoData::PrintData() 
{
  
  for (G4int i = 0; i <numberOfVacancies; i++)
    {
      G4cout << "---- TransitionData for the vacancy nb "
	     <<i
	     <<" ----- "
	     <<G4endl;
      
      for (size_t k = 0; k<=NumberOfTransitions(i); k++)
	{ 
	  G4int id = StartShellId(k,i);
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










