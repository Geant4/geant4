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
// $Id: G4EMDataSet.cc,v 1.2 2001-09-05 12:29:51 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31.07.01   MGP          Created
// 05.09.01 V.Ivanchenko Minor fix in FindValue
//
// -------------------------------------------------------------------

#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "g4std/fstream"
#include "g4std/strstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet::G4EMDataSet(G4int Z,
			 G4DataVector* points, 
			 G4DataVector* values,
			 const G4VDataSetAlgorithm* interpolation,
			 G4double unitE, G4double unitData)
  :z(Z), energies(points), data(values), algorithm(interpolation)
{
  numberOfBins = energies->size();
  unit1 = unitE;
  unit2 = unitData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet:: G4EMDataSet(G4int Z, 
			  const G4String& dataFile,
			  const G4VDataSetAlgorithm* interpolation,
			  G4double unitE, G4double unitData)
  :z(Z), algorithm(interpolation)
{
  energies = new G4DataVector;
  data = new G4DataVector;
  unit1 = unitE;
  unit2 = unitData;  
  LoadData(dataFile);
  numberOfBins = energies->size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EMDataSet::~G4EMDataSet()
{ 
  delete energies;
  delete data;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EMDataSet::FindValue(G4double e, G4int id) const
{
  G4double value;

  // Protections
  size_t bin = FindBinLocation(e);
  if (bin == numberOfBins)
    {
      value = (*data)[bin];
    }
  else if (bin == 0)
    {
      value = (*data)[bin];
    }
  else
    {
      value = algorithm->Calculate(e,bin,*energies,*data);
    }
  
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4EMDataSet::FindBinLocation(G4double energy) const
{
  // Protection against call outside allowed range
  if (energy < (*energies)[0])
    {
      G4cout << z
	     << " - WARNING - G4EMDataSet::FindBinLocation(e) called with e= "
             << energy 
	     << " outside lower limit " << (*energies)[0]
	     << " -> is replaced with lower limit" 
	     << G4endl;
      energy = (*energies)[0];
    }

  size_t lowerBound = 0;
  size_t upperBound = numberOfBins - 1;
  
  // Binary search
  while (lowerBound <= upperBound) 
    {
      size_t midBin = (lowerBound + upperBound)/2;
      if ( energy < (*energies)[midBin] ) upperBound = midBin-1;
      else lowerBound = midBin+1;
  }
  
  return upperBound;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::LoadData(const G4String& fileName)
{
  // Build the complete string identifying the file with the data set
  
  char nameChar[100] = {""};
  G4std::ostrstream ost(nameChar, 100, G4std::ios::out);
  
  ost << fileName << z << ".dat";
  
  G4String name(nameChar);
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep = "G4EMDataSet - G4LEDATA environment variable not set";
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + "/" + name;
  G4std::ifstream file(dirFile);
  G4std::filebuf* lsdp = file.rdbuf();
  
  if (! (lsdp->is_open()) ) {
      G4String excep = "G4EMDataSet - data file: " + dirFile + " not found";
      G4Exception(excep);
    }
  G4double a = 0;
  G4int k = 1;

  do
    {
      file >> a;
      G4int nColumns = 2;
      // The file is organized into two columns:
      // 1st column is the energy
      // 2nd column is the corresponding value
      // The file terminates with the pattern: -1   -1
      //                                       -2   -2
      if (a == -1 || a == -2)
	{
	}
      else
	{
	  if (k%nColumns != 0)
	    {	
	      G4double e = a * unit1;
	      energies->push_back(e);
	      k++;
	    }
	  else if (k%nColumns == 0)
	    {
	      G4double value = a * unit2;
	      data->push_back(value);
	      k = 1;
	    }
	}
      
    } while (a != -2); // end of file
  
  file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4EMDataSet::PrintData() const
{
  size_t size = numberOfBins;
  for (size_t i=0; i<size; i++)
    {
      G4double e = (*energies)[i]  / unit1;
      G4double sigma = (*data)[i] / unit2 ;
      G4cout << "Point: "
	     << e
	     << " - Data value : "
	     << sigma 
	     << G4endl; 
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
