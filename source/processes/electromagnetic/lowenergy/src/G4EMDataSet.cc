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
// $Id: G4EMDataSet.cc,v 1.5 2001-10-08 07:48:57 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
//
// -------------------------------------------------------------------

#include "G4EMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "g4std/fstream"
#include "g4std/strstream"

// Constructor

G4EMDataSet::G4EMDataSet(G4int Z,
			 G4DataVector* points, 
			 G4DataVector* values,
			 G4VDataSetAlgorithm* interpolation,
			 G4double unitE, G4double unitData)
  :z(Z), energies(points), data(values), algorithm(interpolation)
{
  numberOfBins = energies->size();
  unit1 = unitE;
  unit2 = unitData;
  if (interpolation == 0) 
    G4Exception("G4EMDataSet::G4EMDataSet - interpolation algorithm = 0");
}

G4EMDataSet:: G4EMDataSet(G4int Z, 
			  const G4String& dataFile,
			  G4VDataSetAlgorithm* interpolation,
			  G4double unitE, G4double unitData)
      :z(Z), algorithm(interpolation)
{
  energies = new G4DataVector;
  data = new G4DataVector;
  unit1 = unitE;
  unit2 = unitData;  
  LoadData(dataFile);
  numberOfBins = energies->size();
  if (interpolation == 0) 
    G4Exception("G4EMDataSet::G4EMDataSet - interpolation algorithm = 0");

}

// Destructor

G4EMDataSet::~G4EMDataSet()
{ 
  delete algorithm;
  delete energies;
  delete data;
}


G4double G4EMDataSet::FindValue(G4double e, G4int id) const
{
  G4double value;
  G4double e0 = (*energies)[0];
  // Protections
  size_t bin = FindBinLocation(e);
  if (bin == numberOfBins)
    {
      //      G4cout << "WARNING - G4EMDataSet::FindValue: energy outside upper boundary"
      //     << G4endl;
      value = (*data)[bin];
    }
  else if (e <= e0)
    {
      //     G4cout << "WARNING - G4EMDataSet::FindValue: energy outside lower boundary"
      //     << G4endl;
      value = (*data)[0];
    }
  else
    {
      if (algorithm == 0) 
	G4Exception("G4EMDataSet::FindValue - interpolation algorithm = 0");
      value = algorithm->Calculate(e,bin,*energies,*data);
    }
  
  return value;
}

G4int G4EMDataSet::FindBinLocation(G4double energy) const
{
  // Protection against call outside allowed range
  G4double e0 = (*energies)[0];
  if (energy < e0)
    {
      //  G4cout << z
      //     << " - WARNING - G4EMDataSet::FindBinLocation called with argument " 
      //     << energy 
      //    << " outside lower limit " 
      //   << e0
      //   << "; replaced with lower limit" 
      //  << G4endl;
      energy = e0;
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
  
  if (! (lsdp->is_open()) )
	{
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
