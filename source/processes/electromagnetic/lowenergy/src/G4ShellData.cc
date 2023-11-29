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
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created
// 26 Dec 2010 V.Ivanchenko Fixed Coverity warnings   
//
// -------------------------------------------------------------------

#include "G4ShellData.hh"
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <valarray>
#include <functional>
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ShellData::G4ShellData(G4int minZ, G4int maxZ, G4bool isOccupancy)
  : zMin(minZ), zMax(maxZ), occupancyData(isOccupancy)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ShellData::~G4ShellData()
{
  for (auto& pos : idMap)
    {
      std::vector<G4double>* dataSet = pos.second;
      delete dataSet;
    }
  for (auto& pos2 : bindingMap) 
    {
      G4DataVector* dataSet = pos2.second;
      delete dataSet;
    }

  if (occupancyData)
    {
      for (auto& pos3 : occupancyPdfMap)
	{
	  std::vector<G4double>* dataSet = pos3.second;
	  delete dataSet;
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

std::size_t G4ShellData::NumberOfShells(G4int Z) const
{
  G4int z = Z - 1;
  G4int n = 0;

  if (Z>= zMin && Z <= zMax)
    {
      n = nShells[z];
    }
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4double>& G4ShellData::ShellIdVector(G4int Z) const
{
  if (Z < zMin || Z > zMax) {    
    G4Exception("G4ShellData::ShellIdVector","de0001",FatalErrorInArgument, "Z outside boundaries");
  }  
  auto pos = idMap.find(Z);
  std::vector<G4double>* dataSet = (*pos).second;
  return *dataSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const std::vector<G4double>& G4ShellData::ShellVector(G4int Z) const
{
  if (Z < zMin || Z > zMax) 
    G4Exception("G4ShellData::ShellVector()","de0001",JustWarning,"Z outside boundaries");
  auto pos = occupancyPdfMap.find(Z);
  std::vector<G4double>* dataSet = (*pos).second;
  return *dataSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4ShellData::ShellId(G4int Z, G4int shellIndex) const
{
  G4int n = -1;

  if (Z >= zMin && Z <= zMax)
    {
      auto pos = idMap.find(Z);
      if (pos!= idMap.end())
	{
	  std::vector<G4double> dataSet = *((*pos).second);
	  G4int nData = (G4int)dataSet.size();
	  if (shellIndex >= 0 && shellIndex < nData)
	    {
	      n = (G4int) dataSet[shellIndex];
	    }
	}
    }
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ShellData::ShellOccupancyProbability(G4int Z, G4int shellIndex) const
{
  G4double prob = -1.;
  if (Z >= zMin && Z <= zMax)
    {
      auto pos = idMap.find(Z);
      if (pos!= idMap.end())
	{
	  std::vector<G4double> dataSet = *((*pos).second);
	  G4int nData = (G4int)dataSet.size();
	  if (shellIndex >= 0 && shellIndex < nData)
	    {
	      prob = dataSet[shellIndex];
	    }
	}
    }
  return prob;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ShellData::BindingEnergy(G4int Z, G4int shellIndex)  const
{
  G4double value = 0.;

  if (Z >= zMin && Z <= zMax)
    {
      auto pos = bindingMap.find(Z);
      if (pos!= bindingMap.end())
	{
	  G4DataVector dataSet = *((*pos).second);
	  G4int nData = (G4int)dataSet.size();
	  if (shellIndex >= 0 && shellIndex < nData)
	    {
	      value = dataSet[shellIndex];
	    }
	}
    }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ShellData::PrintData() const
{
  for (G4int Z = zMin; Z <= zMax; Z++)
    {
      G4cout << "---- Shell data for Z = "
	     << Z
	     << " ---- "
	     << G4endl;
      G4int nSh = nShells[Z-1];
      auto posId = idMap.find(Z);
      std::vector<G4double>* ids = (*posId).second;
      auto posE = bindingMap.find(Z);
      G4DataVector* energies = (*posE).second;
      for (G4int i=0; i<nSh; ++i)
	{
	  G4int id = (G4int) (*ids)[i];
	  G4double e = (*energies)[i] / keV;
	  G4cout << i << ") ";

	  if (occupancyData) 
	    {
	      G4cout << " Occupancy: ";
	    }
	  else 
	    {
	      G4cout << " Shell id: ";
	    }
	  G4cout << id << " - Binding energy = "
		 << e << " keV ";
	    if (occupancyData)
	      {
		auto posOcc = occupancyPdfMap.find(Z);
                std::vector<G4double> probs = *((*posOcc).second);
                G4double prob = probs[i];
		G4cout << "- Probability = " << prob;
	      }
	    G4cout << G4endl;
	}
      G4cout << "-------------------------------------------------" 
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ShellData::LoadData(const G4String& fileName)
{ 
  // Build the complete string identifying the file with the data set
  std::ostringstream ost;
  ost << fileName << ".dat";
  G4String name(ost.str());
  
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4ShellData::LoadData()");
      G4Exception(excep,"em0006",FatalException,"Please set G4LEDATA");
      return;
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {

      G4String excep = "G4ShellData::LoadData()";
      G4String msg = "data file: " + dirFile + " not found";
      G4Exception(excep, "em0003",FatalException, msg );
      return;
    }

  G4double a = 0;
  G4int k = 1;
  G4int sLocal = 0;
  
  G4int Z = 1;
  G4DataVector* energies = new G4DataVector;
  std::vector<G4double>* ids = new std::vector<G4double>;

  do {
    file >> a;
    G4int nColumns = 2;
    if (a == -1)
      {
	if (sLocal == 0)
	  {
	    // End of a shell data set
	    idMap[Z] = ids;
            bindingMap[Z] = energies;
            G4int n = (G4int)ids->size();
	    nShells.push_back(n);
	    // Start of new shell data set
	    
	    ids = new std::vector<G4double>;
            energies = new G4DataVector;
            Z++;	    
	  }      
	sLocal++;
	if (sLocal == nColumns)
	{
	  sLocal = 0;
	}
      }

    // moved out of the do-while since might go to a leak. 
    //    else if (a == -2)
    //      {
    // End of file; delete the empty vectors created when encountering the last -1 -1 row
    //	delete energies;
    //	delete ids;
    //nComponents = components.size();
    //      }
    else
      {
	// 1st column is shell id
	if(k%nColumns != 0)
	  {	    
	    ids->push_back(a);
	    k++;
	  }
	else if (k%nColumns == 0)
	  {
	    // 2nd column is binding energy
	    G4double e = a * MeV;
	    energies->push_back(e);
	    k = 1;
	  }
      }
  } while (a != -2); // end of file
  file.close();    
  delete energies;
  delete ids;

  // For Doppler broadening: the data set contains shell occupancy and binding energy for each shell
  // Build additional map with probability for each shell based on its occupancy

  if (occupancyData)
    {
      // Build cumulative from raw shell occupancy
      for (G4int ZLocal=zMin; ZLocal <= zMax; ++ZLocal)
	{
	  std::vector<G4double> occupancy = ShellIdVector(ZLocal);

	  std::vector<G4double>* prob = new std::vector<G4double>;
	  G4double scale = 1. / G4double(ZLocal);

	  prob->push_back(occupancy[0] * scale);
	  for (std::size_t i=1; i<occupancy.size(); ++i)
	    {
	      prob->push_back(occupancy[i]*scale + (*prob)[i-1]);
	    }
	  occupancyPdfMap[ZLocal] = prob;
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4ShellData::SelectRandomShell(G4int Z) const
{
  if (Z < zMin || Z > zMax) 
    G4Exception("G4ShellData::SelectrandomShell","de0001",FatalErrorInArgument, "Z outside boundaries");

  G4int shellIndex = 0;    
  std::vector<G4double> prob = ShellVector(Z);
  G4double random = G4UniformRand();

  // Binary search the shell with probability less or equal random
  G4int nShellsLocal = (G4int)NumberOfShells(Z);
  G4int upperBound = nShellsLocal;

  while (shellIndex <= upperBound) 
    {
      G4int midShell = (shellIndex + upperBound) / 2;
      if ( random < prob[midShell] ) 
	upperBound = midShell - 1;
      else 
	shellIndex = midShell + 1;
    }  
  if (shellIndex >= nShellsLocal) shellIndex = nShellsLocal - 1;

  return shellIndex;
}
