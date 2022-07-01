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
//
// -------------------------------------------------------------------

#include "G4DopplerProfile.hh"
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"

#include <fstream>
#include <sstream>
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DopplerProfile::G4DopplerProfile(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ)
{  
  LoadBiggsP("/doppler/p-biggs");

  for (G4int Z=zMin; Z<zMax+1; Z++)
    {
      LoadProfile("/doppler/profile",Z);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Destructor
G4DopplerProfile::~G4DopplerProfile()
{
  for (auto& pos : profileMap)
    {
      G4VEMDataSet* dataSet = pos.second;
      delete dataSet;
      dataSet = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

size_t G4DopplerProfile::NumberOfProfiles(G4int Z) const
{
  G4int n = 0;
  if (Z>= zMin && Z <= zMax) n = nShells[Z-1];
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VEMDataSet* G4DopplerProfile::Profiles(G4int Z) const
{
  if (Z < zMin || Z > zMax) 
    G4Exception("G4DopplerProfile::Profiles",
		    "em1005",FatalException,"Z outside boundaries");
  auto pos = profileMap.find(Z);
  G4VEMDataSet* dataSet = (*pos).second;
  return dataSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4VEMDataSet* G4DopplerProfile::Profile(G4int Z, G4int shellIndex) const
{
  const G4VEMDataSet* profis = Profiles(Z);
  const G4VEMDataSet* profi = profis->GetComponent(shellIndex);
  return profi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DopplerProfile::PrintData() const
{
  for (G4int Z=zMin; Z<zMax; Z++)
    {
      const G4VEMDataSet* profis = Profiles(Z);
      profis->PrintData();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DopplerProfile::LoadBiggsP(const G4String& fileName)
{
  std::ostringstream ost;  
  ost << fileName << ".dat";
  G4String name(ost.str());
  
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    { 
      G4Exception("G4DopplerProfile::LoadBiggsP",
		    "em0006",FatalException,"G4LEDATA environment variable not set");
      return;
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception("G4DopplerProfile::LoadBiggsP",
		    "em0003",FatalException,excep);
    }

  G4double p;
  while(!file.eof()) 
    {
      file >> p;
      biggsP.push_back(p);
    }

  // Make sure that the number of data loaded corresponds to the number in Biggs' paper
  if (biggsP.size() != nBiggs)
    G4Exception("G4DopplerProfile::LoadBiggsP",
		    "em1006",FatalException,"Number of momenta read in is not 31");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DopplerProfile::LoadProfile(const G4String& fileName,G4int Z)
{
  std::ostringstream ost;
  ost << fileName << "-" << Z << ".dat";
  G4String name(ost.str());
  
  const char* path = G4FindDataDir("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4LEDATA environment variable not set");
      G4Exception("G4DopplerProfile::LoadProfile",
		    "em0006",FatalException,excep);
      return;
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception("G4DopplerProfile::LoadProfile",
		    "em0003",FatalException,excep);
    }

  G4double p;
  G4int nShell = 0;

  // Create CompositeDataSet for the current Z
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation;
  G4VEMDataSet* dataSetForZ = new G4CompositeEMDataSet(interpolation,1.,1.,1,1);

  while (!file.eof()) 
    {
      nShell++;
      G4DataVector* profi = new G4DataVector;
      G4DataVector* biggs = new G4DataVector;

      // Read in profile data for the current shell
      for (size_t i=0; i<nBiggs; i++)
	{ 
	  file >> p;
	  profi->push_back(p);
          biggs->push_back(biggsP[i]);
	}

      // Create G4EMDataSet for the current shell
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4EMDataSet(Z, biggs, profi, algo, 1., 1., true);
     
      // Add current shell profile component to G4CompositeEMDataSet for the current Z
      dataSetForZ->AddComponent(dataSet);
    }

  // Fill in number of shells for the current Z
  nShells.push_back(nShell);
 
  profileMap[Z] = dataSetForZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 G4double G4DopplerProfile::RandomSelectMomentum(G4int Z, G4int shellIndex) const
{
  G4double value = 0.;
  const G4VEMDataSet* profis = Profiles(Z);
  value = profis->RandomSelect(shellIndex);
  return value;
}
