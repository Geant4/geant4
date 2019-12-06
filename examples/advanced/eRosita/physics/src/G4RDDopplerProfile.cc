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

#include "G4RDDopplerProfile.hh"
#include "G4DataVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"

#include <fstream>
#include <sstream>
#include "Randomize.hh"

// The following deprecated header is included because <functional> seems not to be found on MGP's laptop
//#include "function.h"

// Constructor

G4RDDopplerProfile::G4RDDopplerProfile(G4int minZ, G4int maxZ)
  : zMin(minZ), zMax(maxZ)
{  
  nBiggs = 31;

  LoadBiggsP("/doppler/p-biggs");

  for (G4int Z=zMin; Z<zMax+1; Z++)
    {
      LoadProfile("/doppler/profile",Z);
    }
}

// Destructor
G4RDDopplerProfile::~G4RDDopplerProfile()
{
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::iterator pos;
  for (pos = profileMap.begin(); pos != profileMap.end(); ++pos)
    {
      G4RDVEMDataSet* dataSet = (*pos).second;
      delete dataSet;
      dataSet = 0;
    }
}


size_t G4RDDopplerProfile::NumberOfProfiles(G4int Z) const
{
  G4int n = 0;
  if (Z>= zMin && Z <= zMax) n = nShells[Z-1];
  return n;
}


const G4RDVEMDataSet* G4RDDopplerProfile::Profiles(G4int Z) const
{
  std::map<G4int,G4RDVEMDataSet*,std::less<G4int> >::const_iterator pos;
  if (Z < zMin || Z > zMax)
    G4Exception("G4RDDopplerProfile::Profiles()", "OutOfRange",
                FatalException, "Z outside boundaries!");
  pos = profileMap.find(Z);
  G4RDVEMDataSet* dataSet = (*pos).second;
  return dataSet;
}


const G4RDVEMDataSet* G4RDDopplerProfile::Profile(G4int Z, G4int shellIndex) const
{
  const G4RDVEMDataSet* profis = Profiles(Z);
  const G4RDVEMDataSet* profi = profis->GetComponent(shellIndex);
  return profi;
}


void G4RDDopplerProfile::PrintData() const
{
  for (G4int Z=zMin; Z<zMax; Z++)
    {
      const G4RDVEMDataSet* profis = Profiles(Z);
      profis->PrintData();
    }
}


void G4RDDopplerProfile::LoadBiggsP(const G4String& fileName)
{
  std::ostringstream ost;  
  ost << fileName << ".dat";
  G4String name(ost.str());
  
  char* path = std::getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4LEDATA environment variable not set!");
      G4Exception("G4RDDopplerProfile::LoadBiggsP()",
                  "InvalidSetup", FatalException, excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("Data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception("G4RDDopplerProfile::LoadBiggsP()",
                  "DataNotFound", FatalException, excep);
    }

  G4double p;
  while(!file.eof()) 
    {
      file >> p;
      biggsP.push_back(p);
    }

  // Make sure that the number of data loaded corresponds to the number in Biggs' paper
  if (biggsP.size() != nBiggs)
    G4Exception("G4RDDopplerProfile::LoadBiggsP()", "InvalidCondition",
                FatalException, "Number of momenta read in is not 31!");
}


void G4RDDopplerProfile::LoadProfile(const G4String& fileName,G4int Z)
{
  std::ostringstream ost;
  ost << fileName << "-" << Z << ".dat";
  G4String name(ost.str());
  
  char* path = std::getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4LEDATA environment variable not set!");
      G4Exception("G4RDDopplerProfile::LoadProfile()",
                  "InvalidSetup", FatalException, excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("Data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception("G4RDDopplerProfile::LoadProfile()",
                  "DataNotFound", FatalException, excep);
    }

  G4double p;
  G4int nShell = 0;

  // Create CompositeDataSet for the current Z
  G4RDVDataSetAlgorithm* interpolation = new G4RDLogLogInterpolation;
  G4RDVEMDataSet* dataSetForZ = new G4RDCompositeEMDataSet(interpolation,1.,1.,1,1);

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
	  //	  if (i == 16) G4cout << "profile = " << p << G4endl;
	}

      // Create G4RDEMDataSet for the current shell
      G4RDVDataSetAlgorithm* algo = interpolation->Clone();
      G4RDVEMDataSet* dataSet = new G4RDEMDataSet(Z, biggs, profi, algo, 1., 1., true);
     
      // Add current shell profile component to G4RDCompositeEMDataSet for the current Z
      dataSetForZ->AddComponent(dataSet);
    }

  // Fill in number of shells for the current Z
  nShells.push_back(nShell);
 
  profileMap[Z] = dataSetForZ;
}


 G4double G4RDDopplerProfile::RandomSelectMomentum(G4int Z, G4int shellIndex) const
{
  G4double value = 0.;
  const G4RDVEMDataSet* profis = Profiles(Z);
  value = profis->RandomSelect(shellIndex);
  return value;
}
