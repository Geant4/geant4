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
// $Id: G4DopplerProfile.cc,v 1.1 2008-03-10 19:10:12 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

// The following deprecated header is included because <functional> seems not to be found on MGP's laptop
//#include "function.h"

// Constructor

G4DopplerProfile::G4DopplerProfile(G4int minZ, G4int maxZ)
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
G4DopplerProfile::~G4DopplerProfile()
{
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::iterator pos;
  for (pos = profileMap.begin(); pos != profileMap.end(); ++pos)
    {
      G4VEMDataSet* dataSet = (*pos).second;
      delete dataSet;
    }

  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::iterator pos2;
  for (pos2 = profilePdfMap.begin(); pos2 != profilePdfMap.end(); ++pos2)
    {
      G4VEMDataSet* dataSet = (*pos2).second;
      delete dataSet;
    }
}


size_t G4DopplerProfile::NumberOfProfiles(G4int Z) const
{
  G4int n = 0;

  if (Z>= zMin && Z <= zMax)
    {
      n = nShells[Z-1];
    }
  return n;
}


const G4VEMDataSet* G4DopplerProfile::Profiles(G4int Z)
{
  std::map<G4int,G4VEMDataSet*,std::less<G4int> >::const_iterator pos;
  if (Z < zMin || Z > zMax) G4Exception("G4DopplerProfile::Profiles - Z outside boundaries");
  pos = profileMap.find(Z);
  G4VEMDataSet* dataSet = (*pos).second;
  return dataSet;
}


const G4VEMDataSet* G4DopplerProfile::Profile(G4int Z, G4int shellIndex)
{
  const G4VEMDataSet* profis = Profiles(Z);
  const G4VEMDataSet* profi = profis->GetComponent(shellIndex);
  return profi;
}


void G4DopplerProfile::PrintData() const
{
  /*
  for (G4int Z = zMin; Z <= zMax; Z++)
    {
      G4cout << "---- Shell data for Z = "
	     << Z
	     << " ---- "
	     << G4endl;
      G4int nSh = nShells[Z-1];
      std::map<G4int,std::vector<G4double>*,std::less<G4int> >::const_iterator posId;
      posId = idMap.find(Z);
      std::vector<G4double>* ids = (*posId).second;
      std::map<G4int,G4DataVector*,std::less<G4int> >::const_iterator posE;
      posE = bindingMap.find(Z);
      G4DataVector* energies = (*posE).second;
      for (G4int i=0; i<nSh; i++)
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
		std::map<G4int,std::vector<G4double>*,std::less<G4int> >::const_iterator posOcc;
		posOcc = occupancyPdfMap.find(Z);
                std::vector<G4double> probs = *((*posOcc).second);
                G4double prob = probs[i];
		G4cout << "- Probability = " << prob;
	      }
	    G4cout << G4endl;
	}
      G4cout << "-------------------------------------------------" 
	     << G4endl;
    }
  */
}


void G4DopplerProfile::LoadData(const G4String& fileName)
{ 
  // Build the complete string identifying the file with the data set
  /*
  std::ostringstream ost;
  
  ost << fileName << ".dat";
  
  G4String name(ost.str());
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4EMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("G4DopplerProfile - data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception(excep);
    }

  G4double a = 0;
  G4int k = 1;
  G4int s = 0;
  
  G4int Z = 1;
  G4DataVector* energies = new G4DataVector;
  std::vector<G4double>* ids = new std::vector<G4double>;

  do {
    file >> a;
    G4int nColumns = 2;
    if (a == -1)
      {
	if (s == 0)
	  {
	    // End of a shell data set
	    idMap[Z] = ids;
            bindingMap[Z] = energies;
            G4int n = ids->size();
	    nShells.push_back(n);
	    // Start of new shell data set
	    ids = new std::vector<G4double>;
            energies = new G4DataVector;
            Z++;	    
	  }      
	s++;
	if (s == nColumns)
	{
	  s = 0;
	}
      }
    else if (a == -2)
      {
	// End of file; delete the empty vectors created when encountering the last -1 -1 row
	delete energies;
	delete ids;
	//nComponents = components.size();
      }
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

  // For Doppler broadening: the data set contains shell occupancy and binding energy for each shell
  // Build additional map with probability for each shell based on its occupancy

  if (occupancyData)
    {
      // Build cumulative from raw shell occupancy

      for (G4int Z=zMin; Z <= zMax; Z++)
	{
	  std::vector<G4double> occupancy = ShellIdVector(Z);

	  std::vector<G4double>* prob = new std::vector<G4double>;
	  G4double scale = 1. / G4double(Z);

	  prob->push_back(occupancy[0] * scale);
	  for (size_t i=1; i<occupancy.size(); i++)
	    {
	      prob->push_back(occupancy[i]*scale + (*prob)[i-1]);
	    }
	  occupancyPdfMap[Z] = prob;

	  
	}
    }
  */
}


G4double G4DopplerProfile::RandomSelectMomentum(G4int Z, G4int ShellIndex) const
{
  /*
  if (Z < zMin || Z > zMax) G4Exception("G4DopplerProfile::RandomSelect - Z outside boundaries");

  G4int shellIndex = 0;    
  std::vector<G4double> prob = ShellVector(Z);
  G4double random = G4UniformRand();

  // std::vector<G4double>::const_iterator pos;
  // pos = lower_bound(prob.begin(),prob.end(),random);

  // Binary search the shell with probability less or equal random

  G4int nShells = NumberOfShells(Z);
  G4int upperBound = nShells;

  while (shellIndex <= upperBound) 
    {      G4int midShell = (shellIndex + upperBound) / 2;
      if ( random < prob[midShell] ) 
	upperBound = midShell - 1;
      else 
	shellIndex = midShell + 1;
    }  
  if (shellIndex >= nShells) shellIndex = nShells - 1;
  */
  return 0.;
}


void G4DopplerProfile::LoadBiggsP(const G4String& fileName)
{
  std::ostringstream ost;
  
  ost << fileName << ".dat";
  
  G4String name(ost.str());
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4EMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("G4DopplerProfile::LoadBiggsP data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception(excep);
    }

  G4double p;
  while(!file.eof()) 
    {
      file >> p;
      biggsP.push_back(p);
      G4cout << "Biggs p = " << p << G4endl;
    }

  if (biggsP.size() != nBiggs)
    G4Exception("G4DopplerProfile::LoadBiggsP - Number of momenta read in is not 31");
}


void G4DopplerProfile::LoadProfile(const G4String& fileName,G4int Z)
{
  std::ostringstream ost;
  
  ost << fileName << "-" << Z << ".dat";
  
  G4String name(ost.str());
  
  char* path = getenv("G4LEDATA");
  if (!path)
    { 
      G4String excep("G4EMDataSet - G4LEDATA environment variable not set");
      G4Exception(excep);
    }
  
  G4String pathString(path);
  G4String dirFile = pathString + name;
  std::ifstream file(dirFile);
  std::filebuf* lsdp = file.rdbuf();

  if (! (lsdp->is_open()) )
    {
      G4String s1("G4DopplerProfile::LoadProfile data file: ");
      G4String s2(" not found");
      G4String excep = s1 + dirFile + s2;
      G4Exception(excep);
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
	  //	  if (i == 16) G4cout << "profile = " << p << G4endl;
	}
      // Create G4EMDataSet for the current shell
      G4VDataSetAlgorithm* algo = interpolation->Clone();
      G4VEMDataSet* dataSet = new G4EMDataSet(Z, biggs, profi, algo, 1., 1.);
      //     dataSet->PrintData();
      // Add current shell profile component to G4CompositeEMDataSet for the current Z
      dataSetForZ->AddComponent(dataSet);
    }

  G4cout << Z << ") has " << nShell << " shells" << G4endl;
  // Fill in number of shells for the current Z
  nShells.push_back(nShell);
 
  profileMap[Z] = dataSetForZ;
}



G4double G4DopplerProfile::IntegrateProfile(const std::vector<G4double>& profileVector)
{
  return 0.;
}
