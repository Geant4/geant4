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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4NuclearLevelManager
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 24 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//      
// -------------------------------------------------------------------

#include "G4NuclearLevelManager.hh"

#include "globals.hh"
#include "G4NuclearLevel.hh"
#include "G4ios.hh"
#include <stdlib.h>
#include "g4std/fstream"
#include "g4std/strstream"
#include "g4std/algorithm"

G4NuclearLevelManager::G4NuclearLevelManager():
  _nucleusA(0), _nucleusZ(0), _levels(0), _levelEnergy(0), _gammaEnergy(0), _probability(0)
{ }

G4NuclearLevelManager::G4NuclearLevelManager(G4int Z, G4int A): _nucleusA(A), _nucleusZ(Z)
{ 


  if (A <= 0 || Z <= 0 || Z > A )
    G4Exception("==== G4NuclearLevelManager ==== (Z,A) <0, or Z>A");

  _levels = 0;

  MakeLevels();
}

G4NuclearLevelManager::~G4NuclearLevelManager()
{ 
  if ( _levels ) {
    if (_levels->size()>0) 
    {
      G4std::vector<G4NuclearLevel*>::iterator pos;
      for(pos=_levels->begin(); pos!=_levels->end(); pos++)
        if (*pos) delete *pos;
      _levels->clear();
    }
    delete _levels;
   _levels = 0;
  }
}

void G4NuclearLevelManager::SetNucleus(G4int Z, G4int A)
{
  if (_nucleusZ != Z || _nucleusA != A)
    {
      _nucleusA = A;
      _nucleusZ = Z;
      MakeLevels();
    }

}

G4bool G4NuclearLevelManager::IsValid(G4int Z, G4int A) const
{
  G4bool valid = true;

  if (A < 0 || Z < 0 || A < Z) valid = false;

  char* env = getenv("G4LEVELGAMMADATA");
  if (env == 0)
    G4Exception("G4NuclearLevelManager - Please set the G4LEVELGAMMADATA environment variable");
  G4String dirName(env);
  char name[100] = {""};
  G4std::ostrstream ost(name, 100, G4std::ios::out);
  ost << dirName << "/" << "z" << Z << ".a" << A;
  G4String file(name); 

#ifdef G4USE_STD_NAMESPACE  
  G4std::ifstream inFile(file, G4std::ios::in);
#else
  ifstream inFile(file, ios::in|ios::nocreate);
#endif
  if (! inFile) valid = false;  

  return valid;
}


G4int G4NuclearLevelManager::NumberOfLevels() const
{
  G4int n = 0;
  if (_levels != 0) n = _levels->size();
  return n;
}


const G4PtrLevelVector* G4NuclearLevelManager::GetLevels() const
{
  return _levels;
}


const G4NuclearLevel* G4NuclearLevelManager::NearestLevel(G4double energy, G4double eDiffMax) const
{
  G4int iNear = -1;
  
  G4double diff = 9999. * GeV;
  if (_levels != 0)
    {
      unsigned int i = 0;
      for (i=0; i<_levels->size(); i++)
	{
	  G4double e = (*_levels)[i]->Energy();
	  G4double eDiff = abs(e - energy);
	  if (eDiff < diff && eDiff <= eDiffMax)
	    { 
	      diff = eDiff; 
	      iNear = i;
	    }
	}
    }
  if (_levels != 0 && iNear >= 0 && iNear < static_cast<G4int>(_levels->size()) )
    { return (*_levels)[iNear]; }
  else
    { return 0; }
}


G4double G4NuclearLevelManager::MinLevelEnergy() const
{
  G4double eMin = 9999.*GeV;
  if (_levels != 0)
    {
      if (_levels->size() > 0) eMin = _levels->front()->Energy(); 
    }
  return eMin;
}


G4double G4NuclearLevelManager::MaxLevelEnergy() const
{
  G4double eMax = 0.;
  if (_levels != 0)
    {
      if (_levels->size() > 0) eMax = _levels->back()->Energy(); 
    }
  return eMax;
}


const G4NuclearLevel* G4NuclearLevelManager::HighestLevel() const
{
  if (_levels!= 0 && _levels->size() > 0) return _levels->front(); 
  else return 0; 
}


const G4NuclearLevel* G4NuclearLevelManager::LowestLevel() const
{
  if (_levels != 0 && _levels->size() > 0) return _levels->back();
  else return 0;
}


G4bool G4NuclearLevelManager::Read(G4std::ifstream& dataFile)
{
  const G4double minProbability = 0.001;

  G4bool result = true;

  if (dataFile >> _levelEnergy)
    {
      dataFile >> _gammaEnergy >> _probability >> _polarity >> _halfLife
	       >> _angularMomentum;
      _levelEnergy *= keV;
      _gammaEnergy *= keV;
      _halfLife *= second;

      // The following adjustment is needed to take care of anomalies in 
      // data files, where some transitions show up with relative probability
      // zero
      if (_probability < minProbability) _probability = minProbability;

      // G4cout << "Read " << _levelEnergy << " " << _gammaEnergy << " " << _probability << G4endl;
    }
  else
    {
      result = false;
    }

  return result;
}


void G4NuclearLevelManager::MakeLevels()
{
  char* env = getenv("G4LEVELGAMMADATA");
  if (env == 0)
    G4Exception("G4NuclearLevelManager: please set the G4LEVELGAMMADATA environment variable");
  G4String dirName(env);
  char name[100] = {""};
  G4std::ostrstream ost(name, 100, G4std::ios::out);
  ost << dirName << "/" << "z" << _nucleusZ << ".a" << _nucleusA;
  G4String file(name); 
#ifdef G4USE_STD_NAMESPACE  
  G4std::ifstream inFile(file, G4std::ios::in);
#else
  ifstream inFile(file, ios::in|ios::nocreate);
#endif
  if (! inFile) 
    {
      //      G4cout << " G4NuclearLevelManager: (" << _nucleusZ << "," << _nucleusA 
      //  	     << ") does not have LevelsAndGammas file" << G4endl;
      return;
    }

  if (_levels != 0)
    {
      if (_levels->size()>0) 
      {
        G4std::vector<G4NuclearLevel*>::iterator pos;
        for(pos=_levels->begin(); pos!=_levels->end(); pos++)
          if (*pos) delete *pos;
        _levels->clear();
      }
      delete _levels;
    }

  _levels = new G4PtrLevelVector;

  G4DataVector eLevel;
  G4DataVector eGamma;
  G4DataVector wGamma;
  G4DataVector pGamma; // polarity
  G4DataVector hLevel; // half life
  G4DataVector aLevel; // angular momentum

  while (Read(inFile))
    {
      eLevel.push_back(_levelEnergy);
      eGamma.push_back(_gammaEnergy);
      wGamma.push_back(_probability);
      pGamma.push_back(_polarity);
      hLevel.push_back(_halfLife);
      aLevel.push_back(_angularMomentum);
    }

  // ---- MGP ---- Don't forget to close the file 
  inFile.close();

  G4int nData = eLevel.size();

  //  G4cout << " ==== MakeLevels ===== " << nData << " data read " << G4endl;

  G4double thisLevelEnergy = eLevel[0];
  G4double thisLevelHalfLife = 0.;
  G4double thisLevelAngMom = 0.;
  G4DataVector thisLevelEnergies;
  G4DataVector thisLevelWeights;
  G4DataVector thisLevelPolarities;

  G4double e = -1.;
  G4int i;
  for (i=0; i<nData; i++)
    {
      e = eLevel[i];
      if (e != thisLevelEnergy)
      {
	//	  G4cout << "Making a new level... " << e << " " 
	//		 << thisLevelEnergies.entries() << " " 
	//		 << thisLevelWeights.entries() << G4endl;
	
	G4NuclearLevel* newLevel = new G4NuclearLevel(thisLevelEnergy,
						      thisLevelHalfLife,
						      thisLevelAngMom,
						      thisLevelEnergies,
						      thisLevelWeights,
						      thisLevelPolarities);
	  _levels->push_back(newLevel);
	  // Reset data vectors
	  thisLevelEnergies.clear();
	  thisLevelWeights.clear();
	  thisLevelPolarities.clear();
	  thisLevelEnergy = e;
      }
      // Append current data
      thisLevelEnergies.push_back(eGamma[i]);
      thisLevelWeights.push_back(wGamma[i]);
      thisLevelPolarities.push_back(pGamma[i]);
      thisLevelHalfLife = hLevel[i];
      thisLevelAngMom = aLevel[i];
    }
  // Make last level
  if (e > 0.)
    {
      G4NuclearLevel* newLevel = new G4NuclearLevel(e,thisLevelHalfLife,
						    thisLevelAngMom,
						    thisLevelEnergies,
						    thisLevelWeights,
						    thisLevelPolarities);
      _levels->push_back(newLevel);
    }
#ifdef G4USE_OSPACE
  G4std::sort(_levels->begin(), _levels->end());
#else
  G4std::stable_sort(_levels->begin(), _levels->end());
#endif
  return;
}


void G4NuclearLevelManager::PrintAll()
{
  G4int nLevels = 0;
  if (_levels != 0) nLevels = _levels->size();

  G4cout << " ==== G4NuclearLevelManager ==== (" << _nucleusZ << ", " << _nucleusA
	 << ") has " << nLevels << " levels" << G4endl
	 << "Highest level is at energy " << MaxLevelEnergy() << " MeV "
	 << G4endl << "Lowest level is at energy " << MinLevelEnergy()
	 << " MeV " << G4endl;

  G4int i = 0;
  for (i=0; i<nLevels; i++)
    { (*_levels)[i]->PrintAll(); }
}


G4NuclearLevelManager::G4NuclearLevelManager(const G4NuclearLevelManager &right)
{
  _levelEnergy = right._levelEnergy;
  _gammaEnergy = right._gammaEnergy;
  _probability = right._probability;
  _polarity = right._polarity;
  _halfLife = right._halfLife;
  _angularMomentum = right._angularMomentum;
  _nucleusA = right._nucleusA;
  _nucleusZ = right._nucleusZ;
  if (right._levels != 0)   
    {
      _levels = new G4PtrLevelVector;
      G4int n = right._levels->size();
      G4int i;
      for (i=0; i<n; i++)
	{
	  _levels->push_back(new G4NuclearLevel(*(right._levels->operator[](i))));
	}
#ifdef G4USE_OSPACE
      G4std::sort(_levels->begin(), _levels->end());
#else
      G4std::stable_sort(_levels->begin(), _levels->end());
#endif
    }
  else 
    {
      _levels = 0;
    }
}









