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
    _nucleusA(0), _nucleusZ(0), _fileName(""), _validity(false), 
    _levels(0), _levelEnergy(0), _gammaEnergy(0), _probability(0)
{ }

G4NuclearLevelManager::G4NuclearLevelManager(const G4int Z, const G4int A, const G4String& filename) :
    _nucleusA(A), _nucleusZ(Z), _fileName(filename)
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
	G4std::for_each(_levels->begin(), _levels->end(), DeleteLevel());
	_levels->clear();
      }
    delete _levels;
    _levels = 0;
  }
}

void G4NuclearLevelManager::SetNucleus(const G4int Z, const G4int A, const G4String& filename)
{
    if (_nucleusZ != Z || _nucleusA != A)
    {
	_nucleusA = A;
	_nucleusZ = Z;
	_fileName = filename;
	MakeLevels();
    }
    
}

G4bool G4NuclearLevelManager::IsValid() const
{
    return _validity;
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


const G4NuclearLevel* G4NuclearLevelManager::
NearestLevel(const G4double energy, const G4double eDiffMax) const
{
    G4int iNear = -1;
  
    G4double diff = 9999. * GeV;
    if (_levels != 0)
    {
	unsigned int i = 0;
	for (i=0; i<_levels->size(); i++)
	{
	    G4double e = _levels->operator[](i)->Energy();
	    G4double eDiff = abs(e - energy);
	    if (eDiff < diff && eDiff <= eDiffMax)
	    { 
		diff = eDiff; 
		iNear = i;
	    }
	}
    }
    if (_levels != 0 && iNear >= 0 && iNear < static_cast<G4int>(_levels->size()) )
    { 
	return _levels->operator[](iNear); 
    }
    else
    { 
	return 0; 
    }
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
	       >> _angularMomentum  >> _totalCC >> _kCC >> _l1CC >> _l2CC 
	       >> _l3CC >> _m1CC >> _m2CC >> _m3CC >> _m4CC >> _m5CC
	       >> _nPlusCC;
      _levelEnergy *= keV;
      _gammaEnergy *= keV;
      _halfLife *= second;
	
      // The following adjustment is needed to take care of anomalies in 
      // data files, where some transitions show up with relative probability
      // zero
      if (_probability < minProbability) _probability = minProbability;
      // the folowwing is to convert icc probability to accumulative ones
      _l1CC += _kCC;
      _l2CC += _l1CC;
      _l3CC += _l2CC;
      _m1CC += _l3CC;
      _m2CC += _m1CC;
      _m3CC += _m2CC;
      _m4CC += _m3CC;
      _m5CC += _m4CC;
      _nPlusCC += _m5CC;
      _kCC /= _nPlusCC;
      _l1CC /= _nPlusCC;
      _l2CC /= _nPlusCC;
      _l3CC /= _nPlusCC;
      _m1CC /= _nPlusCC;
      _m2CC /= _nPlusCC;
      _m3CC /= _nPlusCC;
      _m4CC /= _nPlusCC;
      _m5CC /= _nPlusCC;
      _nPlusCC /= _nPlusCC;  
	
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
  _validity = false;
#ifdef G4USE_STD_NAMESPACE  
  G4std::ifstream inFile(_fileName, G4std::ios::in);
#else
  ifstream inFile(_fileName, ios::in|ios::nocreate);
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
  else 
    {
      _validity = true;
    }

  _levels = new G4PtrLevelVector;
  
  G4DataVector eLevel;
  G4DataVector eGamma;
  G4DataVector wGamma;
  G4DataVector pGamma; // polarity
  G4DataVector hLevel; // half life
  G4DataVector aLevel; // angular momentum
  G4DataVector kConve; //  internal convertion coefficiencies
  G4DataVector l1Conve;
  G4DataVector l2Conve;
  G4DataVector l3Conve;
  G4DataVector m1Conve;
  G4DataVector m2Conve;
  G4DataVector m3Conve;
  G4DataVector m4Conve;
  G4DataVector m5Conve;
  G4DataVector npConve;
  G4DataVector toConve;
 
	
  while (Read(inFile))
    {
      eLevel.push_back(_levelEnergy);
      eGamma.push_back(_gammaEnergy);
      wGamma.push_back(_probability);
      pGamma.push_back(_polarity);
      hLevel.push_back(_halfLife);
      aLevel.push_back(_angularMomentum);
      kConve.push_back(_kCC);
      l1Conve.push_back(_l1CC);
      l2Conve.push_back(_l2CC);
      l3Conve.push_back(_l3CC);
      m1Conve.push_back(_m1CC);
      m2Conve.push_back(_m2CC);
      m3Conve.push_back(_m3CC);
      m4Conve.push_back(_m4CC);
      m5Conve.push_back(_m5CC);
      npConve.push_back(_nPlusCC);
      toConve.push_back(_totalCC);
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
  G4DataVector thisLevelkCC;
  G4DataVector thisLevell1CC;
  G4DataVector thisLevell2CC;
  G4DataVector thisLevell3CC;
  G4DataVector thisLevelm1CC;
  G4DataVector thisLevelm2CC;
  G4DataVector thisLevelm3CC;
  G4DataVector thisLevelm4CC;
  G4DataVector thisLevelm5CC;
  G4DataVector thisLevelnpCC;
  G4DataVector thisLeveltoCC;
 
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
							thisLevelPolarities,
							thisLevelkCC,
							thisLevell1CC,
							thisLevell2CC,
							thisLevell3CC,
							thisLevelm1CC,
							thisLevelm2CC,
							thisLevelm3CC,
							thisLevelm4CC,
							thisLevelm5CC,
							thisLevelnpCC,
							thisLeveltoCC );
							
	  _levels->push_back(newLevel);
	  // Reset data vectors
	  thisLevelEnergies.clear();
	  thisLevelWeights.clear();
	  thisLevelPolarities.clear();
	  thisLevelkCC.clear();
	  thisLevell1CC.clear();
	  thisLevell2CC.clear();
	  thisLevell3CC.clear();
	  thisLevelm1CC.clear();
	  thisLevelm2CC.clear();
	  thisLevelm3CC.clear();
	  thisLevelm4CC.clear();
	  thisLevelm5CC.clear();
	  thisLevelnpCC.clear();
	  thisLeveltoCC.clear();
	  thisLevelEnergy = e;
	}
      // Append current data
      thisLevelEnergies.push_back(eGamma[i]);
      thisLevelWeights.push_back(wGamma[i]);
      thisLevelPolarities.push_back(pGamma[i]);
      thisLevelkCC.push_back(kConve[i]);
      thisLevell1CC.push_back(l1Conve[i]);
      thisLevell2CC.push_back(l2Conve[i]);
      thisLevell3CC.push_back(l3Conve[i]);
      thisLevelm1CC.push_back(m1Conve[i]);
      thisLevelm2CC.push_back(m2Conve[i]);
      thisLevelm3CC.push_back(m3Conve[i]);
      thisLevelm4CC.push_back(m4Conve[i]);
      thisLevelm5CC.push_back(m5Conve[i]);
      thisLevelnpCC.push_back(npConve[i]);
      thisLeveltoCC.push_back(toConve[i]);
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
						      thisLevelPolarities,
						      thisLevelkCC,
						      thisLevell1CC,
						      thisLevell2CC,
						      thisLevell3CC,
						      thisLevelm1CC,
						      thisLevelm2CC,
						      thisLevelm3CC,
						      thisLevelm4CC,
						      thisLevelm5CC,
						      thisLevelnpCC,
						      thisLeveltoCC );

	_levels->push_back(newLevel);
    }

    //    G4std::sort(_levels->begin(), _levels->end());

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
    { _levels->operator[](i)->PrintAll(); }
}


G4NuclearLevelManager::G4NuclearLevelManager(const G4NuclearLevelManager &right)
{
    _levelEnergy = right._levelEnergy;
    _gammaEnergy = right._gammaEnergy;
    _probability = right._probability;
    _polarity = right._polarity;
    _halfLife = right._halfLife;
    _angularMomentum = right._angularMomentum;
    _kCC = right._kCC;
    _l1CC = right._l1CC;
    _l2CC = right._l2CC;
    _l3CC = right._l3CC;
    _m1CC = right._m1CC;
    _m2CC = right._m2CC;
    _m3CC = right._m3CC;
    _m4CC = right._m4CC;
    _m5CC = right._m5CC;
    _nPlusCC = right._nPlusCC;
    _totalCC = right._totalCC;
    _nucleusA = right._nucleusA;
    _nucleusZ = right._nucleusZ;
    _fileName = right._fileName;
    _validity = right._validity;
    if (right._levels != 0)   
      {
	_levels = new G4PtrLevelVector;
	G4int n = right._levels->size();
	G4int i;
	for (i=0; i<n; i++)
	  {
	    _levels->push_back(new G4NuclearLevel(*(right._levels->operator[](i))));
	  }
	
	G4std::sort(_levels->begin(), _levels->end());
      }
    else 
      {
	_levels = 0;
      }
}

