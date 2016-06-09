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
// $Id$
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
//      Creation date: 25 October 1998
//
//      Modifications: 
//        21 Nov. 2001, Fan Lei (flei@space.qinetiq.com)
//              Added K->N+ internal  conversion coefficiencies and their access
//              functions      
//      
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//        02 May 2003,   Vladimir Ivanchenko remove rublic copy contructor
//        06 Oct 2010, M. Kelsey -- Use object storage, not pointers, drop
//		public access to list
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVELMANAGER_HH
#define G4NUCLEARLEVELMANAGER_HH 1

#include <iosfwd>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4PtrLevelVector.hh"

class G4NuclearLevel;


class G4NuclearLevelManager 
{

public:

  //G4NuclearLevelManager();
  G4NuclearLevelManager(G4int Z, G4int A, const G4String& filename);
  G4NuclearLevelManager(const G4NuclearLevelManager & right);  

  ~G4NuclearLevelManager();
  
  void SetNucleus(G4int Z, G4int A, const G4String& filename);

  inline G4bool IsValid() const { return _validity; }

  inline G4int NumberOfLevels() const { return (_levels ? _levels->size() : 0); }

  const G4NuclearLevel* GetLevel(G4int i) const;

  const G4NuclearLevel* NearestLevel(G4double energy,
				     G4double eDiffMax=9999.*CLHEP::GeV) const;

  const G4NuclearLevel* LowestLevel() const;
  const G4NuclearLevel* HighestLevel() const;

  G4double MinLevelEnergy() const;
  G4double MaxLevelEnergy() const;

  void PrintAll();
  
private:  
  const G4NuclearLevelManager& operator=(const G4NuclearLevelManager &right);
  G4bool operator==(const G4NuclearLevelManager &right) const;
  G4bool operator!=(const G4NuclearLevelManager &right) const;

  G4bool Read(std::ifstream& aDataFile);
  G4bool ReadDataLine(std::ifstream& dataFile);
  G4bool ReadDataItem(std::istream& dataFile, G4double& x);
  void ProcessDataLine();

  void MakeLevels();
  void ClearLevels();

  G4NuclearLevel* UseLevelOrMakeNew(G4NuclearLevel* level);
  void AddDataToLevel(G4NuclearLevel* level);
  void FinishLevel(G4NuclearLevel* level);

  G4int _nucleusA;
  G4int _nucleusZ;
  G4String _fileName;
  G4bool _validity;
  G4PtrLevelVector* _levels;

  // Buffers for reading data file
  char buffer[30];		// For doubles in scientific notation
  G4double _levelEnergy;
  G4double _gammaEnergy;
  G4double _probability;
  G4double _polarity;
  G4double _halfLife;
  G4double _angularMomentum;
  G4double _kCC;
  G4double _l1CC;
  G4double _l2CC;
  G4double _l3CC;
  G4double _m1CC;  
  G4double _m2CC;
  G4double _m3CC;
  G4double _m4CC;
  G4double _m5CC;
  G4double _nPlusCC;
  G4double _totalCC;
};

#endif
