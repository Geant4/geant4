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
//      
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVELMANAGER_HH
#define G4NUCLEARLEVELMANAGER_HH

#include "globals.hh"
#include "G4PtrLevelVector.hh"
#include "G4NuclearLevel.hh"
#include "G4ios.hh"
#include <fstream>

class G4NuclearLevelManager 
{

public:

  G4NuclearLevelManager();
  G4NuclearLevelManager(const G4int Z, const G4int A, const G4String& filename);
  G4NuclearLevelManager(const G4NuclearLevelManager & right);  

  ~G4NuclearLevelManager();
  
  void SetNucleus(const G4int Z, const G4int A, const G4String& filename);

  G4bool IsValid() const;

  G4int NumberOfLevels() const;

  const G4PtrLevelVector* GetLevels() const;

  const G4NuclearLevel* NearestLevel(const G4double energy, const G4double eDiffMax=9999.*GeV) const;

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
 
  void MakeLevels();

  G4int _nucleusA;
  G4int _nucleusZ;
  G4String _fileName;
  G4bool _validity;
  G4PtrLevelVector* _levels;
  
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
