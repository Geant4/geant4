// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4NuclearLevelManager
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 25 October 1998
//
//      Modifications: 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added half-life, angular momentum, parity, emissioni type
//              reading from experimental data. 
//      
// -------------------------------------------------------------------

#ifndef G4NUCLEARLEVELMANAGER_HH
#define G4NUCLEARLEVELMANAGER_HH

#include "globals.hh"
#include "G4PtrLevelVector.hh"
#include "G4NuclearLevel.hh"
#include "G4ios.hh"
#include <fstream.h>

class G4NuclearLevelManager 
{

public:

  G4NuclearLevelManager();
  G4NuclearLevelManager(G4int Z, G4int A);

  ~G4NuclearLevelManager();
  
  void SetNucleus(G4int Z, G4int A);

  G4bool IsValid(G4int Z, G4int A) const;

  G4int NumberOfLevels() const;

  const G4PtrLevelVector* GetLevels() const;

  const G4NuclearLevel* NearestLevel(G4double energy, G4double eDiffMax=9999.*GeV) const;

  const G4NuclearLevel* LowestLevel() const;
  const G4NuclearLevel* HighestLevel() const;

  G4double MinLevelEnergy() const;
  G4double MaxLevelEnergy() const;

  void PrintAll();

  G4NuclearLevelManager(const G4NuclearLevelManager &right);  
  
protected:

private:  

  const G4NuclearLevelManager& operator=(const G4NuclearLevelManager &right);
  G4bool operator==(const G4NuclearLevelManager &right) const;
  G4bool operator!=(const G4NuclearLevelManager &right) const;

  G4bool Read(ifstream& aDataFile);
 
  void MakeLevels();

  G4int _nucleusA;
  G4int _nucleusZ;
  G4PtrLevelVector* _levels;
  
  G4double _levelEnergy;
  G4double _gammaEnergy;
  G4double _probability;
  G4double _polarity;
  G4double _halfLife;
  G4double _angularMomentum;
};

#endif














