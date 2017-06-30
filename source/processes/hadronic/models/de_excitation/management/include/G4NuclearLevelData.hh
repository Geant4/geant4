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
// $Id: G4NuclearLevelData.hh 86536 2014-11-13 19:05:21Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 header file 
//
//      File name:     G4NuclearLevelData
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 9 February 2014
//
//      Modifications:
//      
// -------------------------------------------------------------------
//
// Nuclear level data uploaded at initialisation of Geant4 from 
// data files of the G4LEVELGAMMADATA
// 

#ifndef G4NUCLEARLEVELDATA_HH
#define G4NUCLEARLEVELDATA_HH 1

#include "globals.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4Threading.hh"
#include <vector>

class G4LevelReader;
class G4LevelManager;

class G4NuclearLevelData 
{
private:

  G4NuclearLevelData();

  static G4NuclearLevelData* theInstance;

public:

  static G4NuclearLevelData* GetInstance();

  ~G4NuclearLevelData();

  // run time call to access or to create level manager
  const G4LevelManager* GetLevelManager(G4int Z, G4int A);

  // add private data to isotope from master thread
  G4bool AddPrivateData(G4int Z, G4int A, const G4String& filename);

  // access to min/max A in the level DB 
  G4int GetMinA(G4int Z) const;
  G4int GetMaxA(G4int Z) const;

  // check max energy of a level without upload of the data
  G4double GetMaxLevelEnergy(G4int Z, G4int A) const;
  G4float MaxLevelEnergy(G4int Z, G4int A) const;

  // check closest level if the energy is below the max level energy
  G4double GetLevelEnergy(G4int Z, G4int A, G4double energy);

  // check closest level below given energy 
  G4double GetLowEdgeLevelEnergy(G4int Z, G4int A, G4double energy);

  // check if residual excitation energy corresponding to
  // discrete level and if it is the case select closest level
  G4double FindLevel(G4int Z, G4int A, G4double resMass, G4double Mass,
                     G4double partMass, G4double T);

  // access to all model parameters
  G4DeexPrecoParameters* GetParameters();
  
private:

  void InitialiseForIsotope(G4int Z, G4int A);

  G4NuclearLevelData(G4NuclearLevelData &) = delete;
  G4NuclearLevelData & operator=(const G4NuclearLevelData &right) = delete;

  G4DeexPrecoParameters* fDeexPrecoParameters;
  G4LevelReader*    fLevelReader;

  static const G4int ZMAX = 118;
  static const G4int AMIN[ZMAX];
  static const G4int AMAX[ZMAX];
  static const G4int LEVELIDX[ZMAX];

  std::vector<const G4LevelManager*> fLevelManagers[ZMAX];
  std::vector<G4bool> fLevelManagerFlags[ZMAX];

#ifdef G4MULTITHREADED
  static G4Mutex nuclearLevelDataMutex;
#endif
};

#endif
