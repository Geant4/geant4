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
// $Id: G4NuclearLevelData.cc 86536 2014-11-13 19:05:21Z vnivanch $
//
// -------------------------------------------------------------------
//
//      GEANT4 source file 
//
//      File name:     G4NuclearLevelData
//
//      Author:        V.Ivanchenko
// 
//      Creation date: 10 February 2015
//
//      Modifications:
//      
// -------------------------------------------------------------------

#include "G4NuclearLevelData.hh"
#include "G4LevelReader.hh"
#include "G4LevelManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

G4NuclearLevelData* G4NuclearLevelData::theInstance = 0;

const G4int G4NuclearLevelData::AMIN[] = {0,
    1,3,6,7,8,10,12,14,17,18,               // Z= 1-10
   21,22,23,26,27,28,32,32,36,36,           // Z= 11-20
   41,42,45,45,48,49,53,54,57,59,           // Z= 21-30
   61,62,66,68,70,72,74,76,78,80,           // Z= 31-40
   82,84,86,88,92,93,95,98,101,101,         // Z= 41-50
  105,105,109,110,114,118,121,122,125,128,  // Z= 51-60
  131,130,137,138,140,140,141,144,145,151,  // Z= 61-70
  151,154,159,160,161,162,168,168,172,172,  // Z= 71-80
  181,180,187,189,196,196,206,206,212,214,  // Z= 81-90
  229,230,233,236,241,240,240,245,243,249,  // Z= 91-100
  251,251};
const G4int G4NuclearLevelData::AMAX[] = {0,
    3,8,9,12,13,16,21,22,27,32,             // Z= 1-10
   32,36,35,42,43,46,45,48,50,53,           // Z= 11-20
   56,58,60,64,64,68,68,76,73,78,           // Z= 21-30
   81,84,84,88,88,96,96,102,102,108,        // Z= 31-40
  105,110,111,114,115,121,120,130,130,134,  // Z= 41-50
  135,139,139,144,145,148,149,152,151,156,  // Z= 51-60
  155,160,159,164,166,170,169,172,175,178,  // Z= 61-70
  180,184,190,190,190,198,198,204,201,208,  // Z= 71-80,
  210,212,215,218,217,222,227,232,232,234,  // Z= 81-90,
  236,240,242,246,246,249,251,253,254,256,  // Z= 91-100
  251,254};

#ifdef G4MULTITHREADED
G4Mutex G4NuclearLevelData::nuclearLevelDataMutex = G4MUTEX_INITIALIZER;
#endif

G4NuclearLevelData* G4NuclearLevelData::GetInstance()
{
  if (!theInstance)  { 
    static G4NuclearLevelData theData;
    theInstance = &theData; 
  }
  return theInstance;
}   

G4NuclearLevelData::G4NuclearLevelData()
{
  fLevelReader = new G4LevelReader();
  for(G4int Z=0; Z<ZMAX; ++Z) {
    (fLevelManagers[Z]).resize(AMAX[Z]-AMIN[Z]+1,nullptr);
    (fLevelManagerFlags[Z]).resize(AMAX[Z]-AMIN[Z]+1,false);
  }
}

G4NuclearLevelData::~G4NuclearLevelData()
{
  delete fLevelReader;
  for(G4int Z=1; Z<ZMAX; ++Z) {
    size_t nn = (fLevelManagers[Z]).size();
    for(size_t j=0; j<nn; ++j) { delete (fLevelManagers[Z])[j]; }
  }
}

const G4LevelManager* 
G4NuclearLevelData::GetLevelManager(G4int Z, G4int A)
{
  const G4LevelManager* man = nullptr;
  //G4cout << "G4NuclearLevelData: Z= " << Z << " A= " << A << G4endl;  
  if(0 < Z && Z < ZMAX && A >= AMIN[Z] && A <= AMAX[Z]) {
    if(!(fLevelManagerFlags[Z])[A - AMIN[Z]]) {
      InitialiseForIsotope(Z, A);
    }
    man = (fLevelManagers[Z])[A - AMIN[Z]];
  }
  //G4cout << man << G4endl;
  return man;
}

G4bool 
G4NuclearLevelData::AddPrivateData(G4int Z, G4int A, const G4String& filename)
{
  G4bool res = false; 
  if(A >= AMIN[Z] && A <= AMAX[Z]) { 
    const G4LevelManager* newman = 
      fLevelReader->MakeLevelManager(Z, A, filename);
    if(newman) { 
      delete (fLevelManagers[Z])[A - AMIN[Z]]; 
      (fLevelManagers[Z])[A - AMIN[Z]] = newman;
      (fLevelManagerFlags[Z])[A - AMIN[Z]] = true;
      res = true;
    }
  }
  return res;
}

void G4NuclearLevelData::InitialiseForIsotope(G4int Z, G4int A)
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4NuclearLevelData::nuclearLevelDataMutex);
#endif
  if(!(fLevelManagerFlags[Z])[A - AMIN[Z]]) {
    (fLevelManagers[Z])[A - AMIN[Z]] = 
      fLevelReader->CreateLevelManager(Z, A);
    (fLevelManagerFlags[Z])[A - AMIN[Z]] = true;
  }
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4NuclearLevelData::nuclearLevelDataMutex);
#endif
}
