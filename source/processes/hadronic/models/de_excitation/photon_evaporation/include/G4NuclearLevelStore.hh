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
// 04-10-2010  M. Kelsey -- Replace G4String keys with integers (ZZZAAA),
//		            move string operation to GenerateFilename()
// 17-11-2010 V. Ivanchenko - make as a classical singleton. 
//
// 17-10-2011 L. Desorgher -- Allow the user to replace the radioactive 
//               decay data provided in Geant4 by its own data file for 
//               a given isotope
// 06-01-2012 V. Ivanchenko - added nuclear level data structure for 
//               usage in HEP where internal conversion is neglected;
//               cleanup the code; new method GetLevelManager  

#ifndef G4NuclearLevelStore_hh 
#define G4NuclearLevelStore_hh 1

#include "G4NuclearLevelManager.hh"
#include "G4LevelManager.hh"
#include "G4LevelReader.hh"
#include "globals.hh"
#include <map>

class G4NuclearLevelStore
{
private:

  G4NuclearLevelStore();

public:

  static G4NuclearLevelStore* GetInstance();

  G4NuclearLevelManager* GetManager(G4int Z, G4int A);

  G4LevelManager* GetLevelManager(G4int Z, G4int A);

  ~G4NuclearLevelStore();

  void AddUserEvaporationDataFile(G4int Z, G4int A, 
				  const G4String& filename);

private:

  G4int GenerateKey(G4int Z, G4int A) const { return Z*1000+A; }

  G4String GenerateFileName(G4int Z, G4int A) const;

  typedef std::map<G4int,G4NuclearLevelManager*> ManagersMap;
  typedef std::map<G4int,G4LevelManager*> MapForHEP;

  G4LevelReader reader;
  ManagersMap   theManagers;
  MapForHEP     managersForHEP;
  G4String      dirName;

  static G4NuclearLevelStore* theInstance;

  G4bool userFiles;
  std::map<G4int, G4String> theUserDataFiles;
};
#endif
