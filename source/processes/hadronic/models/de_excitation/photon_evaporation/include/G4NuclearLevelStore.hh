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
// $Id: G4NuclearLevelStore.hh,v 1.4 2010-11-17 16:50:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 04-10-2010  M. Kelsey -- Replace G4String keys with integers (ZZZAAA),
//		            move string operation to GenerateFilename()
// 17-11-2010 V. Ivanchenko - make as a classical singleton. 

#ifndef G4NuclearLevelStore_hh 
#define G4NuclearLevelStore_hh 1

#include "G4NuclearLevelManager.hh"
#include "G4String.hh"
#include <map>

class G4NuclearLevelStore
{
private:
  G4NuclearLevelStore();

public:

  static G4NuclearLevelStore* GetInstance();

  G4NuclearLevelManager* GetManager(G4int Z, G4int A);
  ~G4NuclearLevelStore();

private:
  G4int GenerateKey(G4int Z, G4int A) const { return Z*1000+A; }

  G4String GenerateFilename(G4int Z, G4int A) const;

  typedef std::map<G4int,G4NuclearLevelManager*> ManagersMap;

  ManagersMap theManagers;
  G4String dirName;

  static G4NuclearLevelStore* theInstance;
};
#endif
