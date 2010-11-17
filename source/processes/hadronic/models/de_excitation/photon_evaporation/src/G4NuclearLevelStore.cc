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
// $Id: G4NuclearLevelStore.cc,v 1.5 2010-11-17 16:50:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 06-10-2010 M. Kelsey -- Drop static data members.
// 17-11-2010 V. Ivanchenko - make as a classical singleton. 

#include "G4NuclearLevelStore.hh"
#include <sstream>

G4NuclearLevelStore* G4NuclearLevelStore::theInstance = 0;

G4NuclearLevelStore* G4NuclearLevelStore::GetInstance()
{
  if(!theInstance) {
    static G4NuclearLevelStore store;
    theInstance = &store;
  }
  return theInstance;
}

G4NuclearLevelStore::G4NuclearLevelStore()
{
  char* env = getenv("G4LEVELGAMMADATA");
  if (env == 0) 
    {
      G4cout << "G4NuclarLevelStore: please set the G4LEVELGAMMADATA environment variable\n";
      dirName = "";
    }
  else
    {
      dirName = env;
      dirName += '/';
    }
}

G4NuclearLevelStore::~G4NuclearLevelStore()
{
  ManagersMap::iterator i;
  for (i = theManagers.begin(); i != theManagers.end(); ++i)
    delete i->second;
}

G4String 
G4NuclearLevelStore::GenerateFilename(G4int Z, G4int A) const 
{
  std::ostringstream streamName; 
  streamName << 'z' << Z << ".a" << A;
  G4String name(streamName.str());
  return name;
}

G4NuclearLevelManager* 
G4NuclearLevelStore::GetManager(G4int Z, G4int A) 
{
  G4NuclearLevelManager * result = 0; 
  if (A < 1 || Z < 1 || A < Z)
    {
	G4cerr << "G4NuclearLevelStore::GetManager: Wrong values Z = " << Z 
	       << " A = " << A << '\n';
	return result;
    }

  // Generate the key = filename
  G4int key = GenerateKey(Z,A);
    
  // Check if already exists that key
  ManagersMap::iterator idx = theManagers.find(key);
  // If doesn't exists then create it
  if ( idx == theManagers.end() )
    {
      G4String file = GenerateFilename(Z,A);
      result = new G4NuclearLevelManager(Z,A,dirName + file);
      theManagers.insert(std::make_pair(key,result));
    }
  // But if it exists...
  else
    {
      result = idx->second;
    }
    
  return result; 
}
