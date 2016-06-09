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

#include "G4NuclearLevelStore.hh"
#include "g4std/strstream"

G4std::map<G4String,G4NuclearLevelManager*> G4NuclearLevelStore::theManagers;
G4String G4NuclearLevelStore::dirName("");

G4NuclearLevelStore* G4NuclearLevelStore::GetInstance()
{
  static G4NuclearLevelStore theInstance;
  return &theInstance;
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
    G4std::map<G4String,G4NuclearLevelManager*>::iterator i;
    for (i = theManagers.begin(); i != theManagers.end(); ++i)
    {
	if ( (*i).second ) delete (*i).second;
    }
}

G4String G4NuclearLevelStore::GenerateKey(const G4int Z, const G4int A)
{
    char chname[10] = {' '};
    G4std::ostrstream streamName(chname, 10, G4std::ios::out);
    streamName << 'z' << Z << ".a" << A;
    G4String name(chname);
    return name;
}


G4NuclearLevelManager* G4NuclearLevelStore::GetManager(const G4int Z, const G4int A)
{
    G4NuclearLevelManager * result = 0; 
    if (A < 1 || Z < 1 || A < Z)
    {
	G4cerr << "G4NuclearLevelStore::GetManager: Wrong values Z = " << Z 
	       << " A = " << A << '\n';
	return result;
    }
    // Generate the key = filename
    G4String key(this->GenerateKey(Z,A));
    
    // Check if already exists that key
    G4std::map<G4String,G4NuclearLevelManager*>::iterator idx = theManagers.find(key);
    // If doesn't exists then create it
    if ( idx == theManagers.end() )
    {
	result = new G4NuclearLevelManager();
	result->SetNucleus(Z,A,dirName + key);
	theManagers.insert(G4std::make_pair(key,result));
    }
    // But if it exists...
    else
    {
	result = idx->second;
    }
    
    return result; 
}
