// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigitizerModule.hh,v 1.1 1999-01-07 16:14:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VDigitizerModule_H
#define G4VDigitizerModule_H 1

class G4DigiManager;
class G4VDigiCollection;
#include "globals.hh"
#include <rw/tvordvec.h>

class G4VDigitizerModule
{
  public:
    G4VDigitizerModule(G4String modName);
    virtual ~G4VDigitizerModule();
    int operator==(const G4VDigitizerModule &right) const;
    int operator!=(const G4VDigitizerModule &right) const;

  public:
    virtual void Digitize() = 0;

  protected:
    void StoreDigiCollection(G4VDigiCollection* aDC);
    void StoreDigiCollection(G4int DCID,G4VDigiCollection* aDC);

  protected:
    G4DigiManager* DigiManager;
    G4String moduleName;
    RWTValOrderedVector<G4String> collectionName;
    G4int verboseLevel;

  public:
    inline G4int GetNumberOfCollections() const
    { return collectionName.entries(); }
    inline G4String GetCollectionName(G4int i) const
    { return collectionName[i]; }
    inline G4String GetName() const
    { return moduleName; }
    inline void SetVerboseLevel(G4int val)
    { verboseLevel = val; }
};

#endif

