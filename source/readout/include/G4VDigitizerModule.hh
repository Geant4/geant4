// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigitizerModule.hh,v 1.4 1999-12-15 14:53:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4VDigitizerModule_H
#define G4VDigitizerModule_H 1

class G4DigiManager;
class G4VDigiCollection;
#include "globals.hh"
#include "g4rw/tvordvec.h"

// class description:
//
//  This is the abstract base class of the digitizer module. The user's
// digitizer module which generates digits must be derived from this
// class.
//  In the derived class constructor, name(s) of digi collection(s) which
// are made by the digitizer module must be set to "collectionName" string
// vector.

class G4VDigitizerModule
{
  public: // with description
    G4VDigitizerModule(G4String modName);
    //  Constructor. The user's concrete class must use this constructor
    // by the constructor initializer of the derived class. The name of
    // the detector module must be unique.
  public:
    virtual ~G4VDigitizerModule();
    int operator==(const G4VDigitizerModule &right) const;
    int operator!=(const G4VDigitizerModule &right) const;

  public: // with description
    virtual void Digitize() = 0;
    //  The pure virtual method that the derived class must implement.
    // In the concrete implementation of this method, necessary digi
    // collection object must be constructed and set to G4DCofThisEvent
    // by StoreDigiCollection protected method.

  protected:
    void StoreDigiCollection(G4VDigiCollection* aDC);
    void StoreDigiCollection(G4int DCID,G4VDigiCollection* aDC);

  protected:
    G4DigiManager* DigiManager;
    G4String moduleName;
    G4RWTValOrderedVector<G4String> collectionName;
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

