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
// G4VDigitizerModule
//
// class description:
//
// This is the abstract base class of the digitizer module. The user's
// digitizer module which generates digits must be derived from this
// class.
// In the derived class constructor, name(s) of digi collection(s) which
// are made by the digitizer module must be set to "collectionName" string
// vector.

// Author: M.Asai
// --------------------------------------------------------------------
#ifndef G4VDigitizerModule_hh
#define G4VDigitizerModule_hh 1

#include "globals.hh"
#include <vector>

class G4DigiManager;
class G4VDigiCollection;

class G4VDigitizerModule
{
  public: // with description

    G4VDigitizerModule(const G4String& modName);
    // Constructor. The user's concrete class must use this constructor
    // by the constructor initializer of the derived class. The name of
    // the detector module must be unique.

    virtual ~G4VDigitizerModule();
    G4bool operator==(const G4VDigitizerModule& right) const;
    G4bool operator!=(const G4VDigitizerModule& right) const;

    virtual void Digitize() = 0;
    //  The pure virtual method that the derived class must implement.
    // In the concrete implementation of this method, necessary digi
    // collection object must be constructed and set to G4DCofThisEvent
    // by StoreDigiCollection protected method.

  public:

    inline G4int GetNumberOfCollections() const
    { return G4int(collectionName.size()); }
    inline G4String GetCollectionName(G4int i) const
    { return collectionName[i]; }
    inline G4String GetName() const
    { return moduleName; }
    inline void SetVerboseLevel(G4int val)
    { verboseLevel = val; }

  protected:

    void StoreDigiCollection(G4VDigiCollection* aDC);
    void StoreDigiCollection(G4int DCID,G4VDigiCollection* aDC);

  protected:

    G4DigiManager* DigiManager;
    G4String moduleName;
    std::vector<G4String> collectionName;
    G4int verboseLevel;
};

#endif
