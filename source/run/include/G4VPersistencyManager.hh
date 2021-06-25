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
// G4VPersistencyManager
//
// Class description:
//
// This is an abstract base class for persistency management. The user's
// concrete class derived from this class must act as a singleton. The user
// must construct the object of his/her concrete persistency manager in
// the main().
// The virtual methods of Store() and Retrieve() will be invoked from
// G4RunManager if the persistency manager exists.
// Even if the user does not use any ODBMS system, the user can use this
// class especially for Store() methods. Writing an ASCII file for storing
// event information can be delegated to this class, for example.

// Author: Youhei Morita, 2001
// --------------------------------------------------------------------
#ifndef G4VPersistencyManager_hh
#define G4VPersistencyManager_hh 1

#include "globals.hh"

class G4Event;
class G4Run;
class G4VPhysicalVolume;

class G4VPersistencyManager
{
  public:

    static G4VPersistencyManager* GetPersistencyManager();
      // Static method to return the pointer to the singleton object.
      // Note that this method does NOT create the singleton itself.

    virtual ~G4VPersistencyManager();

    virtual G4bool Store(const G4Event* anEvent)          = 0;
    virtual G4bool Store(const G4Run* aRun)               = 0;
    virtual G4bool Store(const G4VPhysicalVolume* world)  = 0;
      // Stores G4Event, G4Run, and geometry tree characterised
      // by the world volume.

    virtual G4bool Retrieve(G4Event*& anEvent)            = 0;
    virtual G4bool Retrieve(G4Run*& aRun)                 = 0;
    virtual G4bool Retrieve(G4VPhysicalVolume*& theWorld) = 0;
      // Restores G4Event, G4Run, and geometry tree characterised
      // by the world volume.

  protected:

    G4VPersistencyManager();

  private:

    static G4ThreadLocal G4VPersistencyManager* fPersistencyManager;
};

#endif
