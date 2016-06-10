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
//
// $Id: G4VPersistencyManager.hh 66892 2013-01-17 10:57:59Z gunter $
//

#ifndef G4VPersistencyManager_h
#define G4VPersistencyManager_h 1

#include "globals.hh"

class G4Event;
class G4Run;
class G4VPhysicalVolume;

// class description:
//
//  This is an abstract base class for persistency management. The user's
// concrete class derived from this class must be a singleton. The user
// must construct the object of his/her concrete persistency manager at
// his/her main().
//  The virtual methods of Store() and Retreive() will be invoked from
// G4RunManager if the persistency manager exists.
//  Even if the user does not use any ODBMS, the user can use this class
// especially for Store() methods. Writing an ASCII file for storing
// event information can be delegated to this class, for example.
//

class G4VPersistencyManager 
{
  public: // with description
      static G4VPersistencyManager* GetPersistencyManager();
      //  Static method to return the pointer to the singleton object.
      // Note that this method does NOT create the singleton object.

  protected:
      G4VPersistencyManager();

  public:
      virtual ~G4VPersistencyManager();

  private: 
      static G4ThreadLocal G4VPersistencyManager * fPersistencyManager;

  public: // with description
      virtual G4bool Store(const G4Event* anEvent)=0;
      virtual G4bool Store(const G4Run* aRun)=0;
      virtual G4bool Store(const G4VPhysicalVolume* theWorld)=0;
      //  Stores G4Event, G4Run, and geometry tree characterized by the world volume.

      virtual G4bool Retrieve(G4Event*& anEvent)=0;
      virtual G4bool Retrieve(G4Run*& aRun)=0;
      virtual G4bool Retrieve(G4VPhysicalVolume*& theWorld)=0;
      //  Restore G4Event, G4Run, and geometry tree characterized by the world volume.

};




#endif

