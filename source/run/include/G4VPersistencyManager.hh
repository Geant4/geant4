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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VPersistencyManager.hh,v 1.4 2001-07-11 10:08:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
      static G4VPersistencyManager * fPersistencyManager;

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

