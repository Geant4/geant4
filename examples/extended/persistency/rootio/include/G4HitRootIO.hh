// $Id: G4HitRootIO.hh,v 1.1 2002-12-04 02:44:28 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4HitRootIO.hh
//
// History:
//   '02.5.7  Youhei Morita  Initial creation

#ifndef HIT_ROOT_IO_HH
#define HIT_ROOT_IO_HH 1

#include "G4HCofThisEvent.hh"
#include "G4PersistencyCenter.hh"
#include "G4RootIOManager.hh"
#include "G4HCIOcatalog.hh"

// Class inherited:
#include "G4VPHitIO.hh"

// Class Description:
//   Manager class to store and retrieve Hit objects.
// 
//   This is a singleton class and should be constructed only
//   by GetG4HitRootIO().

class G4HitRootIO
 : public G4VPHitIO
{
    public: // With description
      G4HitRootIO();
      // Constructor

      virtual ~G4HitRootIO();
      // Destructor

    public: // With description
      static G4HitRootIO* GetG4HitRootIO();
      // Construct a new singleton G4HitRootIO object if it does not exist.

      bool Store(const G4HCofThisEvent* hcevt);
      // Store hit collections.
      // Concrete class of G4VPHitsCollectionIO must be registered
      // with G4VPHitIO::AddHCIOmanager() before calling this method.

      bool Retrieve(G4HCofThisEvent*& hcevt);
      // Retrieve hit collections
      // Concrete class of G4VPHitsCollectionIO must be registered
      // with G4VPHitIO::AddHCIOmanager() before calling this method.

}; // End of class G4HitRootIO

#endif

