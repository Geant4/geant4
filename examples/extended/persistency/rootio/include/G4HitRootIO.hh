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
// $Id: G4HitRootIO.hh,v 1.4 2002-12-13 14:45:41 gunter Exp $
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
//   by GetHitRootIO().

class G4HitRootIO
 : public G4VPHitIO
{
    public: // With description
      G4HitRootIO();
      // Constructor

      virtual ~G4HitRootIO();
      // Destructor

    public: // With description
      static G4HitRootIO* GetHitRootIO();
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

