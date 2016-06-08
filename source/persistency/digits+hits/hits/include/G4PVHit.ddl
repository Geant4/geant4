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
// $Id: G4PVHit.ddl,v 1.7.4.1 2001/06/28 19:11:24 gunter Exp $
// GEANT4 tag $Name:  $
//

// Class Description:
//   This is a persistent version of abstract base class of a detector
// digit.  User should inherit from this class for his/her sensitive
// detector digit, to be added to the correponding digits collection.
//

#ifndef G4PVHit_h
#define G4PVHit_h 1

#include "G4PersistentSchema.hh"

#include "G4VHit.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVHit
 : public HepPersObj, public G4VHit
{

  public: // with description
      G4PVHit();
      // Constructor.
      virtual ~G4PVHit();
      // (Virtual) destructor.

  public:
      int operator==(const G4PVHit &right) const;
      virtual void Draw();
      virtual void Print();
};

#endif

