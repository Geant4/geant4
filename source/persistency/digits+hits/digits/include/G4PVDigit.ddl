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
// $Id: G4PVDigit.ddl,v 1.3.4.1 2001/06/28 19:11:23 gunter Exp $
// GEANT4 tag $Name:  $
//

// Class Description:
//   This is a persistent version of abstract base class of a detector
// digit.  User should inherit from this class for his/her sensitive
// detector digit, to be added to the correponding digits collection.
//

#ifndef G4PVDigit_h
#define G4PVDigit_h 1

#include "G4PersistentSchema.hh"

#include "G4VDigi.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PVDigit
 : public HepPersObj, public G4VDigi
{

  public: // with description
      G4PVDigit();
      // Constructor.
      virtual ~G4PVDigit();
      // (Virtual) destructor.

  public:
      int operator==(const G4PVDigit &right) const;
      virtual void Draw();
      virtual void Print();
};

#endif

