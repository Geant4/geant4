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
// $Id: G4PTubs.ddl,v 1.5 2001/07/11 10:02:24 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $

// Class Description:
//   Persistent version of G4PTubs solid.

// History:
// 19.06.98 A.Kimura Converted G4Tubs.hh

#ifndef G4PTUBS_DDL
#define G4PTUBS_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Tubs;

class G4PTubs
 : public G4PCSGSolid
{
public: // With description
    G4PTubs(const G4Tubs* theTubs);
    virtual ~G4PTubs();
    // Constructor and Destructor

    G4VSolid* MakeTransientObject() const;
    // Creates a transient boolean solid object.

    virtual G4GeometryType  GetEntityType() const { return G4String("G4Tubs"); }
    // Returns the G4GeometryType of this solid.

protected:

    G4double fRMin,fRMax,fDz,fSPhi,fDPhi;

};
   	
#endif

