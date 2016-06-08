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
// $Id: G4PTrd.ddl,v 1.5 2001/07/11 10:02:23 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $

// Class Description:
//   Persistent version of G4PTrd solid.

// History:
// 19.06.98 A.Kimura Converted G4Trd.hh

#ifndef G4PTRD_DDL
#define G4PTRD_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4Trd;
class G4VSolid;

class G4PTrd
 : public G4PCSGSolid
{
public: // With description
    G4PTrd(const G4Trd* theTrd);
    virtual ~G4PTrd();
    // Constructor and Destructor

    G4VSolid* MakeTransientObject() const;
    // Creates a transient boolean solid object.

    virtual G4GeometryType  GetEntityType() const {return G4String("G4Trd");}
    // Returns the G4GeometryType of this solid.

public:
    void CheckAndSetAllParameters (G4double pdx1, G4double pdx2,
                             G4double pdy1, G4double pdy2,
                             G4double pdz);

protected:

    G4double fDx1,fDx2,fDy1,fDy2,fDz;

};

#endif


