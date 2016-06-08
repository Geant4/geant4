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
// $Id: G4PBooleanSolid.ddl,v 1.6 2001/07/11 10:02:21 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $

// Class Description:
//   Persistent base class for solids created by boolean operations
//   between other solids

// History:
// 10.11.99 Y.Morita, Initial Creation

#ifndef G4PBOOLEANSOLID_DDL
#define G4PBOOLEANSOLID_DDL

#include "G4PVSolid.hh"
#include "G4PersistentSchema.hh"

class G4PBooleanSolid
 : public G4PVSolid
{
public: // With description
          G4PBooleanSolid( const G4String& pName,
                           HepRef(G4PVSolid) persSolidA,
                           HepRef(G4PVSolid) persSolidB );
          virtual ~G4PBooleanSolid();
            // Constructor and Destructor

          virtual G4VSolid* MakeTransientBooleanSolid(
                               G4VSolid* aSolidA,
                               G4VSolid* aSolidB ) const = 0;
            // Creates a transient boolean solid object.

public:  // With Description
    // If Solid is made up from a Boolean operation of two solids,
    //   return the corresponding solid (for no=0 and 1)
    // If the solid is not a "Boolean", return 0
    virtual const HepRef(G4PVSolid) GetConstituentSolid(G4int no) const;
    virtual       HepRef(G4PVSolid) GetConstituentSolid(G4int no);

protected:
          d_Ref<G4PVSolid> fPtrSolidA;
          d_Ref<G4PVSolid> fPtrSolidB;

private:
          G4bool  createdDisplacedSolid;
};

#endif

