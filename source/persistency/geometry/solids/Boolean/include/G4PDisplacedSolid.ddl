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
// $Id: G4PDisplacedSolid.ddl,v 1.5 2001/07/11 10:02:21 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $

// Class Description:
// Persistent class describing solid placements for boolean operations

// History:
// 10.11.99 Y.Morita, Initial creation

#ifndef G4PDisplacedSolid_DDL
#define G4PDisplacedSolid_DDL

#include "G4PersistentSchema.hh"
#include "G4PVSolid.hh"
#include "G4PAffineTransform.hh"

class G4VSolid;

class G4PDisplacedSolid
 : public G4PVSolid
{
public: // With description
        G4PDisplacedSolid ( HepRef(G4PVSolid) persCostituentSolid,
                            HepRef(G4PAffineTransform) pDirectTransform );
        virtual ~G4PDisplacedSolid();
        // Constructor and destructor.

        G4VSolid* MakeTransientObject() const;
        G4VSolid* MakeTransientDisplacedSolid(G4VSolid* aSolid) const;
        // Creates a transient G4VSolid object.

        virtual G4GeometryType  GetEntityType() const
        { return G4String("G4DisplacedSolid"); }
        // Returns the G4GeometryType of the solid.

        HepRef(G4PVSolid) GetConstituentMovedSolid();
        // Returns the pointer of persistent GetConstituentMovedSolid object.

protected:
        d_Ref<G4PVSolid> fPtrSolid;
//      fPtrTransform can be created by fDirectTransform.Inverse()
        d_Ref<G4PAffineTransform> fDirectTransform;

private:

};

#endif
