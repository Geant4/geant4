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
// $Id: G4PIntersectionSolid.ddl,v 1.5 2001/07/11 10:02:21 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $

// Class Description:
//   Persistent class for description of intersection of two CSG solids

// History: 
// 10.11.99 Y.Morita, Initial creation

#ifndef G4PINTERSECTIONSOLID_DDL
#define G4PINTERSECTIONSOLID_DDL

#include "G4PersistentSchema.hh"
#include "G4PBooleanSolid.hh"

class G4VSolid;

class G4PIntersectionSolid: public G4PBooleanSolid
{
public: // With description
        G4PIntersectionSolid( const G4String& pName,
                              HepRef(G4PVSolid) persSolidA,
                              HepRef(G4PVSolid) persSolidB );
        virtual ~G4PIntersectionSolid() ;
        // Constructor and destructor

        G4VSolid* MakeTransientObject() const;
        G4VSolid* MakeTransientBooleanSolid(
                               G4VSolid* aSolidA,
                               G4VSolid* aSolidB ) const;
        // Creates a transient G4IntersectionSolid object.

        virtual G4GeometryType GetEntityType() const 
        { return G4String("G4IntersectionSolid"); }
        // Returns the G4GeometryType of this solid.

protected:

private:

};

#endif

