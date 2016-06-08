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
// $Id: G4PSphere.ddl,v 1.7 2001/07/11 10:02:23 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $

// Class Description:
//   Persistent version of G4PSphere solid.

// History:
// 19.06.98 A.Kimura Converted G4Sphere.hh

#ifndef G4PSphere_DDL
#define G4PSphere_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Sphere;

#include "G4ThreeVector.hh"

class G4PSphere
 : public G4PCSGSolid
{
public: // With description
    G4PSphere(const G4Sphere* theSphere);
    virtual ~G4PSphere() ;
    // Constructor and Destructor

    G4VSolid* MakeTransientObject() const;
    // Creates a transient boolean solid object.

    virtual G4GeometryType  GetEntityType() const {return G4String("G4Sphere");}
    // Returns the G4GeometryType of this solid.

private:

    G4double fRmin,fRmax,
             fSPhi,fDPhi,
	     fSTheta,fDTheta;
};
   	
#endif

