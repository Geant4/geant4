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
// $Id: G4PCons.ddl,v 1.3.4.1 2001/06/28 19:11:31 gunter Exp $
// GEANT4 tag $Name:  $
//
// class G4Cons
//
// History:
// 19.06.98 A.Kimura Converted G4Cons.hh

#ifndef G4PCons_DDL
#define G4PCons_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4Cons;
class G4VSolid;

class G4PCons : public G4PCSGSolid {
public:
    G4PCons(const G4Cons* theCons);
    virtual ~G4PCons() ;

    G4VSolid* MakeTransientObject() const;

    // Naming method (pseudo-RTTI : run-time type identification
    virtual G4GeometryType  GetEntityType() const {return G4String("G4Cons");}

private:

    G4double fRmin1, fRmin2,
             fRmax1, fRmax2,
             fDz,
	     fSPhi, fDPhi;
};

#endif
