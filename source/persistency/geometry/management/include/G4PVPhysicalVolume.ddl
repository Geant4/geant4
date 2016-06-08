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
// $Id: G4PVPhysicalVolume.ddl,v 1.7 2001/07/11 10:02:19 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $

// Class Description:
//   Persistent version of G4VPhysicalVolume.

#ifndef G4PVPHYSICALVOLUME_DDL
#define G4PVPHYSICALVOLUME_DDL 1

#include "G4Pglobals.hh"
#include "geomdefs.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4VPhysicalVolume;
class G4LogicalVolume;

//class G4VPVParameterisation;

// forward declaration and #pragma ooclassref is necessary in solving
// the circular dependency between G4PLogicalVolume and G4PVPhysicalVolume
class G4PLogicalVolume;
#pragma ooclassref G4PLogicalVolume "G4PLogicalVolume_ref.hh"

class G4PVPhysicalVolume: public HepPersObj
{
public: // With description
    G4PVPhysicalVolume();
    G4PVPhysicalVolume( const G4VPhysicalVolume *PhysVol,
                        const HepRef(G4PLogicalVolume) persLogVol);
    ~G4PVPhysicalVolume();
    // Constructor and Destructor

    void SetMother(const HepRef(G4PVPhysicalVolume) persMother);
    // Set the persistent mother volume.

    G4VPhysicalVolume* MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother );
    // Creates a transient G4VPhysicalVolume object.

    G4PString GetName();
    void SetName(const G4PString aName);
    // Get and set the name of the volume.

    G4RotationMatrix* GetRotation();
    void SetRotation(const G4RotationMatrix* aRot);
    // Get and set the rotation matrix. 

    G4ThreeVector&  GetTranslation();
    void SetTranslation(G4ThreeVector atrans);
    // Get and set the translation vector.

    HepRef(G4PLogicalVolume) GetLogicalVolume();
    // Returns the pointer to the persistent logical volume.

    G4bool operator == (const G4PVPhysicalVolume& p) const;
    // Define equality by equal addresses only.

    virtual G4bool IsMany() const = 0;
    virtual G4int GetCopyNo() const = 0;
    virtual void  SetCopyNo(G4int CopyNo) = 0;
    // Functions required of subclasses

protected:
    d_Varray<d_Double> frot;
    d_Varray<d_Double> ftrans;

    G4ThreeVector ftransvector;

    d_Ref<G4PLogicalVolume> flogical;
    G4PString fname;	            	// name of the volume
    d_Ref<G4PVPhysicalVolume> fmother;	// The current moher volume

};

#endif

