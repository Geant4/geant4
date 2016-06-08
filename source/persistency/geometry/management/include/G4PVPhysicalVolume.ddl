// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PVPhysicalVolume.ddl,v 1.4 1999/12/15 14:51:23 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
// 
// class G4PVPhysicalVolume

#ifndef G4PVPHYSICALVOLUME_DDL
#define G4PVPHYSICALVOLUME_DDL 1

#include "globals.hh"
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
public:
    G4PVPhysicalVolume();

    G4PVPhysicalVolume( const G4VPhysicalVolume *PhysVol,
                        const HepRef(G4PLogicalVolume) persLogVol);

// Destructor
    ~G4PVPhysicalVolume();
//
    void SetMother(const HepRef(G4PVPhysicalVolume) persMother);

    G4VPhysicalVolume* MakeTransientObject(
                             G4LogicalVolume* aLogical,
                             G4VPhysicalVolume* aMother );

// Access functions
    G4PString GetName();
    void SetName(const G4PString aName);
    G4RotationMatrix* GetRotation();
    G4ThreeVector&  GetTranslation();
    void SetRotation(const G4RotationMatrix* aRot);
    void SetTranslation(G4ThreeVector atrans);

    HepRef(G4PLogicalVolume) GetLogicalVolume();

// Define equality by equal addresses only.
    G4bool operator == (const G4PVPhysicalVolume& p) const;

// Functions required of subclasses
    virtual G4bool IsMany() const = 0;
    virtual G4int GetCopyNo() const = 0;
    virtual void  SetCopyNo(G4int CopyNo) = 0;

protected:
    d_Varray<d_Double> frot;
    d_Varray<d_Double> ftrans;

    G4ThreeVector ftransvector;

    d_Ref<G4PLogicalVolume> flogical;
    G4PString fname;	            	// name of the volume
    d_Ref<G4PVPhysicalVolume> fmother;	// The current moher volume

};

#endif

