// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PHype.ddl,v 1.2.2.1 1999/12/07 20:50:11 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// class G4PHype
//
// History:
// 19.06.98 A.Kimura Converted G4Hype.hh

#ifndef G4PHYPE_DDL
#define G4PHYPE_DDL

#include "G4PersistentSchema.hh"
#include "G4PCSGSolid.hh"

class G4VSolid;
class G4Hype;

class G4PHype : public G4PCSGSolid {
public:

    G4PHype(const G4Hype* theHype);
    virtual ~G4PHype();

    G4VSolid* MakeTransientObject() const;

    virtual G4GeometryType  GetEntityType() const {return G4String("G4Hype");}

protected:
    // values of hype radius at a given Z
    G4double HypeInnerRadius2(double zVal) const {
	return (tanInnerStereo2*zVal*zVal+innerRadius2);
    }
    G4double HypeOuterRadius2(double zVal) const {
	return (tanOuterStereo2*zVal*zVal+outerRadius2);
    }

    G4double innerRadius; // variable names are quite self explanative
    G4double outerRadius;
    G4double halfLenZ;
    G4double innerStereo;
    G4double outerStereo;

    // precalculated parameters, squared quantities
    G4double tanInnerStereo2; // squared tan of Inner Stereo angle
    G4double tanOuterStereo2; // squared tan of Outer Stereo angle
    G4double innerRadius2;    // squared Inner Radius
    G4double outerRadius2;    // squared Outer Radius
    G4double endInnerRadius2; // squared endcap Inner Radius
    G4double endOuterRadius2; // squared endcap Outer Radius
    G4double endInnerRadius;  // endcap Inner Radius
    G4double endOuterRadius;  // endcap Outer Radius
};
   	
#endif



































