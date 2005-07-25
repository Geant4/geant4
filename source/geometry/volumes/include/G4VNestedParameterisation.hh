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
// $Id: G4VNestedParameterisation.hh,v 1.3 2005-07-25 10:02:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VNestedParameterisation
//
// Class description:
//
// Base class for parameterisations that use information from the parent
// volume to compute the material of a copy/instance of this volume. 
// This is in addition to using the current replication number.
// 
// Notes:
//  - Such a volume can be nested inside a placement volume or a parameterised  
//    volume.
//  - The user can modify the solid type, size or transformation using only
//    the replication number of this parameterised volume.
//    He/she is NOT allowed to change these attributes using information of
//    parent volumes - otherwise incorrect results will occur.
//  Also note that the usual restrictions apply: 
//   - the mother volume, in which these copies are placed, must always be
//     of the same dimensions

// History:
// 24.02.05 - J.Apostolakis - First created version.
// --------------------------------------------------------------------
#ifndef G4VNESTEDPARAMETERISATION_HH
#define G4VNESTEDPARAMETERISATION_HH

#include "G4Types.hh"
#include "G4VPVParameterisation.hh" 

class G4VPhysicalVolume;
class G4VTouchable; 
class G4VSolid;
class G4Material;

// CSG Entities which may be parameterised/replicated
//
class G4Box;
class G4Tubs;
class G4Trd;
class G4Trap;
class G4Cons;
class G4Sphere;
class G4Orb;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

class G4VNestedParameterisation: public G4VPVParameterisation
{
  public:  // with description

    virtual G4Material* ComputeMaterial(const G4int repNo, 
                                        G4VPhysicalVolume *currentVol,
                                        const G4VTouchable *parentTouch) = 0;
      // Mandatory method, required as reason for this class

    virtual void ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume *currentPV) const = 0;

    virtual G4VSolid* ComputeSolid(const G4int no, G4VPhysicalVolume *thisVol);

    virtual G4Material* ComputeMaterial(const G4int, G4VPhysicalVolume *);
      // This methood is implemented in terms of the methods above 
      // providing 're-interpretation' to add information of parent

    // Additional standard Parameterisation methods

    virtual void ComputeDimensions(G4Box &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Tubs &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Trd &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Trap &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Cons &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Sphere &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Orb &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Torus &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Para &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Polycone &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Polyhedra &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

    virtual void ComputeDimensions(G4Hype &,
                                   const G4int,
                                   const G4VPhysicalVolume *) const {}

  private:

    void ReportErrorInTouchable(const G4String& method, 
                                const G4VPhysicalVolume* thisVol) const; 

    G4VPhysicalVolume* ObtainParts( G4VPhysicalVolume* thisVolPlus,
                       const G4String& method, 
                       const G4VTouchable** pPtrTouchableParent) const;
    const G4VPhysicalVolume* ObtainParts( const G4VPhysicalVolume* thisVolPlus,
                       const G4String& method, 
                       const G4VTouchable** pPtrTouchableParent) const;
       // Demux the expected PhysicalTouchable into 
       //   - pointer to this physical volume
       //   - pointer to parent touchable
};

#endif
