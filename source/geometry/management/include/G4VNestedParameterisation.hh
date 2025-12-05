//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// G4VNestedParameterisation
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

// Author: John Apostolakis (CERN), 24.02.2005 - First created version
// --------------------------------------------------------------------
#ifndef G4VNESTEDPARAMETERISATION_HH
#define G4VNESTEDPARAMETERISATION_HH

#include "G4Types.hh"
#include "G4VPVParameterisation.hh" 
#include "G4VVolumeMaterialScanner.hh"
#include "G4VTouchable.hh"

class G4VPhysicalVolume;
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
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Polycone;
class G4Polyhedra;
class G4Hype;

/**
 * @brief G4VNestedParameterisation is a base class for parameterisations that
 * use information from the parent volume to compute the material of a
 * copy/instance of this volume. This is in addition to using the current
 * replication number. Such a volume can be nested inside a placement volume
 * or a parameterised volume. The user can modify the solid type, size or
 * transformation using only the replication number of this parameterised
 * volume; it is NOT allowed to change these attributes using information of
 * the parent volumes, otherwise incorrect results will occur.
 */

class G4VNestedParameterisation : public G4VPVParameterisation, 
                                  public G4VVolumeMaterialScanner
{
  public:

    /**
     * Default Constructor & Destructor.
     */
    G4VNestedParameterisation() = default; 
    ~G4VNestedParameterisation() override = default; 

    // Methods required in derived classes
    // -----------------------------------

    /**
     * Computes the material for the 'currentVol' and replica number 'repNo'.
     * It is a required method, as it is the reason for this class.
     * Must cope with parentTouch=nullptr for navigator's SetupHierarchy().
     *  @param[in] currentVol Pointer to the current physical volume.
     *  @param[in] repNo The copy number index.
     *  @param[in] parentTouch Pointer to the touchable of the parent volume.
     *  @returns A pointer to the associated material.
     */
    virtual G4Material* ComputeMaterial(G4VPhysicalVolume* currentVol,
                                 const G4int repNo, 
                                 const G4VTouchable* parentTouch = nullptr)=0;

    /**
     * Accessors needed to define materials for instances of a Nested
     * Parameterisation. Eeach call should return the materials of all
     * instances with the same mother/ancestor volume.
     */
    G4int GetNumberOfMaterials() const override =0;
    G4Material* GetMaterial(G4int idx) const override =0;

    /**
     * Computes the transformation for the 'currentPV' and replica number 'no'.
     * It is a required method, as it is the reason for this class.
     * Must cope with parentTouch=nullptr for navigator's SetupHierarchy().
     *  @param[in] currentPV Pointer to the current physical volume.
     *  @param[in] no The copy number index.
     */
     void ComputeTransformation(const G4int no,
                                G4VPhysicalVolume* currentPV) const override=0;

    // Methods optional in derived classes
    // -----------------------------------

    /**
     * Computes the solid for the 'thisVol' volume and replica number 'no'.
     * To be optionally defined in derived classes, for parameterisation of
     * the solid type.
     *  @param[in] no The copy number index.
     *  @param[in] thisVol Pointer to the current physical volume.
     */
    G4VSolid* ComputeSolid(const G4int no, G4VPhysicalVolume* thisVol) override;

    /**
     * Method implemented in this class in terms of the above ComputeMaterial().
     */
    G4Material* ComputeMaterial(const G4int repNo, 
                           G4VPhysicalVolume* currentVol,
                           const G4VTouchable* parentTouch = nullptr) override;

    /**
     * Methods to identify nested parameterisations. Required in order
     * to enable material scan for nested parameterisations.
     */
    G4bool IsNested() const override;
    G4VVolumeMaterialScanner* GetMaterialScanner() override; 

  // Dispatch methods for the specific solids where parameterisation is allowed
  // --------------------------------------------------------------------------

    void ComputeDimensions(G4Box &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Tubs &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Trd &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Trap &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Cons &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Sphere &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Orb &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Ellipsoid &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Torus &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Para &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Polycone &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Polyhedra &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}

    void ComputeDimensions(G4Hype &,
                           const G4int,
                           const G4VPhysicalVolume *) const override {}
};

#endif
