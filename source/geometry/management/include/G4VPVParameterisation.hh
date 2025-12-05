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
// G4VPVParameterisation
//
// Class description:
//
// Parameterisation abstract base class, able to compute the transformation
// and (indirectly) the dimensions of parameterised volumes, given a
// replication number.

// Author: Paul Kent (CERN), 25.07.1996, P.Kent - Initial stub version
// --------------------------------------------------------------------
#ifndef G4VPVPARAMETERISATION_HH
#define G4VPVPARAMETERISATION_HH

#include "G4Types.hh"
#include "G4VVolumeMaterialScanner.hh"
#include "G4VTouchable.hh"

class G4VPhysicalVolume;
class G4VSolid;
class G4Material;

// Entities which may be parameterised/replicated
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

class G4VVolumeMaterialScanner; 

/**
 * @brief G4VPVParameterisation ia an abstract base class for Parameterisation,
 * able to compute the transformation and (indirectly) the dimensions of
 * parameterised volumes, given a replication number.
 */

class G4VPVParameterisation
{
  public:

    /**
     * Default Constructor & Destructor.
     */
    G4VPVParameterisation() = default;
    virtual ~G4VPVParameterisation() = default;

    /**
     * Computes the transformation for the 'pv' volume and replica number 'no'.
     * It is a required method, as it is the reason for this class.
     *  @param[in] pv Pointer to the current physical volume.
     *  @param[in] no The copy number index.
     */
    virtual void ComputeTransformation(const G4int no,
                                       G4VPhysicalVolume* pv) const = 0;

    /**
     * Computes the solid for the 'pv' volume and replica number 'no'.
     * To be optionally defined in derived classes, for parameterisation of
     * the solid type.
     *  @param[in] no The copy number index.
     *  @param[in] pv Pointer to the current physical volume.
     */
    virtual G4VSolid* ComputeSolid(const G4int no, G4VPhysicalVolume* pv);
				       
    /**
     * Computes the material for the 'currentVol' and replica number 'repNo'.
     * Must cope with 'parentTouch' for navigator's SetupHierarchy() when
     * used for nested parameterisations.
     *  @param[in] currentVol Pointer to the current physical volume.
     *  @param[in] repNo The copy number index.
     *  @param[in] parentTouch Pointer to the touchable of the parent volume.
     *  @returns A pointer to the associated material.
     */
    virtual G4Material* ComputeMaterial(const G4int repNo, 
                                    G4VPhysicalVolume* currentVol,
                                    const G4VTouchable* parentTouch = nullptr);

    /**
     * Methods to identify nested parameterisations. Required in order
     * to enable material scan for nested parameterisations.
     */
    virtual G4bool IsNested() const;
    virtual G4VVolumeMaterialScanner* GetMaterialScanner(); 

    /**
     * Dispatch methods for the specific solids where parameterisation
     * is allowed.
     */
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
    virtual void ComputeDimensions(G4Ellipsoid &,
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
};

#endif
