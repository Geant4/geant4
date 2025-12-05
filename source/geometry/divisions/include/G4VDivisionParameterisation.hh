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
// G4VDivisionParameterisation
//
// Class description:
//
// Base class for parameterisations defining divisions of volumes
// for different kind of CSG and specific solids.

// Author: Pedro Arce (CIEMAT), 09.05.2001 - Initial version
//         Ivana Hrivnacova (Orsay), 08.04.2004 - Implemented reflection
//         Makoto Asai (SLAC), 21.04.2010 - Added gaps
//---------------------------------------------------------------------
#ifndef G4VDIVISIONPARAMETERISATION_HH
#define G4VDIVISIONPARAMETERISATION_HH 1

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4VPVParameterisation.hh"
#include "G4RotationMatrix.hh"

enum DivisionType { DivNDIVandWIDTH, DivNDIV, DivWIDTH };

class G4VPhysicalVolume;
class G4VSolid;

/**
 * @brief G4VDivisionParameterisation is the base class for parameterisations
 * defining divisions of volumes for different kind of CSG and specific solids.
 */

class G4VDivisionParameterisation : public G4VPVParameterisation
{ 
  public:
  
    /**
     * Constructor with number of divisions and width.
     *  @param[in] axis The axis along which do the division.
     *  @param[in] nDiv The number of division copies to replicate.
     *  @param[in] width The witdh of the divided slice along the axis.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @param[in] divType The kind of division type, based on input parameters.
     *  @param[in] motherSolid Pointer to the solid to be divided.
     */
    G4VDivisionParameterisation( EAxis axis,
                                 G4int nDiv,
                                 G4double width,
                                 G4double offset,
                                 DivisionType divType,
                                 G4VSolid* motherSolid = nullptr );

    /**
     * Destructor.
     */
    ~G4VDivisionParameterisation() override;
  
    /**
     * Invokes the base parameterisation to compute the parameterised solid.
     *  @param[in] pv Pointer to the physical volume to divide.
     *  @returns The pointer to the generated solid slice.
     */
    G4VSolid* ComputeSolid(const G4int, G4VPhysicalVolume* pv) override;

    /**
     * Base method to be implemented in derived classes for defining the
     * transformation in the parameterisation.
     */
    void ComputeTransformation(const G4int copyNo,
                               G4VPhysicalVolume* physVol) const override = 0;
  
    /**
     * Accessors and setters.
     */
    inline const G4String& GetType() const;
    inline EAxis GetAxis() const;
    inline G4int GetNoDiv() const;
    inline G4double GetWidth() const;
    inline G4double GetOffset() const;
    inline G4VSolid* GetMotherSolid() const;
    inline void SetType(const G4String& type);
    inline void SetHalfGap(G4double hg);
    inline G4double GetHalfGap() const;

  protected:

    /**
     * Defines and sets rotation on Z in trasformation for given volume.
     *  @param[in,out] physVol Pointer to the physical volume to apply rotation.
     *  @param[in] rotZ The rotation on Z (default is 0).
     */
    void ChangeRotMatrix( G4VPhysicalVolume* physVol,
                          G4double rotZ = 0.0 ) const;
  
    /**
     * Calculates the number of divisions, based on the input parameters.
     *  @param[in] motherDim The width of the volume to divide.
     *  @param[in] width The width of the division slice.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns The number of divisions.
     */
    G4int CalculateNDiv( G4double motherDim, G4double width,
                         G4double offset ) const;

    /**
     * Calculates the slice width, based on the input parameters.
     *  @param[in] motherDim The width of the volume to divide.
     *  @param[in] nDiv The number of division slices to create.
     *  @param[in] offset The optional offset distance from mother's border.
     *  @returns The width of the division slice.
     */
    G4double CalculateWidth( G4double motherDim, G4int nDiv,
                             G4double offset ) const;
  
    /**
     * Checks the validity of parameters given in input, issuing an exception.
     */
    virtual void CheckParametersValidity();

    /**
     * Internal methods for checking parameters.
     */
    void CheckOffset( G4double maxPar );
    void CheckNDivAndWidth( G4double maxPar );

    /**
     * Returns the max width along the axis of division.
     *  @returns The maximum width of the solid to divide along the axis.
     */
    virtual G4double GetMaxParameter() const = 0;

    /**
     * Sets the offset, taking into account potential reflection.
     */
    G4double OffsetZ() const;

  protected:
  
    G4String ftype;    // String for division type
    EAxis faxis;       // Axis of division
    G4int fnDiv = 0;   // Number of divisions
    G4double fwidth = 0.0;       // Width of division slice
    G4double foffset = 0.0;      // Offset distance from mother's border
    DivisionType fDivisionType;  // Division type ID
    G4VSolid* fmotherSolid = nullptr;  // Pointer to solid to divide
    G4bool fReflectedSolid = false;    // Specifies if solid is reflected
    G4bool fDeleteSolid = false;       // Specifies if delete solid to divide

    static G4ThreadLocal G4RotationMatrix* fRot;  // Rotation transformation

    static const G4int verbose;   // Verbosity level

    G4double kCarTolerance;       // Cached surface tolerance

    G4double fhgap = 0.0;         // Cached gap between slices
};

#include "G4VDivisionParameterisation.icc"

#endif
