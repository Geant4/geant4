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
// G4ParameterisationPolycone[Rho/Phi/Z]
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Polycone along one of each axis Rho, Phi, Z.

// Author: Pedro Arce (CIEMAT), 09.05.2001 - Initial version
//         Ivana Hrivnacova (Orsay), 08.04.2004 - Implemented reflection
//---------------------------------------------------------------------
#ifndef G4PARAMETERISATIONPOLYCONE_HH
#define G4PARAMETERISATIONPOLYCONE_HH 1

#include "G4VDivisionParameterisation.hh"
#include "G4Polycone.hh"

class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Ellipsoid;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polyhedra;

/**
 * @brief G4VParameterisationPolycone is the base class for the parameterised
 * positioning equivalent to dividing a G4Polycone along one of each axis Rho,
 * Phi, Z.
 */

class G4VParameterisationPolycone : public G4VDivisionParameterisation
{ 
  public:
  
    /**
     * Initialises a parameterised polycone, given the axis of parameterisation
     * 'axis' and the number of divided slices 'nCopies'.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided slices.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided slice.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4VParameterisationPolycone( EAxis axis, G4int nCopies,
                                 G4double offset, G4double step,
                                 G4VSolid* pSolid, DivisionType divType );
  
    /**
     * Default Destructor.
     */
    ~G4VParameterisationPolycone() override;
};

/**
 * @brief G4ParameterisationPolyconeRho represents the parameterised positioning
 * equivalent to dividing a G4Polycone along Rho axis.
 */

class G4ParameterisationPolyconeRho : public G4VParameterisationPolycone
{ 
  public:

    /**
     * Initialises a parameterised polycone, along the Rho axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided slices.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided slice.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationPolyconeRho( EAxis axis, G4int nCopies,
                                   G4double offset, G4double step,
                                   G4VSolid* pSolid,
                                   DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationPolyconeRho() override;

    /**
     * Checks the validity of parameters given in input, issuing an exception.
     */
    void CheckParametersValidity() override;

    /**
     * Returns the max width along Rho.
     *  @returns The maximum width of the solid to divide along the Rho axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

/**
 * @brief G4ParameterisationPolyconePhi represents the parameterised positioning
 * equivalent to dividing a G4Polycone along Phi axis.
 */

class G4ParameterisationPolyconePhi : public G4VParameterisationPolycone
{ 
  public:

    /**
     * Initialises a parameterised polycone, along the Phi axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided slices.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided slice.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationPolyconePhi( EAxis axis, G4int nCopies,
                                   G4double offset, G4double step,
                                   G4VSolid* pSolid,
                                   DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationPolyconePhi() override;

    /**
     * Returns the max width along Phi.
     *  @returns The maximum width of the solid to divide along the Phi axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

/**
 * @brief G4ParameterisationPolyconeZ represents the parameterised positioning
 * equivalent to dividing a G4Polycone along Z axis.
 */

class G4ParameterisationPolyconeZ : public G4VParameterisationPolycone
{ 
  public:

    /**
     * Initialises a parameterised polycone, along the Z axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided slices.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided slice.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationPolyconeZ( EAxis axis, G4int nCopies,
                                 G4double offset, G4double step,
                                 G4VSolid* pSolid,
                                 DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationPolyconeZ() override;

    /**
     * Checks the validity of parameters given in input, issuing an exception.
     */
    void CheckParametersValidity() override;

    /**
     * Returns the max width along Z.
     *  @returns The maximum width of the solid to divide along the Z axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions( G4Polycone& pcone, const G4int copyNo,
                            const G4VPhysicalVolume* physVol ) const override;

  private:

    /**
     * Internal accessors for the original R parameters of the solid to divide.
     */
    G4double GetR(G4double z, G4double z1, G4double r1,
                  G4double z2, G4double r2) const;
    G4double GetRmin(G4double z, G4int nsegment) const;
    G4double GetRmax(G4double z, G4int nsegment) const;

    // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
  private:

    G4int fNSegment = 0;
    G4PolyconeHistorical* fOrigParamMother = nullptr;
};

#endif
