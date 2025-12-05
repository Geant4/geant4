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
// G4ParameterisationBox[X/Y/Z]
//
// Class description:
//
// These classes represent the parameterised positioning equivalent to 
// dividing a G4Box along one of each axis X, Y, Z.

// Author: Pedro Arce (CIEMAT), 09.05.2001 - Initial version
//         Ivana Hrivnacova (Orsay), 08.04.2004 - Implemented reflection
// --------------------------------------------------------------------
#ifndef G4PARAMETERISATIONBOX_HH
#define G4PARAMETERISATIONBOX_HH 1

#include "G4VDivisionParameterisation.hh"

class G4VSolid;
class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Cons;
class G4Trd;
class G4Sphere;
class G4Orb;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

/**
 * @brief G4VParameterisationBox is the base class for the parameterised
 * positioning equivalent to dividing a G4Box along one of each axis X, Y, Z.
 */

class G4VParameterisationBox : public G4VDivisionParameterisation
{ 
  public:
  
    /**
     * Initialises a parameterised box, given the axis of parameterisation
     * 'axis' and the number of divided copies 'nCopies'.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4VParameterisationBox( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* pSolid, DivisionType divType );
  
    /**
     * Default Destructor.
     */
    ~G4VParameterisationBox() override;
};

/**
 * @brief G4ParameterisationBoxX represents the parameterised positioning
 * equivalent to dividing a G4Box along X axis.
 */

class G4ParameterisationBoxX : public G4VParameterisationBox
{ 
  public:
  
    /**
     * Initialises a parameterised box along X axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationBoxX( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationBoxX() override;

    /**
     * Returns the max width along X.
     *  @returns The maximum width of the solid to divide along the X axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions(G4Box& box, const G4int copyNo,
                           const G4VPhysicalVolume* physVol) const override;

   private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

/**
 * @brief G4ParameterisationBoxY represents the parameterised positioning
 * equivalent to dividing a G4Box along Y axis.
 */

class G4ParameterisationBoxY : public G4VParameterisationBox
{ 
  public:

    /**
     * Initialises a parameterised box along Y axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationBoxY( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationBoxY() override;


    /**
     * Returns the max width along Y.
     *  @returns The maximum width of the solid to divide along the Y axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation( const G4int copyNo,
                                      G4VPhysicalVolume* physVol ) const override;
    void ComputeDimensions(G4Box& box, const G4int copyNo,
                           const G4VPhysicalVolume* physVol) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

/**
 * @brief G4ParameterisationBoxZ represents the parameterised positioning
 * equivalent to dividing a G4Box along Z axis.
 */

class G4ParameterisationBoxZ : public G4VParameterisationBox
{ 
  public:

    /**
     * Initialises a parameterised box along Z axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationBoxZ( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationBoxZ() override;


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
    void ComputeDimensions(G4Box& box, const G4int copyNo,
                           const G4VPhysicalVolume* physVol) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trd&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Sphere&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Orb&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Torus&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Para&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Hype&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Tubs&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polycone&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Polyhedra&,const G4int,
                            const G4VPhysicalVolume*) const override {}
};

#endif
