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
// G4ParameterisationTrd[X/Y/Z]
//
// Class description:
//
// This class represents the parameterised positioning equivalent to 
// dividing a trapezoid along one of each axis X, Y, Z.

// Author: Pedro Arce (CIEMAT), 09.05.2001 - Initial version
//         Ivana Hrivnacova (Orsay), 08.04.2004 - Implemented reflection
// --------------------------------------------------------------------
#ifndef G4PARAMETERISATIONTRD_HH
#define G4PARAMETERISATIONTRD_HH 1

#include <vector>

#include "G4VDivisionParameterisation.hh"
#include "G4VSolid.hh" 

class G4VPhysicalVolume;

// Dummy declarations to get rid of warnings ...
//
class G4Cons;
class G4Box;
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
 * @brief G4VParameterisationTrd is the base class for the parameterised
 * positioning equivalent to dividing a trapezoid along one of each axis X, Y, Z.
 */

class G4VParameterisationTrd : public G4VDivisionParameterisation
{ 
  public:
  
    /**
     * Initialises a parameterised Trd, given the axis of parameterisation
     * 'axis' and the number of divided copies 'nCopies'.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4VParameterisationTrd( EAxis axis, G4int nCopies,
                            G4double offset, G4double step,
                            G4VSolid* pSolid, DivisionType divType );
  
    /**
     * Default Destructor.
     */
    ~G4VParameterisationTrd() override;

  protected:

    G4bool bDivInTrap = false;
};

/**
 * @brief G4ParameterisationTrdX represents the parameterised positioning
 * equivalent to dividing a G4Trd along X axis.
 */

class G4ParameterisationTrdX : public G4VParameterisationTrd
{ 
  public:

    /**
     * Initialises a parameterised Trd along X axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationTrdX( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationTrdX() override;

    /**
     * Returns the max width along X.
     *  @returns The maximum width of the solid to divide along the X axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume* physVol) const override;
    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;
    void ComputeDimensions(G4Trap& trd, const G4int copyNo,
			   const G4VPhysicalVolume* pv) const override;
  
  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
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
 * @brief G4ParameterisationTrdY represents the parameterised positioning
 * equivalent to dividing a G4Trd along Y axis.
 */

class G4ParameterisationTrdY : public G4VParameterisationTrd
{ 
  public:

    /**
     * Initialises a parameterised Trd along Y axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationTrdY( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationTrdY() override;

    /**
     * Returns the max width along Y.
     *  @returns The maximum width of the solid to divide along the Y axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume *physVol) const override;
    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;
    void ComputeDimensions (G4Trap&,const G4int,
                            const G4VPhysicalVolume*) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
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
 * @brief G4ParameterisationTrdZ represents the parameterised positioning
 * equivalent to dividing a G4Trd along Z axis.
 */

class G4ParameterisationTrdZ : public G4VParameterisationTrd
{ 
  public:

    /**
     * Initialises a parameterised Trd along Z axis.
     *  @param[in] axis The axis along which apply the parameterisation.
     *  @param[in] nCopies The total number of divided copies.
     *  @param[in] offset Potential initial offset along the axis.
     *  @param[in] step The width of the divided entity.
     *  @param[in] pSolid Pointer to the original shape to parameterise.
     *  @param[in] divType String identifier for the kind of division.
     */
    G4ParameterisationTrdZ( EAxis axis, G4int nCopies,
                            G4double width, G4double offset,
                            G4VSolid* pSolid, DivisionType divType );

    /**
     * Default Destructor.
     */
   ~G4ParameterisationTrdZ() override;

    /**
     * Returns the max width along Z.
     *  @returns The maximum width of the solid to divide along the Z axis.
     */
    G4double GetMaxParameter() const override;

    /**
     * Concrete methods implementing the parameterisation.
     */
    void ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume* physVol) const override;
    void ComputeDimensions(G4Trd& trd, const G4int copyNo,
                           const G4VPhysicalVolume* pv) const override;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Cons&,const G4int,
                            const G4VPhysicalVolume*) const override {}
    void ComputeDimensions (G4Box&,const G4int,
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
