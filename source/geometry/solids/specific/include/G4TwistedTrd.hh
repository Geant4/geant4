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
// G4TwistedTrd
//
// Class description:
//
// A G4TwistedTrd is a twisted trapezoid with the X and Y dimensions
// varying along Z.
//
//
// Member Data:
//
//   pDx1    Half-length along x at the surface positioned at -dz
//   pDx2    Half-length along x at the surface positioned at +dz
//   pDy1    Half-length along y at the surface positioned at -dz
//   pDy2    Half-length along y at the surface positioned at +dz
//   pDz     Half-length along z axis
//   pPhiTwist Twist angle

// Author: Oliver Link (CERN), 18.03.2005 - Created
// --------------------------------------------------------------------
#ifndef G4TWISTEDTRD_HH
#define G4TWISTEDTRD_HH

#include "G4VTwistedFaceted.hh"

/**
 * @brief G4TwistedTrd is a twisted trapezoid with the X and Y dimensions
 * varying along Z.
 */

class G4TwistedTrd : public G4VTwistedFaceted
{
  public:

    /**
     * Constructs a twisted Trd, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] pDx1 Half-length along X at the surface positioned at -dz.
     *  @param[in] pDx2 Half-length along X at the surface positioned at +dz.
     *  @param[in] pDy1 Half-length along Y at the surface positioned at -dz.
     *  @param[in] pDy2 Half-length along Y at the surface positioned at +dz.
     *  @param[in] pDz Half-length along Z axis.
     *  @param[in] pPhiTwist The twist angle.
     */
    G4TwistedTrd( const G4String& pName,
                        G4double  pDx1,
                        G4double  pDx2,
                        G4double  pDy1,
                        G4double  pDy2,
                        G4double  pDz,
                        G4double  pPhiTwist );

    /**
     * Default destructor.
     */
    ~G4TwistedTrd() override = default;

    /**
     * Accessors.
     */
    inline G4double GetX1HalfLength() const { return GetDx1() ; }
    inline G4double GetX2HalfLength() const { return GetDx3() ; }
    inline G4double GetY1HalfLength() const { return GetDy1() ; }
    inline G4double GetY2HalfLength() const { return GetDy2() ; }
    inline G4double GetZHalfLength()  const { return GetDz()  ; }
    inline G4double GetPhiTwist()     const { return GetTwistAngle() ; }

    /**
     * Returns the type ID, "G4TwistedTrd" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Returning an estimation of the solid volume (capacity) and
     * surface area, in internal units.
     */
    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream&  StreamInfo(std::ostream& os) const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TwistedTrd(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4TwistedTrd(const G4TwistedTrd& rhs);
    G4TwistedTrd& operator=(const G4TwistedTrd& rhs);
} ;

#endif
