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
// G4TwistedTrap
//
// Class description:
//
// A G4TwistedTrap is a general twisted trapezoid: The faces perpendicular to the
// z planes are trapezia, and their centres are not necessarily on
// a line parallel to the z axis.
//
//     pDz     Half-length along the z-axis
//     pTheta  Polar angle of the line joining the centres of the faces
//             at -/+pDz
//     pPhi    Azimuthal angle of the line joining the centre of the face at
//             -pDz to the centre of the face at +pDz
//     pDy1    Half-length along y of the face at -pDz
//     pDx1    Half-length along x of the side at y=-pDy1 of the face at -pDz
//     pDx2    Half-length along x of the side at y=+pDy1 of the face at -pDz
//
//     pDy2    Half-length along y of the face at +pDz
//     pDx3    Half-length along x of the side at y=-pDy2 of the face at +pDz
//     pDx4    Half-length along x of the side at y=+pDy2 of the face at +pDz
//     pAlph   Angle with respect to the y axis from the centre of the side
//
//
// A special regular case of a trapezoid with equal endcaps is available, 
// with polar,azimuthal and tilt angles set to zero.
//
 
// Author: Oliver Link (CERN), 27.10.2004 - Created
// --------------------------------------------------------------------
#ifndef G4TWISTEDTRAP_HH
#define G4TWISTEDTRAP_HH

#include "G4VTwistedFaceted.hh"

/**
 * @brief G4TwistedTrap is a general twisted trapezoid. The faces perpendicular
 * to the Z planes are trapezia, and their centres are not necessarily on
 * a line parallel to the Z axis.
 */

class G4TwistedTrap : public G4VTwistedFaceted
{
  public:

    /**
     * Constructs a twisted trapezoid, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] pPhiTwist The twist angle.
     *  @param[in] pDx1 Half X length at -pDz,-pDy.
     *  @param[in] pDx2 Half X length at -pDz,+pDy.
     *  @param[in] pDy Half-length along the Y axis.
     *  @param[in] pDz Half-length along the Z axis.
     */
    G4TwistedTrap(const G4String& pName,
                        G4double  pPhiTwist,
                        G4double  pDx1,  // half x length at -pDz,-pDy
                        G4double  pDx2,  // half x length at -pDz,+pDy
                        G4double  pDy,
                        G4double  pDz);

    /**
     * Alternative construct of a twisted trapezoid, given its parameters.
     *  @param[in] pName The solid name.
     *  @param[in] pPhiTwist The twist angle.
     *  @param[in] pDz Half-length along the Z axis.
     *  @param[in] pTheta Polar angle of the line joining the centres of the
     *             faces at -/+pDz.
     *  @param[in] pPhi Azimuthal angle of the line joining the centres of the
     *             faces at -/+pDz.
     *  @param[in] pDy1 Half Y length at -pDz.
     *  @param[in] pDx1 Half X length at -pDz, y=-pDy1.
     *  @param[in] pDx2 Half X length at -pDz, y=+pDy1.
     *  @param[in] pDy2 Half Y length at +pDz.
     *  @param[in] pDx3 Half X length at +pDz, y=-pDy2.
     *  @param[in] pDx4 Half X length at +pDz, y=+pDy2.
     *  @param[in] pAlph Angle with respect to the Y axis from centre of side.
     */
    G4TwistedTrap(const G4String& pName,       // Name of instance
                        G4double  pPhiTwist,   // twist angle
                        G4double  pDz,     // half z length
                        G4double  pTheta,  // direction between end planes
                        G4double  pPhi,    // defined by polar and azim. angles
                        G4double  pDy1,    // half y length at -pDz
                        G4double  pDx1,    // half x length at -pDz,-pDy
                        G4double  pDx2,    // half x length at -pDz,+pDy
                        G4double  pDy2,    // half y length at +pDz
                        G4double  pDx3,    // half x length at +pDz,-pDy
                        G4double  pDx4,    // half x length at +pDz,+pDy
                        G4double  pAlph    // tilt angle
                  );
 
    /**
     * Default destructor.
     */
    ~G4TwistedTrap() override = default;

    /**
     * Accessors.
     */
    inline G4double GetY1HalfLength() const { return GetDy1() ; }
    inline G4double GetX1HalfLength() const { return GetDx1() ; }
    inline G4double GetX2HalfLength() const { return GetDx2() ; }
    inline G4double GetY2HalfLength() const { return GetDy2() ; }
    inline G4double GetX3HalfLength() const { return GetDx3() ; }
    inline G4double GetX4HalfLength() const { return GetDx4() ; }
    inline G4double GetZHalfLength()  const { return GetDz()  ; }
    inline G4double GetPhiTwist()     const { return GetTwistAngle() ; }
    inline G4double GetTiltAngleAlpha()    const { return GetAlpha() ; }
    inline G4double GetPolarAngleTheta()   const { return GetTheta() ; }
    inline G4double GetAzimuthalAnglePhi() const { return GetPhi()   ; }

    /**
     * Returns the type ID, "G4TwistedTrap" of the solid.
     */
    G4GeometryType GetEntityType() const override;

    /**
     * Makes a clone of the object for use in multi-treading.
     *  @returns A pointer to the new cloned allocated solid.
     */
    G4VSolid* Clone() const override;

    /**
     * Streams the object contents to an output stream.
     */
    std::ostream &StreamInfo(std::ostream& os) const override;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4TwistedTrap(__void__&);

    /**
     * Copy constructor and assignment operator.
     */
    G4TwistedTrap(const G4TwistedTrap& rhs);
    G4TwistedTrap& operator=(const G4TwistedTrap& rhs); 
} ;

#endif
