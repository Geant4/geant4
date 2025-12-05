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
// G4IntersectingCone
//
// Class description:
//
// Utility class which calculates the intersection
// of an arbitrary line with a fixed cone.

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4INTERSECTINGCONE_HH
#define G4INTERSECTINGCONE_HH

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

/**
 * @brief G4IntersectingCone is a utility class used to calculate the
 * intersection of an arbitrary line with a fixed cone.
 */

class G4IntersectingCone
{
  public:

    /**
     * Constructor given r,z values.
     *  @param[in] r r values.
     *  @param[in] z Z values.
     */
    G4IntersectingCone( const G4double r[2], const G4double z[2] );

    /**
     * Default Destructor.
     */
    ~G4IntersectingCone() = default;
  
    /**
     * Calculates the intersection of a line with the conical surface,
     * ignoring any Phi division.
     */
    G4int LineHitsCone( const G4ThreeVector& p, const G4ThreeVector& v,
                              G4double* s1, G4double* s2 );
  
    /**
     * Checks r or z extent, as appropriate, to see if the point is
     * possibly on the cone.
     */
    G4bool HitOn( const G4double r, const G4double z );
  
    /**
     * Accessors for R and Z bounds of side.
     */
    inline G4double RLo() const { return rLo; }
    inline G4double RHi() const { return rHi; }
    inline G4double ZLo() const { return zLo; }
    inline G4double ZHi() const { return zHi; }
  
    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4IntersectingCone(__void__&);

  private:

    /**
     * Calculating the intersection of a line with the conical surface.
     * Internal methods used by LineHitsCone().
     */
    G4int LineHitsCone1( const G4ThreeVector& p, const G4ThreeVector& v,
                               G4double* s1, G4double* s2 );
    G4int LineHitsCone2( const G4ThreeVector& p, const G4ThreeVector& v,
                               G4double* s1, G4double* s2 );
  private:

    /** Z, R bounds of side. */
    G4double zLo, zHi, rLo, rHi;

    /** True if cone is type 1. */
    G4bool type1 = false; // (std::fabs(z1-z2)>std::fabs(r1-r2))

    /** Cone radius parameters - type 1: r = A + B*z; type 2: z = A + B*r. */
    G4double A, B;
};

#endif
