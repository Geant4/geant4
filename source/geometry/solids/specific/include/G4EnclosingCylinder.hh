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
// G4EnclosingCylinder
//
// Class description:
//
// A utility class for quickly deciding if a point is clearly outside a
// polyhedra or polycone or deciding if a trajectory is clearly going to
// miss those shapes.

// Author: David C. Williams (UCSC), 1998 - Created
// --------------------------------------------------------------------
#ifndef G4ENCLOSINGCYLINDER_HH
#define G4ENCLOSINGCYLINDER_HH

#include "G4Types.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4ReduciblePolygon;

/**
 * @brief G4EnclosingCylinder is a utility class defining an envelope for
 * quickly deciding if a point is clearly outside a polyhedra or polycone or
 * deciding if a trajectory is clearly going to miss those shapes.
 */

class G4EnclosingCylinder
{
  public:

    /**
     * Constructs the envelope, given its parameters.
     *  @param[in] rz Pointer to the polygon structure.
     *  @param[in] phiIsOpen Boolean flag to indicate if it is a section in Phi.
     *  @param[in] startPhi Starting Phi angle.
     *  @param[in] totalPhi Total Phi angle of the section.
     */
    G4EnclosingCylinder( const G4ReduciblePolygon* rz,
                               G4bool phiIsOpen, 
                               G4double startPhi, G4double totalPhi );

    /**
     * Default Destructor.
     */
    ~G4EnclosingCylinder() = default;
  
    /**
     * Decides very rapidly if the point 'p' is outside the cylinder.
     *  @returns If not certain to be outside, return false.
     */
    G4bool MustBeOutside( const G4ThreeVector& p ) const;

    /**
     * Decides very rapidly if the trajectory is going to miss the cylinder.
     *  @returns If not certain to miss, return false.
     */
    G4bool ShouldMiss( const G4ThreeVector& p, const G4ThreeVector& v ) const;

    /**
     * Fake default constructor for usage restricted to direct object
     * persistency for clients requiring preallocation of memory for
     * persistifiable objects.
     */
    G4EnclosingCylinder(__void__&);

  private:

    G4double radius;    // radius of our cylinder
    G4double zLo, zHi;  // z extent
  
    G4bool    phiIsOpen; // true if there is a phi segment
    G4double  startPhi,  // for isPhiOpen==true, starting of phi segment
              totalPhi;  // for isPhiOpen==true, size of phi segment

    G4double rx1=0.0, ry1=0.0,
             dx1=0.0, dy1=0.0;
    G4double rx2=0.0, ry2=0.0,
             dx2=0.0, dy2=0.0;
     
    G4bool   concave;  // true, if x/y cross section is concave
};

#endif
