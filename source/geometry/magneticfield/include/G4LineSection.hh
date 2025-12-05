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
// G4LineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// Author: John Apostolakis (CERN), 1999
// --------------------------------------------------------------------

#ifndef G4LineSection_hh
#define G4LineSection_hh

#include "G4Types.hh" 
#include "G4ThreeVector.hh"

/**
 * @brief G4LineSection is a utility class that calculates the distance
 * of a point from a line section.
 */

class G4LineSection
{
  public:

    /**
     * Constructor for G4LineSection.
     *  @param[in] PntA Coordinates of point A defining the line.
     *  @param[in] PntB Coordinates of point B defining the line.
     */
    G4LineSection( const G4ThreeVector& PntA,
                   const G4ThreeVector& PntB );

    /**
     * Default Destructor.
     */
    ~G4LineSection() = default;

    /**
     * Returns the distance of point 'OtherPnt' from the line.
     */
    G4double Dist( const G4ThreeVector& OtherPnt ) const;

    /**
     * Returns the distance squared.
     */
    inline G4double GetABdistanceSq() const;

    /**
     * Defines line and returns the distance of point 'OtherPnt' from it.
     */
    inline static G4double Distline( const G4ThreeVector& OtherPnt, 
                                     const G4ThreeVector& LinePntA, 
                                     const G4ThreeVector& LinePntB );
  private:

    G4ThreeVector EndpointA;
    G4ThreeVector VecAtoB;
    G4double fABdistanceSq = 0.0;
};

// Inline methods implementations

inline
G4double G4LineSection::GetABdistanceSq() const
{
  return fABdistanceSq;
}

inline
G4double G4LineSection::Distline( const G4ThreeVector& OtherPnt, 
                                  const G4ThreeVector& LinePntA, 
                                  const G4ThreeVector& LinePntB )
{
  G4LineSection LineAB( LinePntA, LinePntB );  // Line from A to B
  return LineAB.Dist( OtherPnt );
}

#endif
