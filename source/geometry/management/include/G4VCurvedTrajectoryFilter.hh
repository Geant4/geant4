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
// G4VCurvedTrajectoryFilter
//
// Class description:
//
// Decides which intermediate points on a curved trajectory merit
// being stored. Defines the compromise between accuracy of
// representation of the curved trajectory and memory use.
//
// Derived classes should implement the filtering algorithm in the
// method TakeIntermediatePoint().
//
// IMPORTANT: This class heap allocates vectors of auxiliary points,
// which it does not delete. The vectors must find their way to a
// subclass of G4VTrajectoryPoint, which must take responsibility for
// deleting them.

// Author: Jacek Generowicz (CERN), 30.10.2002
// ------------------------------------------------------------------------
#ifndef G4VCURVEDTRAJECTORYFILTER_HH
#define G4VCURVEDTRAJECTORYFILTER_HH

#include "G4ThreeVector.hh"
#include <vector>

/**
 * @brief G4VCurvedTrajectoryFilter defines a filter for deciding which
 * intermediate points on a curved trajectory merit being stored. It defines
 * the compromise between accuracy of representation of the curved trajectory
 * and memory use. Derived classes should implement the filtering algorithm.
 */

class G4VCurvedTrajectoryFilter
{
  public:

    /**
     * Default Constructor & Destructor.
     */
    G4VCurvedTrajectoryFilter() = default;
    virtual ~G4VCurvedTrajectoryFilter() = default;

    // Probably do not want these objects to be copied,
    // so make the copy constructor deleted in derived classes.
  
    /**
     * Each segment stores the auxiliary points of a single step.
     */
    void CreateNewTrajectorySegment();

    /**
     * Submits intermediate points for the filter to consider keeping or
     * rejecting. Derived classes should implement the filtering algorithm
     * in this method.
     */
    virtual void TakeIntermediatePoint( G4ThreeVector newPoint ) = 0; 

    /**
     * Returns the vector of points, transferring the ownership.
     */
    std::vector<G4ThreeVector>* GimmeThePointsAndForgetThem();
  
  protected:

    std::vector<G4ThreeVector>* fpFilteredPoints = nullptr;
};

#endif
