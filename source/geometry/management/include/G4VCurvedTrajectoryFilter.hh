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
//
// $Id: G4VCurvedTrajectoryFilter.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// class G4VCurevedTrajectoryFilter
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

// History
//
// - First version: Oct 30, 2002  Jacek Generowicz
// ------------------------------------------------------------------------
#ifndef G4VCurvedTrajectoryFilter_hh
#define G4VCurvedTrajectoryFilter_hh

#include "G4ThreeVector.hh"
#include <vector>

class G4VCurvedTrajectoryFilter
{

public:  // with description

  G4VCurvedTrajectoryFilter();
  virtual ~G4VCurvedTrajectoryFilter();

    // Probably do not want these objects to be copied,
    // so make the copy constructor private
  
  void CreateNewTrajectorySegment( );
    // Each segment stores the auxiliary points of a single step.

  virtual void TakeIntermediatePoint( G4ThreeVector newPoint ) = 0;
    // Submit intermediate points for the filter to consider keeping or
    // rejecting. Derived classes should implement the filtering algorithm
    // in this method.

  std::vector<G4ThreeVector>* GimmeThePointsAndForgetThem();
  
protected:

  std::vector<G4ThreeVector>* fpFilteredPoints;
};

#endif  /* End of ifndef G4VCurvedTrajectoryFilter_hh */
