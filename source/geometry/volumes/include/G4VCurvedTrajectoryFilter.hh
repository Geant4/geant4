//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VCurvedTrajectoryFilter.hh,v 1.4 2003-02-06 18:53:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "g4std/vector"
#include "G4TransportationManager.hh"

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

  G4std::vector<G4ThreeVector>* GimmeThePointsAndForgetThem();
  
protected:

  G4std::vector<G4ThreeVector>* fpFilteredPoints;
};

#endif  /* End of ifndef G4VCurvedTrajectoryFilter_hh */
