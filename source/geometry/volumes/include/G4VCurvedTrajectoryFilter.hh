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
// $Id: G4VCurvedTrajectoryFilter.hh,v 1.1 2002-11-09 00:25:08 jacek Exp $
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
// IMPORTANT: This class heap allocates vectors of auxiliary points,
// which it does not delete. The vectors must find their way to a
// subclass of G4VTrajectoryPoint, which must take responsibility for
// deleting them.

// History
// - Modified:
// - First version: Oct 30, 2002  Jacek Generowicz


// Writing the concrete implementation in the base class initially,
// remember to move implementation to a derived class once happy with
// structure. (jacek 30/10/2002)

#ifndef G4VCurvedTrajectoryFilter_hh
#define G4VCurvedTrajectoryFilter_hh

#include "G4ThreeVector.hh"
#include "g4std/vector"
#include "G4TransportationManager.hh"

class G4VCurvedTrajectoryFilter {
public:

  G4VCurvedTrajectoryFilter(){};
  virtual ~G4VCurvedTrajectoryFilter(){};

  // Probably do not want these objects to be copied, so make the copy
  // constructor prvate (jacek 30/10/2002)
     	 
  
  // Each segment stores the auxiliary points of a single step.
  void CreateNewTrajectorySegment( );

  // Submit intermediate points for the filter to consider keeping or
  // rejecting
  void TakeIntermediatePoint( G4ThreeVector newPoint );
  
  G4std::vector<G4ThreeVector>* GimmeThePointsAndForgetThem();
  
private:
  G4std::vector<G4ThreeVector>* fpFilteredPoints;
};

#endif  /* End of ifndef G4VCurvedTrajectoryFilter_hh */
