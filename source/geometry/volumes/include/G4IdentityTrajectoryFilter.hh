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
// $Id: G4IdentityTrajectoryFilter.hh,v 1.1 2002-11-19 15:23:24 jacek Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class GIdentityTrajectoryFilter
//
// Class description:
//
// Implements a trajectory point filter which accepts all points submitted to it.
//
// IMPORTANT: The base class heap allocates vectors of auxiliary
// points, which it does not delete. The vectors must find their way
// to a subclass of G4VTrajectoryPoint, which must take responsibility
// for deleting them.

// History
// - Modified:
// - First version: Nov 19, 2002  Jacek Generowicz

#ifndef G4IdentityTrajectoryFilter_hh
#define G4IdentityTrajectoryFilter_hh

//#include "g4std/vector"
#include "G4VCurvedTrajectoryFilter.hh"

class G4IdentityTrajectoryFilter : public G4VCurvedTrajectoryFilter {
public:

  G4IdentityTrajectoryFilter(){};
  virtual ~G4IdentityTrajectoryFilter(){};

  // Probably do not want these objects to be copied, so make the copy
  // constructor prvate (jacek 30/10/2002)

  // Submit intermediate points for the filter to consider keeping or
  // rejecting.
  void TakeIntermediatePoint( G4ThreeVector newPoint );
};

#endif  /* End of ifndef G4IdentityTrajectoryFilter_hh */
