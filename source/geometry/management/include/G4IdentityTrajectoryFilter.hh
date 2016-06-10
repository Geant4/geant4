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
// $Id: G4IdentityTrajectoryFilter.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
//
// - First version: Nov 19, 2002  Jacek Generowicz
// ------------------------------------------------------------------------
#ifndef G4IdentityTrajectoryFilter_hh
#define G4IdentityTrajectoryFilter_hh

#include "G4VCurvedTrajectoryFilter.hh"
#include "G4ThreeVector.hh"

class G4IdentityTrajectoryFilter : public G4VCurvedTrajectoryFilter
{

public:  // with description

  G4IdentityTrajectoryFilter();
  virtual ~G4IdentityTrajectoryFilter();

    // Probably do not want these objects to be copied,
    // so make the copy constructor private

  void TakeIntermediatePoint( G4ThreeVector newPoint );
    // Submit intermediate points for the filter
    // to consider keeping or rejecting
};

#endif  /* End of ifndef G4IdentityTrajectoryFilter_hh */
