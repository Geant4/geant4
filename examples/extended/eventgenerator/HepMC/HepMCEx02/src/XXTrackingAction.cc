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

// ====================================================================
//    XXTrackingAction.cc
//    $Id: XXTrackingAction.cc,v 1.1 2002-04-29 20:45:05 asaim Exp $
//
// ====================================================================

#include "XXTrackingAction.hh"
#include "G4TrackingManager.hh"
#include "G4Track.hh"

///////////////////////////////////////////////////////////////////
void XXTrackingAction::PreUserTrackingAction(const G4Track* atrack)
///////////////////////////////////////////////////////////////////
{
  // create trajectory only for primaries
  //if( atrack-> GetParentID()==0 ) { 
  //  fpTrackingManager-> SetStoreTrajectory(true); 
  //} else { 
  //  fpTrackingManager-> SetStoreTrajectory(false); 
  //}
}
