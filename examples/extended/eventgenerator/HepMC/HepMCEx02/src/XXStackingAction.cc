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
//
//   XXStackingAction.cc
//   $Id: XXStackingAction.cc,v 1.1 2002-04-29 20:45:04 asaim Exp $
//
// ====================================================================

#include "XXStackingAction.hh"
#include "G4Track.hh"

////////////////////////////////////
XXStackingAction::XXStackingAction()
////////////////////////////////////
{ 
}

/////////////////////////////////////
XXStackingAction::~XXStackingAction()
/////////////////////////////////////
{ 
}

/////////////////////////////////////////////////////////////
G4ClassificationOfNewTrack XXStackingAction::ClassifyNewTrack
                             (const G4Track* atrack)
/////////////////////////////////////////////////////////////
{
  G4ClassificationOfNewTrack classification= fWaiting;

  // apply eta cut
  //const G4double eta_max= 3.;
  //if(atrack->GetParentID() == 0) {
  //  G4ThreeVector vmom= atrack-> GetMomentum();
  //  if( abs(vmom.eta()) > eta_max ) {
  //    classification= fKill;
  //  }
  //}

  return classification;
}

/////////////////////////////////
void XXStackingAction::NewStage()
/////////////////////////////////
{
}

////////////////////////////////////////
void XXStackingAction::PrepareNewEvent()
////////////////////////////////////////
{ 
}

