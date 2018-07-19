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
// $Id: RE05TrackingAction.hh 98775 2016-08-09 14:30:39Z gcosmo $
//
/// \file RE05/include/RE05TrackingAction.hh
/// \brief Definition of the RE05TrackingAction class
//

#ifndef RE05TrackingAction_h
#define RE05TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class RE05TrackingAction : public G4UserTrackingAction {

  public:
    RE05TrackingAction() : G4UserTrackingAction() {}
    virtual ~RE05TrackingAction() {}
   
    virtual void PreUserTrackingAction(const G4Track*);
};

#endif
