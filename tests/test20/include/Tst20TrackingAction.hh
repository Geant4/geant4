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
// $Id: Tst20TrackingAction.hh,v 1.3 2001-07-11 10:10:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//-----------------------------------------------------------------------------

#ifndef Tst20TrackingAction_h
#define Tst20TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class Tst20TrackingAction : public G4UserTrackingAction {

  public:
    Tst20TrackingAction(){};
    virtual ~Tst20TrackingAction(){};
   
    virtual void PreUserTrackingAction(const G4Track* aTrack);

};

#endif
