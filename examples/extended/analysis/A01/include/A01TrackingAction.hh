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
// $Id: A01TrackingAction.hh,v 1.3 2002-12-13 11:34:29 gunter Exp $
// --------------------------------------------------------------
//

#ifndef A01TrackingAction_h
#define A01TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class A01TrackingAction : public G4UserTrackingAction
{
  public:
    A01TrackingAction();
    virtual ~A01TrackingAction();

    virtual void PreUserTrackingAction(const G4Track*);
};

#endif
