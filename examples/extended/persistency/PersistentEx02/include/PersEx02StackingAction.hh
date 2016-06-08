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
// $Id: PersEx02StackingAction.hh,v 1.2.8.1 2001/06/28 19:07:24 gunter Exp $
// GEANT4 tag $Name:  $
//

#ifndef PersEx02StackingAction_H
#define PersEx02StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"

class G4Track;

class PersEx02StackingAction : public G4UserStackingAction
{
  public:
    PersEx02StackingAction();
    virtual ~PersEx02StackingAction();

  public:
   virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
   virtual void NewStage();
   virtual void PrepareNewEvent();

};

#endif

