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
// $Id: G4SLRawValues.hh,v 1.1 2002-07-10 15:51:04 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4SLRawValues
//
// Class description:
// This struct holds the raw values for estimators
// vased on the track length. It's used for standard scoreing.
// 
// The Weight, Energy, and Velocity is taken from the pre step point.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4SLRawValues_hh
#define G4SLRawValues_hh G4SLRawValues_hh

#include "globals.hh"
#include "G4Step.hh"

struct G4SLRawValues {
  G4SLRawValues(const G4Step &aStep) :
    fStepLength(aStep.GetStepLength()),
    fPre_Weight(aStep.GetPreStepPoint()->GetWeight()),
    fPre_Energy(aStep.GetPreStepPoint()->GetKineticEnergy()),
    fPre_Velocity(aStep.GetPreStepPoint()->GetVelocity())
  {}
 
  G4double fStepLength;
  G4double fPre_Weight;
  G4double fPre_Energy;
  G4double fPre_Velocity;
};


#endif
