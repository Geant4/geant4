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
// $Id: G4VParallelStepper.hh,v 1.2 2002-04-09 16:23:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VParallelStepper
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VParallelStepper_hh
#define G4VParallelStepper_hh

#include "G4PStep.hh"

class G4VParallelStepper
{

public:  // with description

  virtual ~G4VParallelStepper(){}
  virtual G4PStep GetPStep() const = 0;
  virtual void Init(const G4PTouchableKey &aptk) = 0;
  virtual void Update(const G4PTouchableKey &aptk) = 0;
  virtual void UnSetCrossBoundary() = 0;
};

#endif
