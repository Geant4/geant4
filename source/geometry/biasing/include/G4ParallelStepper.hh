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
// $Id: G4ParallelStepper.hh,v 1.2 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelStepper
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelStepper_hh
#define G4ParallelStepper_hh

#include "G4VParallelStepper.hh"
#include "G4PStep.hh"

class G4ParallelStepper : public G4VParallelStepper
{

public:  // with description

  G4ParallelStepper();
  ~G4ParallelStepper();
  G4ParallelStepper(const G4ParallelStepper &);
  G4ParallelStepper &operator=(const G4ParallelStepper &);

  G4PStep GetPStep() const {return *fPStep;}
  void Init(const G4PTouchableKey &aptk);
  void Update(const G4PTouchableKey &aptk);
  void UnSetCrossBoundary();

private:

  void Error(const G4String &m);

private:

  G4PStep *fPStep;
};

#endif
