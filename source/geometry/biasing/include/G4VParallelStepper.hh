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
// $Id: G4VParallelStepper.hh,v 1.5 2002-09-02 13:25:26 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VParallelStepper
//
// Class description:
//
// This interface is used internally by importance sampling and scoring
// in a "parallel" geometry.
// It abstracts a parallel stepper which servs to comunicate
// information about a G4PStep between the "transportation"
// (done by G4ParallelTransport) of a track in a "parallel" geometry
// and a scorer for that geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VParallelStepper_hh
#define G4VParallelStepper_hh G4VParallelStepper_hh

#include "G4PStep.hh"

class G4VParallelStepper
{

public:  // with description

  virtual ~G4VParallelStepper() {}

  virtual G4PStep GetPStep() const = 0;
    // get the current G4PStep

  virtual void Init(const G4GeometryCell &agCell) = 0;
    // initialise the parallel stepper and the G4PStep
    // pre and post G4GeometryCell of the step are set equal

  virtual void Update(const G4GeometryCell &agCell) = 0;
    // to be called when crossing a boundary of the 
    // "parallel" geometry to update the G4PStep

  virtual void UnSetCrossBoundary() = 0;
    // to be called to unset the fCrossBoundary member of the G4PStep

};

#endif
