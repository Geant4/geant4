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
// $Id: G4ParallelStepper.hh,v 1.5 2002-10-14 12:36:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelStepper
//
// Class description:
//
// Used internally by importance sampling and scoring in a "parallel"
// geometry. (See G4VParallelStepper.hh)

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelStepper_hh
#define G4ParallelStepper_hh G4ParallelStepper_hh

#include "G4VParallelStepper.hh"
#include "G4PStep.hh"

class G4ParallelStepper : public G4VParallelStepper
{

public:  // with description

  G4ParallelStepper();
    // imitilisation

  G4ParallelStepper(const G4ParallelStepper &);
    // create new G4PStep

  virtual ~G4ParallelStepper();
    // delete G4Pstep if created

  virtual G4PStep GetPStep() const;
    // get the current G4PStep
  
  virtual void Init(const G4GeometryCell &agCell);
    // initialise the parallel stepper and the G4PStep
    // pre and post G4GeometryCell of the step are set equal

  virtual void Update(const G4GeometryCell &agCell);
    // to be called when crossing a boundary of the 
    // "parallel" geometry to update the G4PStep

  virtual void UnSetCrossBoundary();
    // to be called to unset the fCrossBoundary member of the G4PStep


  G4ParallelStepper &operator=(const G4ParallelStepper &);
    // create new G4PStep

private:

  void Error(const G4String &m);

private:

  G4PStep *fPStep;
};



#endif
