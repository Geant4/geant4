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
//
// $Id: G4VParallelStepper.hh,v 1.8 2006/06/29 18:16:47 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4VParallelStepper
//
// Class description:
//
// This interface is used internally by importance sampling and scoring
// in a "parallel" geometry.
// It abstracts a parallel stepper which servs to comunicate
// information about a G4GeometryCellStep between the "transportation"
// (done by G4ParallelTransport) of a track in a "parallel" geometry
// and a scorer for that geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VParallelStepper_hh
#define G4VParallelStepper_hh G4VParallelStepper_hh

#include "G4GeometryCellStep.hh"

class G4VParallelStepper
{

public:  // with description

  G4VParallelStepper();
  virtual ~G4VParallelStepper();

  virtual G4GeometryCellStep GetPStep() const = 0;
    // get the current G4GeometryCellStep

  virtual void Init(const G4GeometryCell &agCell) = 0;
    // initialise the parallel stepper and the G4GeometryCellStep
    // pre and post G4GeometryCell of the step are set equal

  virtual void Update(const G4GeometryCell &agCell) = 0;
    // to be called when crossing a boundary of the 
    // "parallel" geometry to update the G4GeometryCellStep

  virtual void UnSetCrossBoundary() = 0;
    // to be called to unset the fCrossBoundary member of the G4GeometryCellStep

};

#endif
