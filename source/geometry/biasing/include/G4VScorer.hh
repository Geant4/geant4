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
// $Id: G4VScorer.hh,v 1.1 2002-10-28 09:53:50 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VScorer
//
// Class description:
//
// This interface is known by scoring. Implementation of 
// scorers should inherit from it. 
// A scorer is applicable to a particle type.
// The scorer (inheriting from this interface) is informed
// about every step of a particle of the applicable type.
// The information consists of the G4Step and the G4GeometryCellStep
// (see G4GeometryCellStep.hh).

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VScorer_hh
#define G4VScorer_hh G4VScorer_hh
  
class G4Step;
class G4GeometryCellStep;
  
class G4VScorer
{

public:  // with description

  G4VScorer();

  virtual ~G4VScorer();

  virtual void Score(const G4Step &step, const G4GeometryCellStep &gstep) = 0;
    // perform scoring for the G4Step and G4GeometryCellStep
};
  
#endif
