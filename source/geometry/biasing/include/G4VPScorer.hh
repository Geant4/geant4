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
// $Id: G4VPScorer.hh,v 1.3 2002-04-10 13:13:07 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VPScorer
//
// Class description:
//
// This interface is known by scoring. Implementation of 
// scorers should inherit from it. 
// A scorer is applicable to a particle type.
// The scorer (inheriting from this interface) is informed
// about every step of a particle of a the applicable type.
// The information consists of the G4Step and the G4PStep
// (see G4PStep.hh).

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VPScorer_hh
#define G4VPScorer_hh G4VPScorer_hh
  
class G4Step;
class G4PStep;
  
class G4VPScorer
{

public:  // with description

  virtual ~G4VPScorer() {}
  virtual void Score(const G4Step &step, const G4PStep &pstep) = 0;
    // perform scoring for the G4Step and G4PStep
};
  
#endif





