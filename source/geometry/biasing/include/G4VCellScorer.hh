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
// $Id: G4VCellScorer.hh,v 1.1 2002-08-13 10:03:29 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VCellScorer
//
// Class description:
//
// This is an interface for an object which does scoring for one cell
// of a geometry.
// The cell may be a physical volume or replica in the mass or a 
// parallel geometry.
// 
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4VCellScorer_hh
#define G4VCellScorer_hh G4VCellScorer_hh

class G4Step;
class G4PTouchableKey;

class G4VCellScorer {
public:  
  virtual  ~G4VCellScorer(){}
  virtual void ScoreAnExitingStep(const G4Step &aStep, 
				  G4PTouchableKey ptk) = 0;
  virtual void ScoreAnEnteringStep(const G4Step &aStep, 
				   G4PTouchableKey ptk) = 0;
  virtual void ScoreAnInVolumeStep(const G4Step &aStep, 
				   G4PTouchableKey ptk) = 0;
};

#endif
