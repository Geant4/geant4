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
// $Id: G4VScorer.hh,v 1.2 2006/06/29 18:16:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
