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
// $Id: G4ImportanceSplitExaminer.hh,v 1.8 2006/06/29 18:16:08 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4ImportanceSplitExaminer
//
// Class description:
//
// This class is used internally by importance sampling.
// Implementation of a importance split examiner 
// (see G4VImportanceSplitExaminer).
// This implementation is used for sampling in a "parallel"
// geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ImportanceSplitExaminer_hh
#define G4ImportanceSplitExaminer_hh G4ImportanceSplitExaminer_hh

#include "G4VImportanceSplitExaminer.hh"

class G4VImportanceAlgorithm;
class G4VParallelStepper;
class G4VIStore;

class G4ImportanceSplitExaminer: public G4VImportanceSplitExaminer
{

public:  // with description

  G4ImportanceSplitExaminer(const G4VImportanceAlgorithm &aIalg,
		      const G4VParallelStepper &astepper,
		      const G4VIStore &istore);

  virtual ~G4ImportanceSplitExaminer();

  virtual G4Nsplit_Weight Examine(G4double w) const; 
    // Get  G4Nsplit_Weight for a given mother track weight.
  

private:

  G4ImportanceSplitExaminer(const G4ImportanceSplitExaminer &);
  G4ImportanceSplitExaminer &operator=(const G4ImportanceSplitExaminer &);

private:

  const G4VImportanceAlgorithm &fIalgorithm;
  const G4VParallelStepper &fPStepper;
  const G4VIStore &fIStore;
};


#endif
