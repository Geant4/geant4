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
// $Id: G4ImportanceSplitExaminer.hh,v 1.5 2002-10-14 12:36:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ImportanceFinder.hh"

class G4VImportanceAlgorithm;
class G4VParallelStepper;
class G4VIStore;

class G4ImportanceSplitExaminer: public G4VImportanceSplitExaminer
{

public:  // with description

  G4ImportanceSplitExaminer(const G4VImportanceAlgorithm &aIalg,
		      const G4VParallelStepper &astepper,
		      const G4VIStore &istore);
    // initialisation and construct G4ImportanceFinder

  virtual ~G4ImportanceSplitExaminer();
    // delete G4ImportanceFinder

  virtual G4Nsplit_Weight Examine(G4double w) const; 
    // Get  G4Nsplit_Weight for a given mother track weight.
  

private:

  G4ImportanceSplitExaminer(const G4ImportanceSplitExaminer &);
  G4ImportanceSplitExaminer &operator=(const G4ImportanceSplitExaminer &);

private:

  const G4VImportanceAlgorithm &fIalgorithm;
  const G4VParallelStepper &fPStepper;
  G4ImportanceFinder fIfinder;
};


#endif
