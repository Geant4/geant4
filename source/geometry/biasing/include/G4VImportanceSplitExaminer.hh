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
// $Id: G4VImportanceSplitExaminer.hh,v 1.4 2002-10-14 12:36:01 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VImportanceSplitExaminer
//
// Class description:
//
// This interface is used internally by importance sampling. 
// It delivers G4Nsplit_Weight according to a track weight. 
// An implementation of the interface decides how to obtain 
// remaining necessary information about the ratio of importances
// in the pre and post "cell". 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VImportanceSplitExaminer_hh
#define G4VImportanceSplitExaminer_hh G4VImportanceSplitExaminer_hh

#include "globals.hh"
#include "G4Nsplit_Weight.hh"

class G4VImportanceSplitExaminer
{

public:  // with description
  G4VImportanceSplitExaminer();
  virtual ~G4VImportanceSplitExaminer();
  virtual G4Nsplit_Weight Examine(G4double w) const = 0; 
    // Get  G4Nsplit_Weight for a given mother track weight.
};

#endif
