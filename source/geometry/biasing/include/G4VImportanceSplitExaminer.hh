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
// $Id: G4VImportanceSplitExaminer.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
