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
// $Id: $
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//    A enumeration to communicate from the G4BiasingProcessInterface
//    process to the current G4VBiasingOperator what case of biasing
//    has been applied in the PostStepDoIt(...)
//
//      ----------------G4BiasingAppliedCase ----------------
//
// Author: M.Verderi (LLR), November 2013
//
// --------------------------------------------------------------------

#ifndef G4BiasingAppliedCase_hh
#define G4BiasingAppliedCase_hh

enum G4BiasingAppliedCase
 {
   BAC_None,               // -- not under biasing
   BAC_NonPhysics,         // -- splitting, killing (not a physics process biasing)
   BAC_FinalState,         // -- physics process final state biasing only
   BAC_Occurence           // -- physics process occurence biasing; may come together with a final state biasing
 };

#endif
