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
// $Id: G4ForceCondition.hh 68795 2013-04-05 13:24:46Z gcosmo $
//
//
//---------------------------------------------------------------
//
// G4ForceCondition  
//
// Class Description:
//   This enumaration specifies possible conditions the three
//   types of DoIt can be assinged by physics processes.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4ForceCondition_h
#define G4ForceCondition_h 1

/////////////////////
enum G4ForceCondition  
/////////////////////
{
  InActivated,       
  // This PostStepDoit is inactivated by a user
  Forced,            
    // This PostStepDoIt is forced to invoke if particle is not a state of StopAndKill.              
  NotForced,         
    // This PostStepDoIt is not forced to invoke.
  Conditionally,     
    // This PostStepDoIt is forced to invoke only when corresponding
    // AlongStepDoIt limits the Step.
  ExclusivelyForced, 
    // Only this PostStepDoIt (or AtRestDoIt) is exclusively forced 
    // to invoke - all other DoIt including AlongStepDoIts are ignored.
  StronglyForced
   // This PostStepDoIt is really forced to invoke, anyway.
};

#endif


