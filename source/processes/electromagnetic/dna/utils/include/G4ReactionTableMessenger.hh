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
/*
 * G4ReactionTableMessenger.hh
 *
 *  Created on: Sep 14, 2015
 *      Author: mkaramit
 */

#ifndef SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_UTILS_INCLUDE_G4REACTIONTABLEMESSENGER_HH_
#define SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_UTILS_INCLUDE_G4REACTIONTABLEMESSENGER_HH_

#include <G4UImessenger.hh>

#include <memory>

class G4UIcmdWithAString;
class G4DNAMolecularReactionTable;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;

class G4ReactionTableMessenger : public G4UImessenger
{
public:
  G4ReactionTableMessenger(G4DNAMolecularReactionTable*);
  virtual ~G4ReactionTableMessenger();
  virtual void SetNewValue(G4UIcommand * command,G4String newValue);

protected:
  G4DNAMolecularReactionTable* fpTable;
  std::unique_ptr<G4UIcmdWithoutParameter> fpActivateReactionUI;
  G4UIcmdWithAString* fpAddReaction;
  G4UIcmdWithAString* fpNewDiffContReaction;
//  G4UIcmdWithAString* fpNewPartDiffContReactionByRadius;
//  G4UIcmdWithAString* fpNewPartDiffContReactionByReactionRate;
  G4UIcmdWithoutParameter* fpPrintTable;
};

#endif /* SOURCE_PROCESSES_ELECTROMAGNETIC_DNA_UTILS_INCLUDE_G4REACTIONTABLEMESSENGER_HH_ */
