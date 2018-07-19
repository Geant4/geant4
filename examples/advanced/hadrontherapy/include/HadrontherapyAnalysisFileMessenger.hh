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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef HadrontherapyAnalysisFileMessenger_h
#define HadrontherapyAnalysisFileMessenger_h 1


#include "G4UImessenger.hh"
#include "globals.hh"

class HadrontherapyAnalysisManager; ///< Provides SetanalysisFileName()
class G4UIcmdWithAString; 
class G4UIcmdWithABool; 

/**
 * A messenger object of this class is created by the AnalysisManager.
 * The point of a messenger is to connect the G4UI with the simulation
 * functionality.
 * The messenger needs to contain a command object and to have SetValue
 * method that is called once a command is set.
 * 
 * @see HadrontherapyAnalysisManager
 */
class HadrontherapyAnalysisFileMessenger: public G4UImessenger
{
public:
  HadrontherapyAnalysisFileMessenger(HadrontherapyAnalysisManager*);
  ~HadrontherapyAnalysisFileMessenger();

  /**   
   * Called when new command given.
   * @param command is a pointer to the given command object
   * @param newValue holds the argument given as a G4String
   * @return is void   
   */     
  void SetNewValue(G4UIcommand* command, G4String newValue);
  HadrontherapyAnalysisManager* AnalysisManager; ///< handle to AnalysisManager
    
private:

  /**   
   * G4 user interface command (that takes a string argument) object
   * Constructor requires command name and messenger class(this).
   */ 
  G4UIcmdWithABool *LetCmd, *secondaryCmd; 
  G4UIcmdWithAString *DoseMatrixCmd;

};

#endif
