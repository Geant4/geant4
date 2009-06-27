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
// HadrontherapyStepMaxMessenger.hh
//
// See more at: http://workgroup.lngs.infn.it/geant4lns/
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
/**
* @file HadrontherapyAnalysisFileMessenger.hh. 
* @author Gillis Danielsen, HIP 
* @date 6-10-09  
*/ 

#ifndef HadrontherapyAnalysisFileMessenger_h
#define HadrontherapyAnalysisFileMessenger_h 1

#ifdef ANALYSIS_USE

#include "G4UImessenger.hh"
#include "globals.hh"

class HadrontherapyAnalysisManager; ///< Provides SetAnalysisFileName()
class G4UIcmdWithAString; ///< Provides possibility to add commands

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
    
  private:
    HadrontherapyAnalysisManager* AnalysisManager; ///< handle to AnalysisManager

	/**   
	* G4 user interface command (that takes a string argument) object
	* Constructor requires command name and messenger class(this).
	*/ 
    G4UIcmdWithAString* FileNameCmd;
};

#endif
#endif
