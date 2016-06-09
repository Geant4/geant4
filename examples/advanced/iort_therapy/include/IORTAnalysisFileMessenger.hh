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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef IORTAnalysisFileMessenger_h
#define IORTAnalysisFileMessenger_h 1


#include "G4UImessenger.hh"
#include "globals.hh"

class IORTAnalysisManager; ///< Provides SetanalysisFileName()
class G4UIcmdWithAString; 
class G4UIcmdWithABool; 

/**
 * A messenger object of this class is created by the AnalysisManager.
 * The point of a messenger is to connect the G4UI with the simulation
 * functionality.
 * The messenger needs to contain a command object and to have SetValue
 * method that is called once a command is set.
 * 
 * @see IORTAnalysisManager
 */
class IORTAnalysisFileMessenger: public G4UImessenger
{
public:
  IORTAnalysisFileMessenger(IORTAnalysisManager*);
  ~IORTAnalysisFileMessenger();

  /**   
   * Called when new command given.
   * @param command is a pointer to the given command object
   * @param newValue holds the argument given as a G4String
   * @return is void   
   */     
  void SetNewValue(G4UIcommand* command, G4String newValue);
    
private:
  IORTAnalysisManager* AnalysisManager; ///< handle to AnalysisManager

  /**   
   * G4 user interface command (that takes a string argument) object
   * Constructor requires command name and messenger class(this).
   */ 
  G4UIcmdWithABool *secondaryCmd; 
  G4UIcmdWithAString *DoseMatrixCmd;
#ifdef G4ANALYSIS_USE_ROOT
  G4UIcmdWithAString *FileNameCmd;
#endif
};

#endif
