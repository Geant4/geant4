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
 * G4PhysChemIO.hh
 *
 *  Created on: 3 févr. 2017
 *      Author: matkara
 */
#ifndef G4PHYSCHEMIO_HH_
#define G4PHYSCHEMIO_HH_

#include "G4VPhysChemIO.hh"

//------------------------------------------------------------------------------
namespace G4PhysChemIO{
  
class FormattedText: public G4VPhysChemIO
{
public:
  FormattedText();
  virtual ~FormattedText();
  
  virtual void InitializeMaster(){}
  virtual void InitializeThread(){}
  virtual void InitializeFile();
  
  virtual void NewRun(){}
  virtual void NewEvent(){}
  
  /**
   * When DNA physics model create a water molecule, you'll get a notification
   * through this method.
   * The ElectronicModification is a flag telling whether the molecule
   * is ionized or excited, the electronic level is calculated by the
   * model and the IncomingTrack is the track responsible for the creation
   * of this molecule (electron, proton...)
   */
  virtual void CreateWaterMolecule(G4int electronicModif,
                                   G4int /*electronicLevel*/,
                                   G4double energy,
                                   const G4Track* /*theIncomingTrack*/);
  
  /**
   * Same idea as the previous method but for solvated electron.
   * This method should be used by the physics model of the ElectronSolvatation
   * process.
   */
  virtual void CreateSolvatedElectron(const G4Track* /*theIncomingTrack*/,
                                      G4ThreeVector* finalPosition = 0);
  
  //============================================================================
  // FILE OPERATIONS
  //============================================================================
  
  /**
   * Tells the chemistry manager to write into a file
   * the position and electronic state of the water molecule
   * and the position thermalized or not of the solvated electron
   */
  virtual void WriteInto(const G4String&,
                         std::ios_base::openmode mode = std::ios_base::out);
  virtual void AddEmptyLineInOuputFile();
  
  /**
   * Close the file specified with WriteInto
   */
  virtual void CloseFile();
  
protected:
  G4int fRunID; // unused
  G4int fEventID; // unused
  G4bool fFileInitialized;
  std::ofstream fOfstream;
};

//------------------------------------------------------------------------------

class G4Analysis: public G4VPhysChemIO
{
public:
  G4Analysis(G4VAnalysisManager*);
  virtual ~G4Analysis();
  
  virtual void InitializeMaster(){}
  virtual void InitializeThread(){}
  virtual void InitializeFile();
  
  virtual void NewRun(){}
  virtual void NewEvent(){}
  
  /**
   * Method used by DNA physics model to create a water molecule.
   * The ElectronicModification is a flag telling wheter the molecule
   * is ionized or excited, the electronic level is calculated by the
   * model and the IncomingTrack is the track responsible for the creation
   * of this molecule, for instance an electron.
   */
  virtual void CreateWaterMolecule(G4int electronicModif,
                                   G4int /*electronicLevel*/,
                                   G4double energy,
                                   const G4Track* /*theIncomingTrack*/);
  
  /**
   * Same idea as the previous method but for solvated electron.
   * This method should be used by the physics model of the ElectronSolvatation
   * process.
   */
  virtual void CreateSolvatedElectron(const G4Track* /*theIncomingTrack*/,
                                      G4ThreeVector* finalPosition = 0);
  
  //============================================================================
  // FILE OPERATIONS
  //============================================================================
  
  /**
   * Tells the chemMan to write into a file
   * the position and electronic state of the water molecule
   * and the position thermalized or not of the solvated electron
   */
  virtual void WriteInto(const G4String&, std::ios_base::openmode mode =
                         std::ios_base::out);
  virtual void AddEmptyLineInOuputFile(){}
  
  /**
   * Close the file specified with WriteInto
   */
  virtual void CloseFile();
  
protected:
  G4VAnalysisManager* fpAnalysisManager;
  int fNtupleID;
  G4bool fFileInitialized;
};
  
}

#endif // G4PHYSCHEMIO_HH_
