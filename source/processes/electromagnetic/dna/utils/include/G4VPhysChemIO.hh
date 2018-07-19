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
 * G4VPhysChemIO.hh
 *
 *  Created on: 3 févr. 2017
 *      Author: matkara
 */

#ifndef G4VPHYSCHEMIO_HH_
#define G4VPHYSCHEMIO_HH_

#include <fstream>
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Track;
class G4VAnalysisManager;

class G4VPhysChemIO
{
public:
  G4VPhysChemIO();
  virtual ~G4VPhysChemIO();

  virtual void InitializeMaster(){}
  virtual void InitializeThread(){}
  virtual void InitializeFile() = 0;
  
  virtual void NewRun() = 0;
  virtual void NewEvent() = 0;

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
                                    const G4Track* /*theIncomingTrack*/) = 0;

   /**
    * Same idea as the previous method but for solvated electron.
    * This method should be used by the physics model of the ElectronSolvatation
    * process.
    */
   virtual void CreateSolvatedElectron(const G4Track* /*theIncomingTrack*/,
                                       G4ThreeVector* finalPosition = 0) = 0;

  //============================================================================
  // FILE OPERATIONS
  //============================================================================

  /**
   * Tells the chemistry manager to write into a file
   * the position and electronic state of the water molecule
   * and the position thermalized or not of the solvated electron
   */
  virtual void WriteInto(const G4String&, std::ios_base::openmode mode =
      std::ios_base::out) = 0;
  virtual void AddEmptyLineInOuputFile(){};

  /**
   * Close the file specified with WriteInto
   */
  virtual void CloseFile() = 0;
};

#endif // G4PHYSCHEMIO_HH_
