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

// The extra MPI worker class defines actions on an extra MPI worker
// which does not perform event processing.
// It calls the BeginOfRunAction() and EndOfRunAction().
// Currently used only for ntuple merging.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#include "G4MPIextraWorker.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4ios.hh"

void G4MPIextraWorker::BeamOn()
{
  G4Run dummyRun;

  G4cout << "G4MPIextraWorker: call BeginOfRunAction()" << G4endl;
  fRunAction->BeginOfRunAction(&dummyRun);

  G4cout << "G4MPIextraWorker: call EndOfRunAction()" << G4endl;
  fRunAction->EndOfRunAction(&dummyRun);
}
