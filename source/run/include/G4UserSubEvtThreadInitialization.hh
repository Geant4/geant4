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
// Author: Makoto Asai (JLAB)
//
// class description:
//
//      This class is used for multi-threaded Geant4.
//      It encapsulates the mechanism of starting/stopping threads.

#ifndef G4UserSubEvtThreadInitialization_hh
#define G4UserSubEvtThreadInitialization_hh

#include "G4Threading.hh"
#include "G4UserTaskThreadInitialization.hh"

class G4UserSubEvtThreadInitialization : public G4UserTaskThreadInitialization
{
  public:  // with description
    G4UserSubEvtThreadInitialization() = default;
    ~G4UserSubEvtThreadInitialization() override = default;

    // Called by StartThread function to create a run-manager implementing worker
    // behvior. User should re-implemtn this function in derived class to
    // instantiate his/her user-defined WorkerRunManager. By default this method
    // instantiates G4WorkerSubEvtRunManager object.
    G4WorkerRunManager* CreateWorkerRunManager() const override;
};

#endif  // G4UserSubEvtThreadInitialization_hh
