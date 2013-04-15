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
// $Id$
//
// Author: Ivana Hrivnacova, 10/04/2013  (ivana@ipno.in2p3.fr)
//
// A default worker initialization using G4VUserApplication

#ifndef G4WorkerInitialization_hh
#define G4WorkerInitialization_hh

#include "G4VUserWorkerInitialization.hh"

class G4VUserApplication;

class G4WorkerInitialization : public G4VUserWorkerInitialization {

  public: // with description
    G4WorkerInitialization(G4VUserApplication* userApplication);
    virtual ~G4WorkerInitialization();

    virtual void WorkerStart() const;
    //    Default implementation of the instantiation of user actions
    //on worker with use of G4VUserApplication interface.

  private:
    G4VUserApplication* fUserApplication;
};

#endif //G4WorkerInitialization_hh

