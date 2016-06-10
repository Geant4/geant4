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
// Class Descrition:
//   This class implements TBB specific worker initialization.
//   It is a sub-class of G4UserWorkerThreadInitialization.
//   Its role is to instantiate a tbbWorkerRunManager to be used by
//   tbb tasks. 
//
// Equivalent in traditional MT:
//   G4UserWorkerThreadInitialization
//
// History:
//    Oct 31st, 2013 A. Dotti - First Implementation

#ifndef TBBUSERWORKERINITIALIZATION_HH
#define TBBUSERWORKERINITIALIZATION_HH

#include "G4UserWorkerThreadInitialization.hh"

class tbbUserWorkerInitialization : public G4UserWorkerThreadInitialization {
public:
    virtual ~tbbUserWorkerInitialization();
    virtual G4WorkerRunManager* CreateWorkerRunManager() const;
    virtual void JoinWorker(G4Thread*);
    virtual G4Thread* CreateAndStartWorker(G4WorkerThread*);
};

#endif //TBBUSERWORKERINITIALIZATION_HH
