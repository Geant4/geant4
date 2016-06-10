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
// $Id: P6DPhysicsList.hh 78019 2013-12-02 15:52:19Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/P6DPhysicsList.hh
/// \brief Definition of the P6DPhysicsList class

#ifndef P6DPhysicsList_h
#define P6DPhysicsList_h 1

#include "G4VUserPhysicsList.hh"

/// The physics list class with Pythia6 decayer

class P6DPhysicsList: public G4VUserPhysicsList
{
  public:

    P6DPhysicsList();
    ~P6DPhysicsList();

    // Construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();
 
    virtual void SetCuts();
   
  private:

    // these methods Construct physics processes and register them
    void ConstructDecay();
    void ConstructEM();
};

// ----------------------------------------------------------------------------

#endif
