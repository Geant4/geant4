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
////
//---------------------------------------------------------------
//
// G4BuilderInterface.hh
//
// Class Description:
//      Provides the common interface to all types of builders.
//
// Creation date: 09.04.2016 adotti
// Modifications:
//
//---------------------------------------------------------------
#ifndef G4BUILDERINTERFACE_HH
#define G4BUILDERINTERFACE_HH

#include "globals.hh"


class G4PhysicsBuilderInterface {
  public:
    G4PhysicsBuilderInterface() = default;
    virtual ~G4PhysicsBuilderInterface() {} 
    virtual void Build() {
      G4Exception("G4PhysicsBuilderInterface::Build","PHYSBLD001",FatalException,
          "Called based class method. Should be implemented in"\
          " inherited class");;
    }
  virtual void RegisterMe(G4PhysicsBuilderInterface*) {
    G4Exception("G4PhysicsBuilderInterface::RegisterMe","PHYSBLD001",FatalException,
                "Called based class method. Should be implemented in"\
                " inherited class, or wrong type of parameter passed.");;
  }
    virtual void SetMinEnergy(G4double) {
      G4Exception("G4PhysicsBuilderInterface::SetMinEnergy","PHYSBLD001",
          FatalException,
          "Called based class method. Should be implemented in"\
          " inherited class");;
    }
    virtual void SetMaxEnergy(G4double) {
      G4Exception("G4PhysicsBuilderInterface::SetMaxEnergy","PHYSBLD001",
          FatalException,
          "Called based class method. Should be implemented in"\
          " inherited class");;
    }
};

#endif
