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
/// \file runAndEvent/RE04/include/RE04PhysicsList.hh
/// \brief Definition of the RE04PhysicsList class
//
// $Id: $
//
#ifndef RE04PhysicsList_h
#define RE04PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

//
/// Physics list class
///
/// - constructor
///     registers EM, synchrotron radiation, GN, decays, hadron
///     elastic scattering, hadron, stopping and ion physics
///
/// - void ConstructProcess()
///     invokes AddParallelWorldProcess() and
///     constructs processes with a G4VPhysicsConstructor vector
///
/// - AddParallelWorldProcess()
///     adds a parallel world process, "paraWorldProc"
///
/// - void SetCuts()
///     invokes default SetCuts methods with SetCutsWithDefault()
//
class RE04PhysicsList: public G4VModularPhysicsList
{
  public:
    RE04PhysicsList(G4String& parWorldName);
    virtual ~RE04PhysicsList();

  public:
    // Construct particle and physics
    virtual void ConstructProcess();
    virtual void SetCuts();

  private:
    void AddParallelWorldProcess();

  private:
    G4String fpWorldName;

};

#endif



