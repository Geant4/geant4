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
//
/// \file field/field04/include/F04PhysicsList.hh
/// \brief Definition of the F04PhysicsList class
//

#ifndef F04PhysicsList_h
#define F04PhysicsList_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class F04PhysicsListMessenger;

class F04StepMax;

class F04PhysicsList: public G4VModularPhysicsList
{
public:

    F04PhysicsList(G4String);
    virtual ~F04PhysicsList();

    void SetStepMax(G4double);
    F04StepMax* GetStepMaxProcess();
    void AddStepMax();
/*
    /// Remove specific physics from physics list.
    void RemoveFromPhysicsList(const G4String&);

    /// Make sure that the physics list is empty.
    void ClearPhysics();
*/
    virtual void ConstructParticle();
    virtual void ConstructProcess();

private:

    G4double fMaxChargedStep;
    static G4ThreadLocal F04StepMax* fStepMaxProcess;

    F04PhysicsListMessenger* fMessenger;

};

#endif
