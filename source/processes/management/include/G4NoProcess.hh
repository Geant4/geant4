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

#ifndef G4NoProcess_h
#define G4NoProcess_h 1

#include "G4VProcess.hh"

class G4NoProcess : public G4VProcess {
  public:

    G4NoProcess() : G4VProcess( "NoProcess", fGeneral ) {};

    virtual ~G4NoProcess() {};

    // This process should not be set to any particle
    virtual G4bool IsApplicable(const G4ParticleDefinition&) override { return false; }

    //  no operations in any GPIL or DoIt
    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double,
                             G4ForceCondition*
                            ) override { return -1.0; };

    virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& ,
                             const G4Step&
                            ) override {return nullptr;};

    virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4ForceCondition*
                            ) override { return -1.0; };

    virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& ,
                             const G4Step&
                            ) override {return nullptr;};

    virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
                             G4double  ,
                             G4double  ,
                             G4double& ,
                             G4GPILSelection*
                            ) override { return -1.0; };

    virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& ,
                             const G4Step&
                            ) override {return nullptr;};

  private:

    // hide copy constr and assignment operator as private
    G4NoProcess(G4NoProcess&);
    G4NoProcess& operator=(const G4NoProcess& right);

};

#endif
