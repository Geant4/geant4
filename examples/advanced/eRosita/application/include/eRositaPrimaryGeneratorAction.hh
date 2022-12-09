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

#ifndef eRositaPrimaryGeneratorAction_h
#define eRositaPrimaryGeneratorAction_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class G4ParticleGun;

class eRositaPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    explicit eRositaPrimaryGeneratorAction();

    ~eRositaPrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;

private:
    G4ParticleGun* particleGun;

    G4double positionX; // x position of a vertex

    G4double positionY; // y position of a vertex

    G4double positionZ; // z position of a vertex

    static constexpr auto INITIAL_MOMENTUM_DIRECTION_X{0.0};
    G4double momentumDirectionX{INITIAL_MOMENTUM_DIRECTION_X}; // x component of initial momentum vector

    static constexpr auto INITIAL_MOMENTUM_DIRECTION_Y{-0.5};
    G4double momentumDirectionY{INITIAL_MOMENTUM_DIRECTION_Y}; // y component of initial momentum vector

    static constexpr auto INITIAL_MOMENTUM_DIRECTION_Z{-1.0};
    G4double momentumDirectionZ{INITIAL_MOMENTUM_DIRECTION_Z}; // z component of initial momentum vector
};
#endif
