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
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelPrimaryGeneratorAction  ------
//           by  G.Santin, F.Longo & R.Giannitrapani (13 nov 2000)
//
// ************************************************************

#include "G4RunManager.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPrimaryGeneratorMessenger.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::GammaRayTelPrimaryGeneratorAction() {
    detector = static_cast<const GammaRayTelDetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    // create a messenger for this class

    gunMessenger = new GammaRayTelPrimaryGeneratorMessenger(this);

    constexpr auto NUMBER_OF_PARTICLES{1};
    particleGun = new G4ParticleGun(NUMBER_OF_PARTICLES);

    // default particle kinematic

    auto *particleTable = G4ParticleTable::GetParticleTable();
    auto *particle = particleTable->FindParticle("e-");
    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

    constexpr auto PARTICLE_ENERGY{30. * MeV};
    particleGun->SetParticleEnergy(PARTICLE_ENERGY);

    auto position = 0.5 * (detector->GetWorldSizeZ());
    particleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, position));
    particleSource = new G4GeneralParticleSource();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelPrimaryGeneratorAction::~GammaRayTelPrimaryGeneratorAction() {
    delete particleGun;
    delete particleSource;
    delete gunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelPrimaryGeneratorAction::GeneratePrimaries(G4Event *event) {
    if (sourceGun) {
        G4cout << "Using G4ParticleGun... " << G4endl;

        // this function is called at the begining of event
        //
        G4double x0 = 0. * cm;
        G4double y0 = 0. * cm;
        G4double z0 = 0.5 * (detector->GetWorldSizeZ());

        G4ThreeVector pos0;
        auto vertex0 = G4ThreeVector(x0, y0, z0);
        auto momentumDirection0 = G4ThreeVector(0., 0., -1.);

        G4double theta;
        G4double phi;
        G4double y = 0.;
        G4double f = 0.;
        G4double theta0 = 0.;
        G4double phi0 = 0.;

        switch (sourceType) {
        case 0:
            particleGun->SetParticlePosition(vertex0);
            particleGun->SetParticleMomentumDirection(momentumDirection0);
            break;
        case 1:
            // GS: Generate random position on the 4PIsphere to create a uniform distribution
            // GS: on the sphere
            phi = G4UniformRand() * twopi;
            do {
                y = G4UniformRand() * 1.0;
                theta = G4UniformRand() * pi;
                f = std::sin(theta);
            } while (y > f);
            vertex0 = G4ThreeVector(1., 0., 0.);
            vertex0.setMag(vertexRadius);
            vertex0.setTheta(theta);
            vertex0.setPhi(phi);
            particleGun->SetParticlePosition(vertex0);

            momentumDirection0 = G4ThreeVector(1., 0., 0.);

            do {
                phi = G4UniformRand() * twopi;
                do {
                    y = G4UniformRand() * 1.0;
                    theta = G4UniformRand() * pi;
                    f = std::sin(theta);
                } while (y > f);
                momentumDirection0.setPhi(phi);
                momentumDirection0.setTheta(theta);
            } while (vertex0.dot(momentumDirection0) >= -0.7 * vertex0.mag());

            particleGun->SetParticleMomentumDirection((G4ParticleMomentum) momentumDirection0);

            break;
        case 2:
            // GS: Generate random position on the upper semi-sphere z > 0 to create a uniform distribution
            // GS: on a plane
            phi = G4UniformRand() * twopi;

            do {
                y = G4UniformRand() * 1.0;
                theta = G4UniformRand() * halfpi;
                f = std::sin(theta) * std::cos(theta);
            } while (y > f);

            vertex0 = G4ThreeVector(1., 0., 0.);

            auto xy = detector->GetWorldSizeXY();
            auto z = detector->GetWorldSizeZ();

            if (vertexRadius > xy * 0.5) {
                G4cout << "vertexRadius too big " << G4endl;
                G4cout << "vertexRadius set to " << xy * 0.45 << G4endl;
                vertexRadius = xy * 0.45;
            }

            if (vertexRadius > z * 0.5) {
                G4cout << "vertexRadius too high " << G4endl;
                G4cout << "vertexRadius set to " << z * 0.45 << G4endl;
                vertexRadius = z * 0.45;
            }

            vertex0.setMag(vertexRadius);
            vertex0.setTheta(theta);
            vertex0.setPhi(phi);

            // GS: Get the user defined direction for the primaries and
            // GS: Rotate the random position according to the user defined direction for the particle

            momentumDirection0 = particleGun->GetParticleMomentumDirection();
            if (momentumDirection0.mag() > 0.001) {
                theta0 = momentumDirection0.theta();
                phi0 = momentumDirection0.phi();
            }

            if (theta0 != 0.) {
                G4ThreeVector rotationAxis(1., 0., 0.);
                rotationAxis.setPhi(phi0 + halfpi);
                vertex0.rotate(theta0 + pi, rotationAxis);
            }
            particleGun->SetParticlePosition(vertex0);
            break;
        }

        constexpr auto INITIAL_PARTICLE_ENERGY{100 * MeV};
        G4double particleEnergy = INITIAL_PARTICLE_ENERGY;

        switch (spectrumType) {
        case 0: // Uniform energy (1 GeV - 10 GeV)
            y = G4UniformRand();
            particleEnergy = y * 9.0 * GeV + 1.0 * GeV;
            G4cout << "Particle energy: " << particleEnergy / GeV << " LIN" << G4endl;
            break;
        case 1: // Logarithmic energy
            y = G4UniformRand();
            particleEnergy = std::pow(10, y) * GeV;
            G4cout << "Particle energy: " << particleEnergy / GeV << " LOG" << G4endl;
            break;
        case 2: // Power law (-4)
            do {
                y = G4UniformRand() * 100000.0;
                particleEnergy = G4UniformRand() * 10. * GeV;
                f = std::pow(particleEnergy * (1 / GeV), -4.);
            } while (y > f);
            // particleGun->SetParticleEnergy(particleEnergy);
            break;
        case 3: // Monochromatic
            particleEnergy = particleGun->GetParticleEnergy();
            // 100 MeV;
            G4cout << "Particle energy: " << particleEnergy << " MONOCHROMATIC" << G4endl;
            break;
        }
        particleGun->SetParticleEnergy(particleEnergy);
        G4cout << "Particle: " << particleGun->GetParticleDefinition()->GetParticleName() << G4endl;
        particleGun->GeneratePrimaryVertex(event);
    } else {
        particleSource->GeneratePrimaryVertex(event);
    }
}
