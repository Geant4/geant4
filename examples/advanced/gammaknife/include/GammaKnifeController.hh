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


#ifndef GAMMAKNIFECONTROLLER_HH
#define GAMMAKNIFECONTROLLER_HH

#define GAMMAKNIFE_SOURCES 201

#include <vector>
#include <string>
#include <map>

#include "G4ScoringManager.hh"
#include "G4VScoringMesh.hh"

#include "GammaKnifeDetectorConstruction.hh"
#include "GammaKnifeMessenger.hh"

class GammaKnifeMessenger;

class GammaKnifeController
{
public:
    GammaKnifeController( GammaKnifeDetectorConstruction* );

    ~GammaKnifeController();

    // Run G4int events for all 201 positions
    void BeamOn( G4int );

    // Load positions from an external file
    void ReadFile( std::string fileName );

private:
    void StoreHits();

    void PrepareHitsAccumulation();

    void AccumulateAllHits();

    void RotateForward( G4int );

    void RotateBack( G4int );

    std::vector<G4double> phiAngles;

    std::vector<G4double> thetaAngles;

    GammaKnifeDetectorConstruction* detector;

    GammaKnifeMessenger* messenger;

    using MeshScoreMap = G4VScoringMesh::MeshScoreMap;
    MeshScoreMap* scoreMaps;
};

#endif // GAMMAKNIFECONTROLLER_HH
