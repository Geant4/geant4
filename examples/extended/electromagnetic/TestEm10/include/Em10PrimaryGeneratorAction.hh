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
/// \file electromagnetic/TestEm10/include/Em10PrimaryGeneratorAction.hh
/// \brief Definition of the Em10PrimaryGeneratorAction class
//
//
// $Id: Em10PrimaryGeneratorAction.hh 73033 2013-08-15 09:24:45Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10PrimaryGeneratorAction_h
#define Em10PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class Em10DetectorConstruction;
class Em10PrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Em10PrimaryGeneratorAction(Em10DetectorConstruction*);
   ~Em10PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    void SetRndmFlag(G4String val) { rndmFlag = val;}
    void Setxvertex(G4double x);
    void Setyvertex(G4double y);
    void Setzvertex(G4double z);

    static G4String GetPrimaryName();

  private:
    G4ParticleGun*                particleGun; //pointer a to G4 service class
//    Em10DetectorConstruction*      Em10Detector; //pointer to the geometry

    Em10PrimaryGeneratorMessenger* gunMessenger; //messenger of this class
    G4String                      rndmFlag;  //flag for a random impact point

    static G4String thePrimaryParticleName;
    G4double xvertex,yvertex,zvertex;
    G4bool vertexdefined;

};

#endif
