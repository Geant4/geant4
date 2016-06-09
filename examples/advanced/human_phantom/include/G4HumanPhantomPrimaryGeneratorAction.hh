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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#ifndef G4HumanPhantomPrimaryGeneratorAction_h
#define G4HumanPhantomPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include <vector>

//class G4HumanPhantomConstruction;
class G4HumanPhantomPrimaryGeneratorMessenger;
class G4ParticleGun;
class G4Event;

class G4HumanPhantomPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    G4HumanPhantomPrimaryGeneratorAction();
    ~G4HumanPhantomPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
    void GenerateBeamAlongZ();
    void GenerateBeamAlongX();
    void GenerateBeamAlongY();
    void GenerateIsotropicFlux();
    void SetBeam(G4String);

  private:
    G4ParticleGun* particleGun;

    G4double x0;
    G4double y0;
    G4double z0;
    
    std::vector<G4double> probability;
   
    G4HumanPhantomPrimaryGeneratorMessenger* messenger;
    
    G4String beamKind;
    G4double worldLength;
};

#endif


