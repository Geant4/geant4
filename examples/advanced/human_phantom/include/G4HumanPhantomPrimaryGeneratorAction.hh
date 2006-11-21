//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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


