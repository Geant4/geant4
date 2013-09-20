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
// $Id: G4RTPrimaryGeneratorAction.hh 66241 2012-12-13 18:34:42Z gunter $
//

#ifndef G4RTPrimaryGeneratorAction_h
#define G4RTPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
class G4Event;
class G4ParticleDefinition;
class G4VPhysicalVolume;
#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4RTPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    G4RTPrimaryGeneratorAction();
    virtual ~G4RTPrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event* anEvent);
    void SetUp();

  private:
    G4ParticleDefinition* particle_definition;
    G4double              particle_energy;
    G4double              particle_time;
    G4ThreeVector         particle_polarization;

    G4VPhysicalVolume*    pWorld;
    EInside               whereisit;

    G4int nColumn;
    G4int nRow;

    G4ThreeVector eyePosition;
    G4ThreeVector eyeDirection;
    G4ThreeVector up;
    G4double headAngle;
    G4double viewSpan;
    G4double stepAngle;
    G4double viewSpanX;
    G4double viewSpanY;

    G4bool distortionOn;
};

#endif


