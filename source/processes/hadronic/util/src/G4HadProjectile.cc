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
#include "G4HadProjectile.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"

G4HadProjectile::G4HadProjectile() 
{
  theMat = 0;
  theDef = 0;
  theTime = 0.0;
  theBoundEnergy = 0.0;
}

G4HadProjectile::G4HadProjectile(const G4Track &aT) 
{
  Initialise(aT);
}

G4HadProjectile::G4HadProjectile(const G4DynamicParticle &aT) 
  : theMat(NULL),
    theOrgMom(aT.Get4Momentum()),
    theDef(aT.GetDefinition())
{
  G4LorentzRotation toZ;
  toZ.rotateZ(-theOrgMom.phi());
  toZ.rotateY(-theOrgMom.theta());
  theMom = toZ*theOrgMom;
  toLabFrame = toZ.inverse();
  theTime = 0.0;
  theBoundEnergy = 0.0;
}

G4HadProjectile::~G4HadProjectile()
{}

void G4HadProjectile::Initialise(const G4Track &aT)
{
  theMat = aT.GetMaterial();
  theOrgMom = aT.GetDynamicParticle()->Get4Momentum();
  theDef = aT.GetDefinition();

  G4LorentzRotation toZ;
  toZ.rotateZ(-theOrgMom.phi());
  toZ.rotateY(-theOrgMom.theta());
  theMom = toZ*theOrgMom;
  toLabFrame = toZ.inverse();

  //VI time of interaction starts from zero
  //   not global time of a track
  theTime = 0.0;
  theBoundEnergy = 0.0;
}

