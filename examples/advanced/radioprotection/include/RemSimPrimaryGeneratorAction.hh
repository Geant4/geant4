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
// $Id: RemSimPrimaryGeneratorAction.hh,v 1.14 2006-06-29 16:23:07 gunter Exp $// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli, guatelli@ge.infn.it
//
#ifndef RemSimPrimaryGeneratorAction_h
#define RemSimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4DataVector;
class RemSimPrimaryGeneratorMessenger;

class RemSimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  RemSimPrimaryGeneratorAction();
  ~RemSimPrimaryGeneratorAction();

public:
  G4double GetInitialEnergy();
  void GeneratePrimaries(G4Event* anEvent);
  void ReadProbability(G4String);

private:  
  G4ParticleGun* particleGun;
  RemSimPrimaryGeneratorMessenger* messenger;
  G4bool readFile;

  G4double* cumulate;
  G4double* energy;

 //stores the probability given in input
  G4DataVector* data;
 //stores the energy data 
  G4DataVector* energies;

  G4int size;
};
#endif


