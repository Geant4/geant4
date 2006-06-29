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
// $Id: Sc01PrimaryGeneratorAction.hh,v 1.2 2006-06-29 18:53:47 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPrimaryGeneratorAction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Sc01PrimaryGeneratorAction_h
#define Sc01PrimaryGeneratorAction_h 1

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;
class Sc01DetectorConstruction;


class Sc01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Sc01PrimaryGeneratorAction(Sc01DetectorConstruction*);
    ~Sc01PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
  private:
    G4ParticleGun* particleGun;
		G4ThreeVector GetRandomPolarization( G4ThreeVector direction );
		G4ThreeVector GetRandomDirection( );
		G4ThreeVector GetRandomPosition( );
  Sc01DetectorConstruction* fDetector;
};

#endif


