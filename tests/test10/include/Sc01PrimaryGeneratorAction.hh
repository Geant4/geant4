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
// $Id: Sc01PrimaryGeneratorAction.hh,v 1.1 2004-01-27 11:19:22 grichine Exp $
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


