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
// $Id: Tst10PrimaryGeneratorAction.hh,v 1.2 2001-07-11 10:09:48 gunter Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This class is a class derived from G4VUserPrimaryGeneratorAction
//      for constructing all particles and processes.
//
//	History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10PrimaryGeneratorAction_h
#define Tst10PrimaryGeneratorAction_h 1

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;

class Tst10PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    Tst10PrimaryGeneratorAction();
    ~Tst10PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
  private:
    G4ParticleGun* particleGun;
		G4ThreeVector GetRandomPolarization( G4ThreeVector direction );
		G4ThreeVector GetRandomDirection( );
		G4ThreeVector GetRandomPosition( );
};

#endif


