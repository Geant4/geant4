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
//
// $Id: MyPrimaryGeneratorAction.hh,v 1.5 2005-02-04 16:28:46 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyPrimaryGeneratorAction_h
#define MyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

// #define GPS // For General Particle Source.

#ifdef GPS
class G4GeneralParticleSource;
#else
class G4ParticleGun;
#endif
class G4Event;

class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MyPrimaryGeneratorAction();
    ~MyPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);
#ifdef GPS
    G4GeneralParticleSource* GetGeneralParticleSource();
#else
    G4ParticleGun* GetParticleGun();
#endif

  private:
#ifdef GPS
    G4GeneralParticleSource* particleGun;
#else
    G4ParticleGun* particleGun;
#endif
};

#endif


