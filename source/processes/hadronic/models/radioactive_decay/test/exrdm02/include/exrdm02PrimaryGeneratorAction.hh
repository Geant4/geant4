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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: exrdm02PrimaryGeneratorAction.hh,v 1.1 2003-10-08 16:31:50 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef exrdm02PrimaryGeneratorAction_h
#define exrdm02PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4IonTable.hh"
#include "RadioactiveDecayGun.hh"

//class G4ParticleGun;

class G4Event;

class exrdm02PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    exrdm02PrimaryGeneratorAction();
    ~exrdm02PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  protected:

     RadioactiveDecayGun *theParticleGun;

};

#endif


