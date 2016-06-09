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
//    **********************************
//    *                                *
//    *    RemSimBasicGenerator.hh     *
//    *                                *
//    **********************************
//
// $Id: RemSimBasicGenerator.hh,v 1.7 2004/05/27 10:33:11 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 

#ifndef RemSimBasicGenerator_h
#define RemSimBasicGenerator_h 1

#include "RemSimVPrimaryGeneratorFactory.hh"
#include "globals.hh"
#include "G4DataVector.hh"
class G4ParticleGun;
class G4Event;
class RemSimRunAction;

class RemSimBasicGenerator : public RemSimVPrimaryGeneratorFactory
{
public:
  RemSimBasicGenerator();
  ~RemSimBasicGenerator();

public:
  void GeneratePrimaries(G4Event* anEvent); 
  G4double  GetInitialEnergy();
  void SetMoon(G4bool){;};

private:
  G4ParticleGun* particleGun;
  G4String randomDirection;
  G4String spectrum;
  G4double energy;
  RemSimRunAction* run;
};
#endif


