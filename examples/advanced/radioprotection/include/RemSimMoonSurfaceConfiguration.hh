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
// $Id: RemSimMoonSurfaceConfiguration.hh,v 1.1 2004-05-17 10:34:56 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RemSimMoonSurfaceConfiguration_h
#define RemSimMoonSurfaceConfiguration_h 1

#include "RemSimVPrimaryGeneratorFactory.hh"
#include "globals.hh"
#include "G4DataVector.hh"
class G4ParticleGun;
class G4Event;
class RemSimRunAction;

class RemSimMoonSurfaceConfiguration : public RemSimVPrimaryGeneratorFactory
{
public:
  RemSimMoonSurfaceConfiguration();
  ~RemSimMoonSurfaceConfiguration();

public:
  void GeneratePrimaries(G4Event* anEvent);
   G4double  GetInitialEnergy();

private:
  G4ParticleGun* particleGun;
  G4String randomDirection;
  G4String spectrum;
  G4double energy;
  RemSimRunAction* run;
};
#endif


