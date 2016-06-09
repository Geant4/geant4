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
// $Id$
//
//
//
//    ********************************************
//    *                                          *
//    *      BrachyPrimaryGeneratorActionI.hh     *
//    *                                          *
//    ********************************************
//Primary particles of the Iodium source. They are delivered from a random
//point of the radionuclide  with random direction. The energy spectrum
// of the gamma delivered is activated

#ifndef BrachyPrimaryGeneratorActionI_h
#define BrachyPrimaryGeneratorActionI_h 1

#include "globals.hh"
#include <vector>
#include "G4VUserPrimaryGeneratorAction.hh"
#include "BrachyPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Run;
class G4Event;
class BrachyAnalysisManager;
class BrachyPrimaryGeneratorAction;

class BrachyPrimaryGeneratorActionI : public  G4VUserPrimaryGeneratorAction
{
 public:
      BrachyPrimaryGeneratorActionI();
      ~BrachyPrimaryGeneratorActionI();

 public:
      void GeneratePrimaries(G4Event* anEvent);
      G4double GetEnergy();

 private:
      G4ParticleGun* particleGun;
      G4double primaryParticleEnergy;
      std::vector<G4double> energySpectrum;
};

#endif

