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
// Implementation of the HETC88 code into Geant4.
// Evaporation and De-excitation parts
// T. Lampen, Helsinki Institute of Physics, May-2000

#ifndef G4BEChargedChannel_h
#define G4BEChargedChannel_h 1

#include "globals.hh"
#include "G4BertiniEvaporationChannel.hh"
#include "Randomize.hh"

class G4BEChargedChannel : public G4BertiniEvaporationChannel
{
public:
  G4BEChargedChannel();
  virtual ~G4BEChargedChannel(); 
  
  virtual void calculateProbability();
  virtual G4DynamicParticle * emit() = 0;
  virtual G4double coulombFactor() = 0;
  G4double coulombFactorForProton();
  G4double qmFactorForProton();
  G4double qmFactorForAlpha();
  G4double sampleKineticEnergy();
  
protected:  
  G4double A;
  G4double spin;
};


#endif
