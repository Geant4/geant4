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

#ifndef G4BertiniEvaporationChannel_h
#define G4BertiniEvaporationChannel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4NucleiProperties.hh" 

class G4BertiniEvaporationChannel
{
public:
  G4BertiniEvaporationChannel();
  virtual ~G4BertiniEvaporationChannel(); 

  void setVerboseLevel( G4int verbose );
  void setNucleusA( G4int inputA );
  void setNucleusZ( G4int inputZ );
  G4int getParticleA();
  G4int getParticleZ();
  void setExcitationEnergy( G4double inputE );
  void setQ( G4double inputQ );
  void setPairingCorrection( G4int isCorrection );
  G4double getLevelDensityParameter();
  G4String getName();
  
  virtual G4double getProbability() ;
  virtual void setProbability( G4double newProb ) ;
  virtual void calculateProbability() = 0 ;
  virtual G4double qmFactor ();
  virtual G4double getQ();
  virtual G4double getCoulomb();
  virtual G4double getThresh();

  virtual G4DynamicParticle* emit() = 0;
  
protected:  
  G4String name;
  G4int verboseLevel;
  G4int nucleusA;
  G4int nucleusZ;
  G4int particleA;
  G4int particleZ;
  G4double exmass;
  G4double emissionProbability;
  G4double rho;
  G4double correction;
  G4double excitationEnergy;
  //  G4int massInNeutronMasses;
  G4int spin;
  G4double Q( G4double a, G4double z );
  G4double pairingEnergyProtons( G4int A );
  G4double pairingEnergyNeutrons( G4int N );
  G4double cameron( G4double a, G4double z );
  G4double cameronShellCorrectionP( G4int Z );
  G4double cameronShellCorrectionN( G4int N );
  void isotropicCosines( G4double&,
			 G4double&,
			 G4double& );
};


#endif
