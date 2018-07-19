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
/// \file exoticphysics/dmparticle/include/TestParameters.hh
/// \brief Definition of the TestParameters class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HistoManager
//
// Description: Singleton class to make analysis and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              HistoManager::GetPointer() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

#ifndef TestParameters_h
#define TestParameters_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4ParticleDefinition.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TestParameters
{

public:
  // With description

  static TestParameters* GetPointer();

private:

  TestParameters();

public: // Without description

  ~TestParameters();

  void SetMaxEnergy(G4double value);

  G4double GetMaxEnergy() const;

  void SetNumberBins(G4int value);

  G4int GetNumberBins() const;

  void SetPositionZ(G4double val);

  G4double GetPositionZ() const;

  void SetBeamEnergy(G4double val);

  G4double GetBeamEnergy() const;

  void SetAlphaFactor(G4double val);

  G4double GetAlphaFactor() const;

  void SetBeamParticle(const G4ParticleDefinition*);

  const G4ParticleDefinition* GetBeamParticle() const;

private:

  static TestParameters* fManager;

  G4double fMaxEnergy;
  G4double fPositionZ;
  G4double fBeamEnergy;
  G4double fEpsilon;

  G4int fBinsE; 

  const G4ParticleDefinition* fParticle;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

