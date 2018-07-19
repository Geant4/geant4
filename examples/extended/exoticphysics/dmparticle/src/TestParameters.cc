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
/// \file exoticphysics/dmparticle/src/TestParameters.cc
/// \brief Implementation of the TestParameters class
//
// $Id: TestParameters.cc 78707 2014-01-17 16:02:01Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   TestParameters
//
// Author:      V.Ivanchenko 01.09.2010
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TestParameters.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TestParameters* TestParameters::fManager = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TestParameters* TestParameters::GetPointer()
{
  if(!fManager) {
    fManager = new TestParameters();
  }
  return fManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TestParameters::TestParameters()
{
  fMaxEnergy   = 100.*GeV;
  fBinsE       = 100;
  fPositionZ   = 0.0;
  fBeamEnergy  = 0.0;
  fEpsilon     = 100.0;

  fParticle = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TestParameters::~TestParameters()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetMaxEnergy(G4double value)
{
  fMaxEnergy = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TestParameters::GetMaxEnergy() const
{
  return fMaxEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetNumberBins(G4int value)
{
  fBinsE = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TestParameters::GetNumberBins() const
{
  return fBinsE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetPositionZ(G4double val)
{
  fPositionZ = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TestParameters::GetPositionZ() const
{
  return fPositionZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetBeamEnergy(G4double val)
{
  fBeamEnergy = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TestParameters::GetBeamEnergy() const
{
  return fBeamEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetAlphaFactor(G4double val)
{
  fEpsilon = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TestParameters::GetAlphaFactor() const
{
  return fEpsilon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TestParameters::SetBeamParticle(const G4ParticleDefinition* ptr)
{
  fParticle = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4ParticleDefinition* TestParameters::GetBeamParticle() const
{
  return fParticle;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
