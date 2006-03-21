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
// $Id: G4CoulombScattering.cc,v 1.1 2006-03-21 15:40:54 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4CoulombScattering
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 22.08.2004
//
// Modifications:
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4CoulombScattering.hh"
#include "G4CoulombScatteringModel.hh"
#include "G4eCoulombScatteringModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4CoulombScattering::G4CoulombScattering(const G4String& name)
  : G4VEmProcess(name),thetaMin(0.5),thetaMax(pi),q2Max(DBL_MAX),isInitialised(false)
{
  SetLambdaBinning(80);
  SetMinKinEnergy(1.0*keV);
  SetMaxKinEnergy(100.0*GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CoulombScattering::~G4CoulombScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  if(!isInitialised) {
    isInitialised = true;
    //    SetVerboseLevel(3);
    SetBuildTableFlag(true);
    SetStartFromNullFlag(false);
    SetLambdaFactor(0.8);
    SetSecondaryParticle(0);
    G4double emin = MinKinEnergy();
    G4double emax = MaxKinEnergy();
    G4CoulombScatteringModel* model = 
      new G4CoulombScatteringModel(thetaMin,thetaMax,true,q2Max);
    model->SetLowEnergyLimit(emin);
    model->SetHighEnergyLimit(emax);
    AddEmModel(1, model);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4CoulombScattering::PrintInfo()
{
  G4cout << "      Coulomb scattering with ThetaMin(degree)= " << thetaMin/degree
	 << "; ThetaMax(degree)= " << thetaMin/degree
         << "; q2Max(GeV^2)= " << q2Max/(GeV*GeV)
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
