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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4EmElementXS
//
// Author:        V. Ivanchenko 
// 
// Creation date: 22 August 2024
//
// -------------------------------------------------------------------
//

#include "G4EmElementXS.hh"
#include "G4ElementData.hh"
#include "G4ElementDataRegistry.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4AutoLock.hh"
#include "G4SystemOfUnits.hh"

namespace
{
  G4Mutex elementXSMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4EmElementXS::G4EmElementXS(G4int zmin, G4int zmax, const G4String& name,
			     const G4String& subname)
  : Zmin(zmin - 1), Zmax(zmax), fSubName(subname)
{
  fParameters = G4EmParameters::Instance();
  auto reg = G4ElementDataRegistry::Instance();
  fData = reg->GetElementDataByName(name);
  if (nullptr == fData) {
    fData = new G4ElementData(Zmax - Zmin);
    fData->SetName(name);
    reg->RegisterMe(fData);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4EmElementXS::Retrieve(G4int ZZ) const
{
  G4int Z = std::min(ZZ, Zmax);
  auto v = fData->GetElementData(Z);
  if (nullptr == v) {
    G4AutoLock l(&elementXSMutex);
    v = fData->GetElementData(Z);
    if (nullptr == v) {
      v = new G4PhysicsFreeVector(false);
      std::ostringstream ost;
      ost << fParameters->GetDirLEDATA() << fSubName << Z << ".dat";
      std::ifstream fin(ost.str().c_str());
      if (!fin.is_open()) {
	G4ExceptionDescription ed;
	ed << "G4EmElementXS: data file <" << ost.str().c_str() << "> for Z=" << Z
	   << " is not opened!" << G4endl;
	G4Exception("G4EmElementXS::Retrieve()", "em0003", FatalException, ed,
		    "G4LEDATA version should be checked");
      } else {
	v->Retrieve(fin, true);
	v->ScaleVector(CLHEP::MeV, CLHEP::barn);
      }
      fData->InitialiseForElement(Z, v);
      l.unlock();
    }
  }
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4EmElementXS::GetXS(G4int Z, G4double ekin) const
{
  auto v = Retrieve(Z);
  return (nullptr != v) ? v->Value(ekin) : 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

