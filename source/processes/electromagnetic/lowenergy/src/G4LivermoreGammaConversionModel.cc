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
// Author: Sebastien Incerti
//         22 January 2012
//         on base of G4LivermoreGammaConversionModel (original version)
//         and G4LivermoreRayleighModel (MT version)
//
// Modifications: Zhuxin Li@CENBG
//                11 March 2020
//                derives from G4PairProductionRelModel
// -------------------------------------------------------------------

#include "G4LivermoreGammaConversionModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4AutoLock.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"

namespace { G4Mutex LivermoreGammaConversionModelMutex = G4MUTEX_INITIALIZER; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreGammaConversionModel::lowEnergyLimit = 2.*CLHEP::electron_mass_c2;
G4PhysicsFreeVector* G4LivermoreGammaConversionModel::data[] = {nullptr};

G4LivermoreGammaConversionModel::G4LivermoreGammaConversionModel
(const G4ParticleDefinition* p, const G4String& nam)
: G4PairProductionRelModel(p,nam),fParticleChange(nullptr),maxZ(100)
{
  verboseLevel = 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods
  if(verboseLevel > 0)
  {
    G4cout << "G4LivermoreGammaConversionModel is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreGammaConversionModel::~G4LivermoreGammaConversionModel()
{
  if(IsMaster()) {
    for(G4int i = 0; i <= maxZ; ++i)
    {
      if(data[i])
      {
        delete data[i];
        data[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::Initialise(
                                const G4ParticleDefinition* particle,
				const G4DataVector& cuts)
{ G4PairProductionRelModel::Initialise(particle, cuts);
  if (verboseLevel > 1)
  {
    G4cout << "Calling Initialise() of G4LivermoreGammaConversionModel."
	   << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV isMater: " << IsMaster()
	   << G4endl;
  }

  if(fParticleChange == nullptr) {
    fParticleChange = GetParticleChangeForGamma();
  }

  if(IsMaster())
  {
    // Initialise element selector
    InitialiseElementSelectors(particle, cuts);

    // Access to elements
    const char* path = G4FindDataDir("G4LEDATA");
    const G4ElementTable* elemTable = G4Element::GetElementTable();
    size_t numElems                 = (*elemTable).size();
    for(size_t ie = 0; ie < numElems; ++ie)
    {
      const G4Element* elem = (*elemTable)[ie];
      const G4int Z         = std::min(maxZ, elem->GetZasInt());
      if(data[Z] == nullptr)
      {
        ReadData(Z, path);
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversionModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1)
  {
    G4cout << "Calling ReadData() of G4LivermoreGammaConversionModel"
	   << G4endl;
  }

  if(data[Z]!= nullptr) { return; }

  const char* datadir = path;

  if(datadir == nullptr)
  {
    datadir = G4FindDataDir("G4LEDATA");
    if(datadir == nullptr)
    {
      G4Exception("G4LivermoreGammaConversionModel::ReadData()",
		  "em0006",FatalException,
		  "Environment variable G4LEDATA not defined");
      return;
    }
  }

  std::ostringstream ost;
  if(G4EmParameters::Instance()->LivermoreDataDir() == "livermore"){
    data[Z] = new G4PhysicsFreeVector(true);
    ost << datadir << "/livermore/pair/pp-cs-" << Z <<".dat";
  }else{
    data[Z] = new G4PhysicsFreeVector();
    ost << datadir << "/epics2017/pair/pp-cs-" << Z <<".dat";
  }

  std::ifstream fin(ost.str().c_str());

  if( !fin.is_open())
  {
    G4ExceptionDescription ed;
    ed << "G4LivermoreGammaConversionModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermoreGammaConversionModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW8.0 or later.");
    return;
  }
  else
  {

    if(verboseLevel > 1) { G4cout << "File " << ost.str()
	     << " is opened by G4LivermoreGammaConversionModel" << G4endl;}

    data[Z]->Retrieve(fin, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreGammaConversionModel::ComputeCrossSectionPerAtom(
           const G4ParticleDefinition* particle,
	   G4double GammaEnergy, G4double Z, G4double, G4double, G4double)
{
  if (verboseLevel > 1)
  {
    G4cout << "G4LivermoreGammaConversionModel::ComputeCrossSectionPerAtom() Z= "
	   << Z << G4endl;
  }

  if (GammaEnergy < lowEnergyLimit) { return 0.0; }

  G4double xs = 0.0;

  G4int intZ = std::max(1, std::min(G4lrint(Z), maxZ));

  G4PhysicsFreeVector* pv = data[intZ];

  // if element was not initialised
  // do initialisation safely for MT mode
  if(pv == nullptr)
  {
    InitialiseForElement(particle, intZ);
    pv = data[intZ];
    if(pv == nullptr) { return xs; }
  }
  // x-section is taken from the table
  xs = pv->Value(GammaEnergy);

  if(verboseLevel > 0)
    {
      G4cout  <<  "*** Gamma conversion xs for Z=" << Z << " at energy E(MeV)="
	      << GammaEnergy/MeV <<  "  cs=" << xs/millibarn << " mb" << G4endl;
    }
  return xs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LivermoreGammaConversionModel::InitialiseForElement(
				      const G4ParticleDefinition*,
				      G4int Z)
{
  G4AutoLock l(&LivermoreGammaConversionModelMutex);
  if(data[Z] == nullptr) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
