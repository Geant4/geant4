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
// Author: Zhuxin Li@CENBG
//         11 March 2020
//         on the base of G4LivermoreGammaConversionModel 
//         derives from G4BetheHeitler5DModel              
// -------------------------------------------------------------------

#include "G4LivermoreGammaConversion5DModel.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4EmParameters.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Exp.hh"
#include "G4AutoLock.hh"

namespace { G4Mutex LivermoreGammaConversion5DModelMutex = G4MUTEX_INITIALIZER; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4int G4LivermoreGammaConversion5DModel::maxZ;
G4double G4LivermoreGammaConversion5DModel::lowEnergyLimit = 2.*CLHEP::electron_mass_c2;
G4PhysicsFreeVector* G4LivermoreGammaConversion5DModel::data[] = {nullptr};

G4LivermoreGammaConversion5DModel::G4LivermoreGammaConversion5DModel(const G4ParticleDefinition* p, 
								     const G4String& nam)
  : G4BetheHeitler5DModel(p, nam), fParticleChange(nullptr)
{
  verboseLevel = 0;
  // Verbosity scale for debugging purposes:
  // 0 = nothing 
  // 1 = calculation of cross sections, file openings...
  // 2 = entering in methods
  if(verboseLevel > 0) 
  {
    G4cout << "G4LivermoreGammaConversion5DModel is constructed " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreGammaConversion5DModel::~G4LivermoreGammaConversion5DModel()
{
  if(IsMaster()) {
    for(G4int i=0; i<maxZ; ++i) {
      if(data[i]) { 
	delete data[i];
	data[i] = nullptr;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermoreGammaConversion5DModel::Initialise( const G4ParticleDefinition* particle,
				                                     const G4DataVector& cuts)
{
 G4BetheHeitler5DModel::Initialise(particle, cuts);
   
 if (verboseLevel > 1) 
   {
     G4cout << "Calling Initialise() of G4LivermoreGammaConversion5DModel." 
	    << G4endl
	    << "Energy range: "
	    << LowEnergyLimit() / MeV << " MeV - "
	    << HighEnergyLimit() / GeV << " GeV isMater: " << IsMaster() 
	    << G4endl;
   }
 
 if(!fParticleChange) {
   fParticleChange = GetParticleChangeForGamma();
  }
 
 if(IsMaster()) 
   {
     // Initialise element selector
     InitialiseElementSelectors(particle, cuts);
     // Access to elements
     const char* path = G4FindDataDir("G4LEDATA");
     G4ProductionCutsTable* theCoupleTable =
       G4ProductionCutsTable::GetProductionCutsTable();
     G4int numOfCouples = G4int(theCoupleTable->GetTableSize());
     for(G4int i=0; i<numOfCouples; ++i) 
       {
	 const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
	 SetCurrentCouple(couple);
	 const G4Material* mat = couple->GetMaterial();
	 const G4ElementVector* theElementVector = mat->GetElementVector();
	 std::size_t nelm = mat->GetNumberOfElements();
	 for (std::size_t j=0; j<nelm; ++j) 
	   {
	     G4int Z = std::max(1, std::min((*theElementVector)[j]->GetZasInt(), maxZ));
	     if(!data[Z]) { ReadData(Z, path); }
	   }
       }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreGammaConversion5DModel::ReadData(size_t Z, const char* path)
{
  if (verboseLevel > 1) 
    {
      G4cout << "Calling ReadData() of G4LivermoreGammaConversion5DModel" 
	   << G4endl;
    }
  
  if(data[Z]) { return; }
  const char* datadir = path;
  if(!datadir) 
    {
    datadir = G4FindDataDir("G4LEDATA");
    if(!datadir) 
    {
      G4Exception("G4LivermoreGammaConversion5DModel::ReadData()",
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
    ed << "G4LivermoreGammaConversion5DModel data file <" << ost.str().c_str()
       << "> is not opened!" << G4endl;
    G4Exception("G4LivermoreGammaConversion5DModel::ReadData()",
		"em0003",FatalException,
		ed,"G4LEDATA version should be G4EMLOW8.0 or later.");
    return;
  }   
  else 
    {
      if(verboseLevel > 1) { G4cout << "File " << ost.str() 
				    << " is opened by G4LivermoreGammaConversion5DModel" << G4endl;} 
      data[Z]->Retrieve(fin, true);
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreGammaConversion5DModel::ComputeCrossSectionPerAtom(
   const G4ParticleDefinition* particle, G4double GammaEnergy, G4double Z, 
   G4double, G4double, G4double)
{
  if (verboseLevel > 1) 
  {
    G4cout << "G4LivermoreGammaConversion5DModel::ComputeCrossSectionPerAtom() Z= " 
	   << Z << G4endl;
  }
  G4double xs = 0.0;
  if (GammaEnergy < lowEnergyLimit) { return xs; } 
  
  G4int intZ = std::max(1, std::min(G4lrint(Z), maxZ));
  G4PhysicsFreeVector* pv = data[intZ];
  // if element was not initialised
  // do initialisation safely for MT mode
  if(!pv) 
  {
    InitialiseForElement(particle, intZ);
    pv = data[intZ];
    if(!pv) { return xs; }
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

void G4LivermoreGammaConversion5DModel::InitialiseForElement(
				      const G4ParticleDefinition*, 
				      G4int Z)
{
  G4AutoLock l(&LivermoreGammaConversion5DModelMutex);
  if(!data[Z]) { ReadData(Z); }
  l.unlock();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
