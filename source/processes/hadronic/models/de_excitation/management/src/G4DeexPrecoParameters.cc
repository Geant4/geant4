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
// $Id: G4DeexPrecoParameters.cc 68028 2013-03-13 13:48:15Z gcosmo $
//
// 15.03.2016 V.Ivanchenko 
//
// List of parameters of the pre-compound model
// and the deexcitation module
//

#include "G4DeexPrecoParameters.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4DeexParametersMessenger.hh"

#ifdef G4MULTITHREADED
G4Mutex G4DeexPrecoParameters::deexPrecoMutex = G4MUTEX_INITIALIZER;
#endif

G4DeexPrecoParameters::G4DeexPrecoParameters() 
{
  SetDefaults();
}

G4DeexPrecoParameters::~G4DeexPrecoParameters() 
{
  delete theMessenger;
}

void G4DeexPrecoParameters::SetDefaults()
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4DeexPrecoParameters::deexPrecoMutex);
#endif
  fStateManager = G4StateManager::GetStateManager();
  theMessenger = new G4DeexParametersMessenger(this);

  fLevelDensity = 0.10/CLHEP::MeV;
  fR0 = 1.5*CLHEP::fermi;
  fTransitionsR0 = 0.6*CLHEP::fermi;
  fFermiEnergy = 35.0*CLHEP::MeV; 
  fPrecoLowEnergy = 0.1*CLHEP::MeV; 
  fPhenoFactor = 1.0; 
  fMinExcitation = 10*CLHEP::eV;
  fMaxLifeTime = 1000*CLHEP::second;
  fMinExPerNucleounForMF = 100*CLHEP::GeV;
  fMinZForPreco = 3;
  fMinAForPreco = 5;
  fPrecoType = 3;
  fDeexType = 3;
  fTwoJMAX = 10;
  fNeverGoBack = false;
  fUseSoftCutoff = false;
  fUseCEM = true;
  fUseGNASH = false;
  fUseHETC = false;
  fUseAngularGen = true;
  fPrecoDummy = false;
  fCorrelatedGamma = false;
  fStoreAllLevels = false;
  fInternalConversion = true;
  fDeexChannelType = fEvaporation;
  fInternalConversionID = 
    G4PhysicsModelCatalog::Register("e-InternalConvertion");
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4DeexPrecoParameters::deexPrecoMutex);
#endif
}

void G4DeexPrecoParameters::SetLevelDensity(G4double val)
{
  if(IsLocked()) { return; }
  fLevelDensity = val/CLHEP::MeV;
}

void G4DeexPrecoParameters::SetR0(G4double val)
{
  if(IsLocked()) { return; }
  fR0 = val;
}

void G4DeexPrecoParameters::SetTransitionsR0(G4double val)
{
  if(IsLocked()) { return; }
  fTransitionsR0 = val;
}

void G4DeexPrecoParameters::SetFermiEnergy(G4double val)
{
  if(IsLocked()) { return; }
  fFermiEnergy = val;
}

void G4DeexPrecoParameters::SetPrecoLowEnergy(G4double val)
{
  if(IsLocked()) { return; }
  fPrecoLowEnergy = val;
}

void G4DeexPrecoParameters::SetPhenoFactor(G4double val)
{
  if(IsLocked()) { return; }
  fPhenoFactor = val;
}

void G4DeexPrecoParameters::SetMinExcitation(G4double val)
{
  if(IsLocked()) { return; }
  fMinExcitation = val;
}

void G4DeexPrecoParameters::SetMaxLifeTime(G4double val)
{
  if(IsLocked()) { return; }
  fMaxLifeTime = val;
}

void G4DeexPrecoParameters::SetMinExPerNucleounForMF(G4double val)
{
  if(IsLocked()) { return; }
  fMinExPerNucleounForMF = val;
}

void G4DeexPrecoParameters::SetMinZForPreco(G4int n)
{
  if(IsLocked() && n < 2) { return; }
  fMinZForPreco = n;
}

void G4DeexPrecoParameters::SetMinAForPreco(G4int n)
{
  if(IsLocked() && n < 0) { return; }
  fMinAForPreco = n;
}

void G4DeexPrecoParameters::SetPrecoModelType(G4int n)
{
  if(IsLocked() && n < 0) { return; }
  fPrecoType = n;
}

void G4DeexPrecoParameters::SetDeexModelType(G4int n)
{
  if(IsLocked() && n < 0) { return; }
  fDeexType = n;
}

void G4DeexPrecoParameters::SetTwoJMAX(G4int n)
{
  if(IsLocked() && n < 0) { return; }
  fTwoJMAX = n;
}

void G4DeexPrecoParameters::SetNeverGoBack(G4bool val)
{
  if(IsLocked()) { return; }
  fNeverGoBack = val;
}

void G4DeexPrecoParameters::SetUseSoftCutoff(G4bool val)
{
  if(IsLocked()) { return; }
  fUseSoftCutoff = val;
}

void G4DeexPrecoParameters::SetUseCEM(G4bool val)
{
  if(IsLocked()) { return; }
  fUseCEM = val;
}

void G4DeexPrecoParameters::SetUseGNASH(G4bool val)
{
  if(IsLocked()) { return; }
  fUseGNASH = val;
}

void G4DeexPrecoParameters::SetUseHETC(G4bool val)
{
  if(IsLocked()) { return; }
  fUseHETC = val;
}

void G4DeexPrecoParameters::SetUseAngularGen(G4bool val)
{
  if(IsLocked()) { return; }
  fUseAngularGen = val;
}

void G4DeexPrecoParameters::SetPrecoDummy(G4bool val)
{
  if(IsLocked()) { return; }
  fPrecoDummy = val;
  fDeexChannelType = fDummy;  
}

void G4DeexPrecoParameters::SetCorrelatedGamma(G4bool val)
{
  if(IsLocked()) { return; }
  fCorrelatedGamma = val; 
}

void G4DeexPrecoParameters::SetStoreICLevelData(G4bool val)
{
  if(IsLocked()) { return; }
  fStoreAllLevels = val;
}

void G4DeexPrecoParameters::SetStoreAllLevels(G4bool val)
{
  SetStoreICLevelData(val);
}

void G4DeexPrecoParameters::SetInternalConversionFlag(G4bool val)
{
  if(IsLocked()) { return; }
  fInternalConversion = val;
}

void G4DeexPrecoParameters::SetDeexChannelsType(G4DeexChannelType val)
{
  if(IsLocked()) { return; }
  fDeexChannelType = val;
}

std::ostream& G4DeexPrecoParameters::StreamInfo(std::ostream& os) const
{
  static const G4String namm[4] = {"Evaporation","GEM","Evaporation+GEM","Dummy"};
  static const G4int nmm[4] = {8, 68, 68, 0};
  size_t idx = (size_t)fDeexChannelType;

  G4int prec = os.precision(5);
  os << "=======================================================================" << "\n";
  os << "======       Pre-compound/De-excitation Physics Parameters     ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Type of pre-compound inverse x-section              " << fPrecoType << "\n";
  os << "Pre-compound model active                           " << (!fPrecoDummy) << "\n";
  os << "Pre-compound low energy (MeV)                       " 
     << fPrecoLowEnergy/CLHEP::MeV << "\n";
  os << "Type of de-excitation inverse x-section             " << fDeexType << "\n";
  os << "Type of de-excitation factory                       " << namm[idx] << "\n";
  os << "Number of de-excitation channels                    " << nmm[idx] << "\n";
  os << "Min excitation energy (keV)                         " 
     << fMinExcitation/CLHEP::keV << "\n";
  os << "Min energy per nucleon for multifragmentation (MeV) " 
     << fMinExPerNucleounForMF/CLHEP::MeV << "\n";
  os << "Level density (1/MeV)                               " 
     << fLevelDensity*CLHEP::MeV << "\n";
  os << "Time limit for long lived isomeres (ns)             " 
     << fMaxLifeTime/CLHEP::ns << "\n";
  os << "Internal e- conversion flag                         " 
     << fInternalConversion << "\n";
  os << "Store e- internal conversion data                   " << fStoreAllLevels << "\n";
  os << "Electron internal conversion ID                     " 
     << fInternalConversionID << "\n";
  os << "Correlated gamma emission flag                      " << fCorrelatedGamma << "\n";
  os << "Max 2J for sampling of angular correlations         " << fTwoJMAX << "\n";
  os << "=======================================================================" << "\n";
  os.precision(prec);
  return os;
}

void G4DeexPrecoParameters::Dump() const
{
  if (G4Threading::IsMasterThread()) { StreamInfo(G4cout); }
}

std::ostream& operator<< (std::ostream& os, const G4DeexPrecoParameters& par)
{
  return par.StreamInfo(os);
}

G4bool G4DeexPrecoParameters::IsLocked() const
{
  return (!G4Threading::IsMasterThread() ||
	  (fStateManager->GetCurrentState() != G4State_PreInit));
}
