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

#ifdef G4MULTITHREADED
G4Mutex G4DeexPrecoParameters::deexPrecoMutex = G4MUTEX_INITIALIZER;
#endif

G4DeexPrecoParameters::G4DeexPrecoParameters() 
{
  fStateManager = G4StateManager::GetStateManager();
  SetDefaults();
}

void G4DeexPrecoParameters::SetDefaults()
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4DeexPrecoParameters::deexPrecoMutex);
#endif
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
  fNeverGoBack = false;
  fUseSoftCutoff = false;
  fUseCEM = true;
  fUseGNASH = false;
  fUseHETC = false;
  fUseAngularGen = false;
  fUseLongFiles = true;
  fCorrelatedGamma = false;
  fStoreAllLevels = false;
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
  if(IsLocked()) { return; }
  fMinZForPreco = n;
}

void G4DeexPrecoParameters::SetMinAForPreco(G4int n)
{
  if(IsLocked()) { return; }
  fMinAForPreco = n;
}

void G4DeexPrecoParameters::SetPrecoModelType(G4int n)
{
  if(IsLocked()) { return; }
  fPrecoType = n;
}

void G4DeexPrecoParameters::SetDeexModelType(G4int n)
{
  if(IsLocked()) { return; }
  fDeexType = n;
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

void G4DeexPrecoParameters::SetUseFilesNEW(G4bool)
{}

void G4DeexPrecoParameters::SetCorrelatedGamma(G4bool val)
{
  if(IsLocked()) { return; }
  fCorrelatedGamma = val; 
}

void G4DeexPrecoParameters::SetStoreAllLevels(G4bool val)
{
  if(IsLocked()) { return; }
  fStoreAllLevels = val;
}

void G4DeexPrecoParameters::SetDeexChannelsType(G4DeexChannelType val)
{
  if(IsLocked()) { return; }
  fDeexChannelType = val;
}

std::ostream& G4DeexPrecoParameters::StreamInfo(std::ostream& os) const
{
  G4int prec = os.precision(5);
  os << "=======================================================================" << "\n";
  os << "======       Pre-compound/De-excitation Physics Parameters     ========" << "\n";
  os << "=======================================================================" << "\n";
  os << "Type of pre-compound inverse x-section              " << fPrecoType << "\n";
  os << "Type of de-excitation inverse x-section             " << fDeexType << "\n";
  os << "Min excitation energy (keV)                         " 
     << fMinExcitation/CLHEP::keV << "\n";
  os << "Level density (1/MeV)                               " 
     << fLevelDensity*CLHEP::MeV << "\n";
  os << "Time limit for long lived isomeres (ns)             " 
     << fMaxLifeTime/CLHEP::ns << "\n";
  os << "Use new data files                                  " << fUseLongFiles << "\n";
  os << "Use complete data files                             " << fStoreAllLevels << "\n";
  os << "Correlated gamma emission flag                      " << fCorrelatedGamma << "\n";
  os << "Electron internal conversion ID                     " << fInternalConversionID << "\n";
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
